import logging
import os
import re
import subprocess
import copy
from ete3 import Tree
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import sys
sys.path.append(os.path.abspath('./')) ## assumes running from repo root
import argparse

#UTILITY
def reformat(input, out_filename):
    """
    :Description: This is the function to transform the aligned fasta file from multiple line version to single line version, which is usually the input format of mega

    :param input: read in a fasta file
    :param out_filename: a string, which is name of the output file
    :return: NA. the output is saved as out_filename
    """

    entries = input.split('>')[1:]
    out = open(out_filename, 'w')

    for l in entries:
        id = re.sub(r'(\:|\'|\(|\)|\,|\?|\;)', '', l.split('\n')[0])
        seq = ''.join(l.split('\n')[1:])
        out.write('>%s\n%s\n' % (id, seq))

    out.close()


###REFACTORED PART BELOW
def align_ks_domains(reference_alignment, ks_names, ks_seq):
    """Function that aligns a number of query KS domain sequences to the
    reference alignment of KS domains.
    """
    #Set file names and write query domains to temp input file
    if not os.path.isdir(os.path.join(os.getcwd(), "processing")):
        os.mkdir(os.path.join(os.getcwd(), "processing"))
    in_temp = os.path.join(os.getcwd(), "processing", "in_seq" + ks_names[0] + ".fasta")
    in_temp_aligned = os.path.join(os.getcwd(), "processing", "in_seq_aligned" + ks_names[0] + ".fasta")
    out_temp = os.path.join(os.getcwd(), "processing", "out_seq" + ks_names[0] + ".fasta")
    alignment_file = os.path.join(os.getcwd(), "processing", "aligned" + ks_names[0] + ".fasta")
    with open(in_temp, "w") as tmp_input:
        for name, seq in zip(ks_names, ks_seq):
            tmp_input.write("%s\n%s\n" % (name, seq))

    #Generate alignment of query sequences
    muscle_cmd = str(MuscleCommandline(input=in_temp, out=in_temp_aligned, quiet='True'))
    sys.stdout.write("\tStarting muscle alignment of KSs...")
    result = subprocess.run(muscle_cmd.split(" "))
    sys.stdout.write("done\n")
    if result.returncode == 1 or 'error' in str(result.stderr).lower():
        logging.error("Alignment of query KS sequences with Muscle failed. Check if Muscle is installed appropriately.")
        sys.exit(1)

    #Align the query alignment to the reference alignment using muscle --profile
    muscle_cmd = str(MuscleCommandline(profile='True', in1=reference_alignment, in2=in_temp_aligned, out=out_temp, quiet='True'))
    sys.stdout.write("\tStarting muscle alignment of KSs to profile...")
    result = subprocess.run(muscle_cmd.split(" "))
    sys.stdout.write("done\n")
    if result.returncode == 1 or 'error' in str(result.stderr).lower():
        logging.error("Alignment of query+reference KS sequences with Muscle failed. Check if Muscle is installed appropriately.")
        sys.exit(1)
    else:
        f_temp_input = open(out_temp, 'r').read()
        reformat(input=f_temp_input, out_filename=alignment_file)

    #Remove temporary files
    for f in [in_temp, out_temp, in_temp_aligned]:
        os.remove(f)

    assert os.path.isfile(alignment_file)
    return alignment_file


def run_pplacer(reference_alignment, alignment_file, data_dir, ks_names):
    """Function that uses the reference tree with the new alignment to place
    query domains onto reference tree.
    """
    #Locations of files
    reference_tree = os.path.join(data_dir, "reference_tree.tre")
    pplacer_json = os.path.join(os.getcwd(), "processing", "pplacer_tree" + ks_names[0] + ".jplace")

    pplacer_cmd = ["pplacer", "-t", reference_tree, "-r", reference_alignment, "-o", pplacer_json, "-c", os.path.join(data_dir, "pplacer_reference"), alignment_file]
    sys.stdout.write("\tRunning pplacer...")
    result = subprocess.run(pplacer_cmd)
    sys.stdout.write("done\n")
    if result.returncode == 1:
        logging.error("Running pplacer failed. Check if the program is installed appropriately.")
        sys.exit(1)
    guppy_cmd = ["guppy", "sing", "--out-dir", os.path.join(os.getcwd(), "processing"), pplacer_json]
    sys.stdout.write("\tRunning guppy...")
    result = subprocess.run(guppy_cmd)
    sys.stdout.write("done\n")
    if result.returncode == 1:
        logging.error("Running guppy failed. Check if the program is installed appropriately.")
        sys.exit(1)
    assert os.path.isfile(os.path.join(os.getcwd(),"processing", "pplacer_tree" + ks_names[0] + ".sing.tre"))
    return os.path.join(os.getcwd(), "processing", "pplacer_tree" + ks_names[0] + ".sing.tre")

def get_leaf2clade(leaf2cladetbl):
    leaf2clade = {}
    clade2ann = {}
    with open(leaf2cladetbl) as c:
        for ln in c.read().splitlines():
            if ln[0] is '#':
                # if comment
                continue
            ksname, clade, ann = ln.split("\t")
            leaf2clade[ksname] = clade
            clade2ann[clade] = ann
    return leaf2clade, clade2ann

def get_clade(q, tree_hits, t, funClades):
    cladelist = {}
    qname = tree_hits[q].name
    newtree = copy.deepcopy(t)
    pref = re.sub(r"^(.+)_#\d+_M=\d+?\.?\d*$", "\g<1>", qname)
    keep = []
    lcount = 0
    for leaf in newtree.get_leaves():
        lcount += 1
        ln = leaf.name.split("_")
        if re.match("^#\d+$", ln[-2]) is None:
            keep.append(leaf)
        elif ln[-2] == '#'+q:
            keep.append(leaf)
    newtree.prune(keep)
    grandparent, parent, elder = None, None, None
    for leaf in newtree.get_leaves():
        if leaf.name == qname:
            parent = leaf.up
            grandparent = leaf.up.up
            break
    #print "\nSiblings: "+str(len(parent.get_leaves()))
    #print "Cousins: "+str(len(grandparent.get_leaves()))
    if len(parent.get_leaves()) > 2:
        elder = parent
    else:
        elder = grandparent
    for leaf in elder.get_leaves():
        if leaf.name == qname:
            continue
        elif leaf.name in funClades:
            if funClades[leaf.name] == 'query_seq':
                continue
            elif funClades[leaf.name] in cladelist.keys():
                cladelist[funClades[leaf.name]] += 1
            else:
                cladelist[funClades[leaf.name]] = 1
        else:
            if 'unknown_clade' in cladelist.keys():
                cladelist['unknown_clade'] += 1
            else:
                cladelist['unknown_clade'] = 1
    clade_assignment = 'clade_not_conserved'
    if len(cladelist.keys()) == 1:
        for k in cladelist:
            clade_assignment = k
    return pref, clade_assignment

def parse_pplacer(pplacer_tree, data_dir, masscutoff, mode, funClades):
    monoclade = {}
    with open(pplacer_tree) as f:
        for ln in f.read().splitlines():
            t = Tree(ln)
            tree_hits = {}
            ## Identify the list of pplacer placements
            for leaf in t.get_leaves():
                ln = leaf.name.split("_")
                if re.match("^#\d+$", ln[-2]) is not None:
                    n = re.sub(r"^#(\d+)$", "\g<1>", ln[-2])
                    tree_hits[n] = leaf
                    funClades[leaf.name] = 'query_seq'
            ## Look to see when threshold is met
            if mode == 'top_down':
                to_check = []
                running_mass = 0
                for placement in sorted(tree_hits): ## NOTE: unsure about behavior when n>=10
                    running_mass += float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                    to_check.append(placement)
                    if(running_mass >= masscutoff):
                        break
                cl = {}
                for q in to_check:
                    pref, clade_assignment = get_clade(q, tree_hits, t, funClades)
                    if clade_assignment in cl.keys():
                        cl[clade_assignment] += 1
                    else:
                        cl[clade_assignment] = 1;
                if len(cl.keys()) == 1:
                    for k in cl:
                        monoclade[pref] = k
                else:
                    monoclade[pref] = 'clade_not_conserved'
            elif mode == 'additive':
                totalmass = {}
                pref = ''
                for placement in tree_hits:
                    mass = float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                    pref, clade_assignment = get_clade(placement, tree_hits, t, funClades)
                    if clade_assignment in totalmass:
                        totalmass[clade_assignment] += mass
                    else:
                        totalmass[clade_assignment] = mass
                best_clade = max(totalmass, key=totalmass.get)
                if totalmass[best_clade] >= masscutoff:
                    monoclade[pref] = best_clade
                else:
                    monoclade[pref] = 'clade_not_conserved'
            else:
                print(mode+'  UNRECOGNIZED MODE...(will be an error in future)')
    return monoclade

def write2fasta(ks_names, ks_seqs, outfile):
    f = open(outfile, "w")
    for i in range(0,len(ks_names),1):
        f.write(ks_names[i]+"\n"+ks_seqs[i]+"\n")
    f.close()
    return outfile


def run_pipeline_pplacer(reference_alignment, alignment_file, data_dir, masscutoff, leaf2clade, ks_names):
    logging.info("Phylogenetic analysis: predicting substrate specificity of KS")
    pplacer_tree = run_pplacer(reference_alignment, alignment_file, data_dir, ks_names)
    callmode = 'additive'
    sys.stdout.write("\tParsing ppclacer results...")
    clade_assignment = parse_pplacer(pplacer_tree, data_dir, masscutoff, callmode, leaf2clade)
    sys.stdout.write("done\n")

    return clade_assignment


def main():
    parser = argparse.ArgumentParser(description='Generate transPACT dendrograms.')
    parser.add_argument('-data_dir', help='The directory containing the reference data files.', default=os.path.join(os.getcwd(), 'transPACT_data'))
    parser.add_argument('-reference_alignment', help='Path to the reference alignment file.', default = os.path.join(os.getcwd(),'transPACT_data','reference_alignment.fasta'))
    parser.add_argument('-clade_annotation', help='Path to clade annotation file.', default = os.path.join(os.getcwd(), 'transPACT_data', 'clade_annotation.txt'))
    parser.add_argument('-output', help='The output directory',default=os.path.join(os.getcwd(),'output'), type=str)
    parser.add_argument('-input', help='Path to input fasta file.')
    args = parser.parse_args()

    reference_alignment = args.reference_alignment
    data_dir = args.data_dir
    leaf2cladetbl = args.clade_annotation
    leaf2clade, clade2ann = get_leaf2clade(leaf2cladetbl)
    if os.path.isfile(args.input):
        fasta = SeqIO.parse(open(args.input),'fasta')
        output_file_name = os.path.splitext(args.input)[0]
    else:
        print('Did not detect a fasta file. Exiting...')
        return

    ks_names = []
    ks_seqs = []
    for rec in fasta:
        ks_names.append('>'+rec.id) ## NOTE: this assumes this header will be unique when compared to the reference set. if not, will return OUTGROUP regardless of the call. If running such non-uniue seqs, change the header to '>query_'+rec.id or similar make it unique
        ks_seqs.append(rec.seq)
    sys.stdout.write("Starting aligning KS sequences...\n")
    alignment_file = align_ks_domains(reference_alignment, ks_names, ks_seqs)
    sys.stdout.write("done\nStarting clade annotation...\n")
    # Individualize output_file_name
    clade = run_pipeline_pplacer(reference_alignment, alignment_file, data_dir, 0.6, leaf2clade, ks_names)
    sys.stdout.write("done\n")
    #Remove temporary files

    for f in ['pplacer_tree' + ks_names[0] + '.jplace', 'pplacer_tree' + ks_names[0] + '.sing.tre', 'precomp_funcAnnotate_perKS.txt']:
        try:
            os.remove(os.path.join(os.getcwd(), "processing", f))
        except OSError as error:
            pass

    with open(''.join((output_file_name, '_funcAnnotate_perKS.txt')), 'w') as out:
        for ks in clade:
            if clade[ks] in clade2ann:
                print("\t".join([str(ks), str(clade[ks]), clade2ann[clade[ks]]]))
                out.write(f"{ks}: {clade[ks]}|{clade2ann[clade[ks]]}\n")
            else:
                print("\t".join([str(ks), str(clade[ks]), 'NA']))
                out.write(f"{ks}: NA\n")
    return clade


if __name__ == '__main__':
    main()
