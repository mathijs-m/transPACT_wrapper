# -*- coding: utf-8 -*-
"""
A script to rerun the entire transPACT pipeline automatically on a folder of
trans-AT PKS genbanks and produce the trees as published in the original paper
"""
from Bio import SeqIO
from collections import defaultdict
import argparse
import multiprocessing as mp
from pathlib import Path
import subprocess
import sys
import re
import os
import pickle as pkl
try:
    from ete3 import Tree, RectFace, TextFace, NodeStyle
    from ete3.treeview import TreeStyle
except ImportError:
    print('Error loading the ete3 modules')


def find_asdomains(filepath, domain_type=None):
    """
    Finds aSDomains and groups them by CDS.
    :param filepath: String with path to the genbank file
    :param domain_type: The type of aSDomain to look for. If None, all aSDomains are returned
    :return: A dictionary of CDSs and their aSDomains
    """
    records = SeqIO.parse(filepath, 'genbank')

    grouped_asdomains = defaultdict(list)
    for record in records:
        if domain_type is None:
            asdomains = [feature for feature in record.features if feature.type == "aSDomain"]
        else:
            asdomains = [feature for feature in record.features if
                         feature.type == "aSDomain" and feature.qualifiers['aSDomain'][0] == domain_type]
        cdss = [feature for feature in record.features if feature.type == "CDS"]
        pairs = [(cds, asdomain) for asdomain in asdomains for cds in cdss if asdomain.location.start in cds]
        for k, v in pairs:
            grouped_asdomains[k].append(v)

    return grouped_asdomains


def find_NRPS_asdomains(filepath, domain_type=None):
    """
    Finds aSDomains and groups them by CDS.
    The direction of the PKSs is determined using the weighted direction of all
    target domain-containing CDSs
    
    :param filepath: String with path to the genbank file
    :param domain_type: The type of aSDomain to look for. If None, all aSDomains are returned
    :return: A dictionary of CDSs and their aSDomains
    """
        
    records = SeqIO.parse(filepath, 'genbank')

    grouped_asdomains = defaultdict(list)
    for record in records:
        if domain_type is None:
            asdomains = [feature for feature in record.features if feature.type == "aSDomain"]
        else:
            asdomains = [feature for feature in record.features if
                         feature.type == "aSDomain" and feature.qualifiers['aSDomain'][0] == domain_type]
        cdss = [feature for feature in record.features if feature.type == "CDS"]
        asmodules = [feature for feature in record.features if feature.type == "aSModule"]
        pairs = [(cds, asmodule, asdomain) for asdomain in asdomains for asmodule in asmodules for cds in cdss
                 if asdomain.location.start in cds and asdomain.location.start in asmodule]

        direction = -1 if sum(
            [cds.location.strand * abs(cds.location.start - cds.location.end) for cds in cdss]) < 0 else 1

        if direction == 1:
            for cds, module, domain in pairs:
                grouped_asdomains[cds].append(domain)
        elif direction == -1:
            for cds, module, domain in reversed(pairs):
                grouped_asdomains[cds].append(domain)

    return grouped_asdomains


def find_transAT_asdomains(filepath, domain_type=None):
    """
    Finds aSDomains and groups them by CDS.
    The direction of the PKSs is determined using the weighted direction of all
    KS-containing CDSs
    
    :param filepath: String with path to the genbank file
    :param domain_type: The type of aSDomain to look for. If None, all aSDomains are returned
    :return: A dictionary of CDSs and their aSDomains
    """
    records = SeqIO.parse(filepath, 'genbank')

    grouped_asdomains = defaultdict(list)
    for record in records:
        if domain_type is None:
            asdomains = [feature for feature in record.features if feature.type == "aSDomain"]
        else:
            asdomains = [feature for feature in record.features if
                         feature.type == "aSDomain" and feature.qualifiers['aSDomain'][0] == domain_type]
        ATdomains = [feature for feature in record.features if
                     feature.type == 'aSDomain' and feature.qualifiers['aSDomain'][0] == 'PKS_AT']
        cdss = [feature for feature in record.features if feature.type == "CDS"]
        asmodules = [feature for feature in record.features if feature.type == "aSModule"]
        pairs = [(cds, asmodule, asdomain) for asdomain in asdomains for asmodule in asmodules for cds in cdss
                 if asdomain.location.start in cds and asdomain.location.start in asmodule and
                 not any([ATdomain.location.start in asmodule for ATdomain in ATdomains])]

        direction = -1 if sum(
            [cds.location.strand * abs(cds.location.start - cds.location.end) for cds in cdss]) < 0 else 1

        if direction == 1:
            for cds, module, domain in pairs:
                grouped_asdomains[cds].append(domain)
        elif direction == -1:
            for cds, module, domain in reversed(pairs):
                grouped_asdomains[cds].append(domain)

    return grouped_asdomains


def generate_color_dict(file='./data/transAT_mc.txt'):
    """
    Generate a color dictionary for the KS architecture visualization
    :param file: String with path to the file containing the color information 
    :return: A dictionary with the KS architecture as key and the color as value
    """
    
    txt = open('./data/transAT_mc.txt', 'r').readlines()
    # txt = [line.replace('amino_acids','amino acids').replace('aMe_bOH', 'a_Me_b_OH').replace('aMe_red/bOH/bketo', '')
    # for line in txt]

    color_dict = {line.split('\t')[1].replace('|', '/'): line.split('\t')[-1].replace('\n', '') for line in txt}
    return color_dict


def generate_ks_dict(clade_output):
    """
    Generate a dictionary with the BGC identifies as key and the KSS clade sequence as value
    :param clade_output: String with path to the file containing the KSS clade information
    :return: A dictionary with the BGC identifier as key and the KSS clade sequence as value
    """
    txt = open(clade_output, 'r').readlines()[1:]
    ks_dict = dict()
    for line in txt:
        cluster = line.split('\t')[0].replace('|c1', '')
        KS_clades = [clade.replace('\n', '').split('|')[0] for clade in line.split('\t')[1:]]
        ks_dict[cluster] = KS_clades
    return ks_dict


def generate_simplified_color_dict():
    """
    Generate a dictionary with the simplified clade names as key and the color as value. These values are used in
    the figure in the manuscript
    :return: A dictionary with the simplified clade names as key and the color as value
    """
    simplified_colors = open(r'./transPACT_data/Simplified_colors.txt', 'r').readlines()
    simplified_colors = {line.split('\t')[0]: line.split('\t')[1].replace('\n', '') for line in simplified_colors}

    simplified_clades = open(r'./transPACT_data/simplified_clades.txt', 'r').readlines()[1:]
    simplified_clades = {line.split('\t')[1]: line.split('\t')[-1].replace('\n', '') for line in simplified_clades}

    simplified_colors_dict = {clade: simplified_colors[simplified_clades[clade]] for clade in simplified_clades if
                              simplified_clades[clade] in simplified_colors}
    return simplified_colors_dict


def extract_organism_and_taxonomy_from_gb(file):
    """
    Extracts the organism and taxonomy from a genbank file
    :param file: String with path to the genbank file
    :return: The accession, organism and taxonomy 
    """ 
    records = SeqIO.parse(file, 'genbank')
    accession = file.name.split('.gbk')[0]
    organism = list()
    taxonomy = list()
    for i, record in enumerate(records):
        try:
            organism.append(record.annotations['organism'])
            taxonomy.append(record.annotations['taxonomy'])
        except KeyError:
            organism.append('')
            taxonomy.append('')
            print('No organism or taxonomy found in file %s' % file)
        if i > 1:
            print('More than one record found in file %s' % file)
    return accession, organism, taxonomy


def generate_accession_organism_dict(input_dirs, cpus=10):
    """
    Generate a dictionary with that contains the organism and taxonomy for each accession number that is used in the analysis
    :param input_dirs: List of directories that contain the genbank files
    :param cpus: Number of cpus to use
    :return: A dictionary with the accession number as key and the organism and taxonomy as value
    """
    gbk_files = []
    gbk_files.extend([file for input_dir in input_dirs for file in Path(input_dir).rglob('*.gbk')])

    accession_dict = defaultdict(dict)
    with mp.Pool(min(cpus, mp.cpu_count())) as p:
        results = p.map(extract_organism_and_taxonomy_from_gb, gbk_files)
    for accession, organism, taxonomy in results:
        if len(organism) > 1:
            print(f'More than one organism found for {accession}. Only using the first one: {organism[0]}. Best check the file!')
        accession_dict[accession]['organism'] = organism[0]
        accession_dict[accession]['taxonomy'] = taxonomy[0]

    return accession_dict


def generate_MIBiG_ref_dict():
    """
    Generate a dictionary with the MIBiG identifier as key and the compound produced by that cluster as value
    :return: A dictionary with the MIBiG identifier as key and the compound produced by that cluster as value
    """
    with open(r'./transPACT_data/MiBIG_references.txt',
              'r') as MiBIG_ref_file:
        MiBIG_ref_file = MiBIG_ref_file.readlines()
        MiBIG_ref_dict = {line.split('\t')[0]: line.split('\t')[1].replace('\n', '') for line in MiBIG_ref_file}
    return MiBIG_ref_dict


def add_name_face_to_tree(node, accession_dict, MIBiG_refs, font_size=10):
    """
    Generate an ete3 TextFace with the name of the organism that can be added to the node of the tree
    :param node: The ete3 node object for which the TextFace should be generated
    :param accession_dict: Dictionary with the accession number and organisms and taxonomies, generated by
                           generate_accession_organism_dict
    :param MIBiG_refs: Dictionary with MIBiG references, generated by generate_MIBiG_ref_dict
    :param font_size: The font size of the text in the TextFace
    :return: An ete3 TextFace with the accession, name of the organism and, if applicable, the MIBiG reference
    """
    # Add the name of the organism to the tree
    if 'BGC0' in node.name:
        # Annotate leaves from MIBiG
        try:
            name_face = TextFace(
                f"{accession_dict[node.name]['organism']} - {node.name.split('.region')[0]} - {MIBiG_refs[node.name.split('.region')[0]].title()}",
                fsize=font_size, fgcolor='teal', ftype='Arial')
        except KeyError:
            print(f"MIBiG reference not found for {node.name}.\n")
            name_face = TextFace(f"{accession_dict[node.name]['organism']} - {node.name}", fsize=font_size, ftype='Arial')
    elif 'region' in node.name:
        # Annotate leaves from antiSMASH
        name_face = TextFace(f"{accession_dict[node.name]['organism']} - {node.name}", fsize=font_size, ftype='Arial')
    else:
        print(f"Node name {node.name} does not contain 'BGC0' or 'region'.\n")
        name_face = TextFace(node.name, fsize=font_size)
    return name_face


def label_Acidobacteria(node, name_face, accession_dict, font_size, block_height):
    """
    Label the nodes of the tree that belong to the phylum Acidobacteria
    Label the nodes belonging to the pho and phb clusters specifically
    :param node: ete3 node object
    :param name_face: ete3 TextFace with the name of the organism
    :param accession_dict: Dictionary with the accession number and organisms and taxonomies
    :param font_size: Float with the font size of the text in the TextFace
    :param block_height: Float with the height of the block in the tree
    :return: the ete3 node and TextFace with the name of the organism
    """
    if 'CP071793.1.region015' in node.name or 'JAFREP010000003.1.region004' in node.name:
        if 'BGC0' in node.name:
            name_face = TextFace(f"{accession_dict[node.name]['organism']} - {node.name}",
                                 fsize=font_size, fgcolor='royalblue')
        else:
            name_face = TextFace(f"{accession_dict[node.name]['organism']} - {node.name}", fsize=font_size,
                                 fgcolor='royalblue', penwidth=2)
        node_style = NodeStyle()
        node_style['hz_line_color'] = 'royalblue'
        node_style['hz_line_width'] = 10
        node.set_style(node_style)
        #node.add_face(RectFace(width=block_height, height=500, fgcolor='royalblue', bgcolor='royalblue'),
        #              column=-1,
        #              position='aligned')
    elif 'Acidobacteria' in accession_dict[node.name]['taxonomy']:
        if 'BGC0' in node.name:
            name_face = TextFace(
                f"{accession_dict[node.name]['organism']} - {node.name}",
                fsize=font_size, fgcolor='Tomato')
        else:
            name_face = TextFace(f"{accession_dict[node.name]['organism']} - {node.name}", fsize=font_size,
                                 fgcolor='Tomato', penwidth=2)
        node_style = NodeStyle()
        node_style['hz_line_color'] = 'Tomato'
        node_style['hz_line_width'] = 10
        node.set_style(node_style)
        #node.add_face(RectFace(width=block_height, height=500, fgcolor='Tomato', bgcolor='Tomato'),
        #              column=-1,
        #              position='aligned')

    return node, name_face


def add_radial_KS_distribution(node, ks_dict, color_dict, block_width, block_height):
    """
    Add the colored radial KS distribution to a leaf
    :param node: ete3 node object
    :param ks_dict: Dictionary with the KS sequences for each BGC
    :param color_dict: Dictionary with the colors for each KS type
    :param block_width: Float with the width of the block in the tree
    :param block_height: Float with the height of the block in the tree
    :return: ete3 node object
    """

    # Add the radial KS distribution to a leaf
    for i, ks in enumerate(ks_dict[node.name.split('|c')[0]]):
        if all([ks == 'NA' for ks in ks_dict[node.name.split('|c')[0]][i:]]):
            break
        ks = ks.split(' - ')[0]
        if ks == 'NA':
            node.add_face(RectFace(width=block_width, height=block_height, fgcolor='lightgray', bgcolor='white',
                                   label={'text': '', 'fontsize': '10', 'color': 'gray'}), column=1 + i,
                          position='aligned')
        elif ks in color_dict.keys():
            node.add_face(
                RectFace(width=block_width, height=block_height, fgcolor=color_dict[ks], bgcolor=color_dict[ks]),
                column=1 + i, position='aligned')
        else:
            sys.stdout.write(f"Error for {ks}\n")
            node.add_face(RectFace(width=block_width, height=block_height, fgcolor='black', bgcolor='black'),
                          column=1 + i, position='aligned')
    return node


def convert_leaves_to_ultrametric(tree, accession_dict):
    """
    Convert the branches that correspond to members of the Acidobaceriaceae to an ultrametric branches.
    This allows nice coloring of those branches all the way to the outside of the tree.
    :param tree: ete3 tree object
    :param accession_dict: Dictionary with the accession number and organisms and taxonomies
    :return: Altered ete3 tree object
    """
    most_distant_leaf, tree_length = tree.get_farthest_leaf()
    current_dist = 0
    for postorder, node in tree.iter_prepostorder():
        if postorder:
            current_dist -= node.dist
        else:
            if node.is_leaf():
                try:
                    if 'Acidobacteria' in accession_dict[node.name]['taxonomy']:
                        node.dist += tree_length - (current_dist + node.dist)
                except KeyError:
                    print(f"KeyError while setting branch lengths for {node.name}!\n")
            elif node.up:  # node is internal
                current_dist += node.dist
    return tree


def get_organism_and_taxonomy(file):
    """
    Get the organism and taxonomy from the genbank file
    :param file: genbank file
    :return: organism and taxonomy
    """
    gbs = SeqIO.parse(file, 'genbank')
    accession = file.split('.gbk')[0]
    accession_dict = defaultdict(dict)
    organisms = []
    taxonomies = []
    print_ = True
    for gb in gbs:
        organisms.append(gb.annotations['organism'])
        taxonomies.append(gb.annotations['taxonomy'])
    if len(set(organisms)) > 1 and print_:
        print(f"More than one organism found for {accession}!")
        print_ = False
    accession_dict[accession]['organism'] = organisms[0]
    accession_dict[accession]['taxonomy'] = taxonomies[0]
    return accession_dict


def make_treestyle(mode='c'):
    """
    Make the tree style
    :param mode: String with the mode of the tree, c for circular, r for radial
    :return: ete3 TreeStyle object
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    if mode == 'c':
        ts.mode = "c"
        ts.arc_start = -180
        ts.arc_span = 360
    ts.allow_face_overlap = True
    ts.force_topology = False
    ts.draw_guiding_lines = True
    ts.root_opening_factor = 0
    ts.guiding_lines_type = 0
    ts.guiding_lines_color = 'Gainsboro'
    ts.show_leaf_name = False
    ts.show_scale = False
    return ts


def process_node(node, accession_dict, MiBIG_refs, ks_dict, color_dict, font_size, block_width, block_height,
                 node_style, names=False):
    """
    Process the node and add the faces
    :param node: node to process
    :param accession_dict: dictionary with the accession information
    :param MiBIG_refs: dictionary with the MIBiG references
    :param ks_dict: dictionary with the KS distribution
    :param color_dict: dictionary with the colors for the KS

    :param font_size: font size for the labels
    :param block_width: width of the KS blocks
    :param block_height: height of the KS blocks
    :param node_style: style for the nodes
    :param names: boolean to add the names to the leaves

    :return: node with the faces
    """
    # Hide node circles
    node.img_style['size'] = 0
    node.set_style(node_style)
    if node.dist == 0.0:
        # Neighbouring nodes with distance 0 don't render properly. Hence, give them a negligible distance.
        node.dist = 1e-6

    # Iterate over the nodes and add the labels to the leaves
    if node.is_leaf():
        # Construct the leaf naming
        try:
            name_face = add_name_face_to_tree(node, accession_dict, MiBIG_refs, font_size=font_size)
            # Label all the Acidobacteria
            node, name_face = label_Acidobacteria(node, name_face, accession_dict, font_size, block_height)
        except KeyError:
            sys.stdout.write(f"KeyError for organism for {node.name}.\n")

        # Add the label to the leaf
        if names:
            node.add_face(name_face, column=0, position='aligned')

        # Generate the radial KS plots
        try:
            node = add_radial_KS_distribution(node, ks_dict, color_dict, block_width, block_height)
            pass
        except KeyError:
            print(f"KeyError for {node.name.split(' ')[0]}!\n")
            node.add_face(RectFace(width=block_height, height=block_height, fgcolor='gray', bgcolor='gray'), column=-1,
                          position='aligned')
    return node


def generate_tree(input_dirs, output_file, tree_output, clade_output, generate_new_accession_dict=True):
    """
    Generate a phylogenetic tree from the transPACT output
    :param input_dirs: List of directories with the Genbanks from which the transPACT output was generated
    :param output_file: String with the path to the transPACT output Newick file
    :param tree_output: String with the path for the output tree that is generated  by this function
    :param clade_output: String with the path to the transPACT output file, containing the KS clading matrix
    :param generate_new_accession_dict: Boolean to generate a new accession dictionary. Generating a new dictionary is
                                        necessary if the input directories have changed.
    :return: None
    """
    # Generate a phylogenetic tree
    accession_dict_path = r'./transPACT_data/accession_dict.pkl'
    if generate_new_accession_dict:
        print('Generating a new accession dictionary...')
        accession_dict = dict(generate_accession_organism_dict(input_dirs))
        pkl.dump(accession_dict, open(accession_dict_path, 'wb'))
        print(f'Done!\n\tSaved to {accession_dict_path}.\n')
    else:
        accession_dict = pkl.load(open(accession_dict_path, 'rb'))

    # Generate some dictionaries to look up data for the tree
    color_dict = generate_simplified_color_dict()
    MiBIG_refs = generate_MiBIG_ref_dict()
    ks_dict = generate_ks_dict(clade_output)

    # Generate the tree
    newick = open(output_file, 'r').read()
    tree = Tree(newick, format=1)
    tree.ladderize()
    block_height = 2000  # The scaling of the blocks and fonts is a little weird.
    block_width = 20000 # The block width is responsible for scaling further outwards.
    font_size = 700
    mode = 'c'

    # Some notes on the output file:
    # PNG/JPG are the most compact file formats and open significantly quicker than vector formats. However, for large
    # trees, zooming and checking the names might not be possible
    # PDFs allow zooming to see the names and individual cluster architectures. However. editting pdfs in vector
    # graphics software is very slow and seems to require a lot of memory. Opens relatively quickly in Adobe Acrobat,
    # however.
    # SVGs can be opened by fewer programs, but Illustrator can readily open SVGs of even very large trees to edit them.

    # Some style settings for the nodes
    node_style = NodeStyle()
    node_style['vt_line_width'] = 0.1
    node_style['hz_line_width'] = 0.1
    node_style['size'] = 0

    for node in tree.traverse():
        node = process_node(node, accession_dict, MiBIG_refs, ks_dict, color_dict, font_size, block_width,
                            block_height, node_style, names=False)

    # Convert to ultrametric
    tree = convert_leaves_to_ultrametric(tree, accession_dict)

    # Render the tree
    _ = tree.render(str(tree_output), tree_style=make_treestyle(mode), w=600, h=600, units='mm', dpi=300)


def run_pplacer_on_single_KS(KS_name, KS_seq, processing_dir):
    """
    Run the pplacer algorithm on a single, omitted, KS sequence
    :param KS_name: String with the name of the KS sequence to be queried
    :param KS_seq: String with the KS amino acid sequence of the KS to be queried
    :param processing_dir: String with the path to the directory where the temporary files are stored
    :return: List with the KS names and the corresponding substrate names
    """
    # Run the pplacer algorithm on a single, omitted, KS sequence
    filename = processing_dir.joinpath(f"{KS_name.replace('|', '_')}.fasta")
    with open(filename, 'w') as file:
        file.write(f">{KS_name}\n{KS_seq}")
    subprocess.Popen(f"python substrate_from_faa.py -input {filename}", shell=True,
                     stdout=subprocess.PIPE).communicate()  # Change the script to accept string
    txt = open(os.path.splitext(filename)[0] + '_funcAnnotate_perKS.txt', 'r').read()
    os.remove(filename)
    os.remove(os.path.splitext(filename)[0] + '_funcAnnotate_perKS.txt')
    return txt.split(': ')[0], txt.split(': ')[1].replace('\n', '')


def save_fastas(aS_domains, processing_dir, output_dir, block_size=100):

    """
    # Function to save the sequences in fasta files of with block_size sequences per file
    # Breaking the sequences up in blocks allows for parallelization of the alignment, which
    # significantly increases speed.

    :param aS_domains: List of tuples with the file name and the CDSs of the file
    :param processing_dir: String with the path to the directory where the temporary files are stored
    :param output_dir: String with the path to the directory where the fasta files are stored
    :param block_size: Integer with the number of sequences per fasta file
    :return: List with the paths to the fasta files and a string with the path to the fasta file of all KS sequences
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(processing_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout.write('Starting parsing ')

    sequences = list()
    for file, cdss in aS_domains:
        counter = 1
        for cds in cdss:
            for aS_domain in cdss[cds]:
                sequences.append((f">{file.stem}_{aS_domain.location.start.real}|KS{counter}",
                                  aS_domain.qualifiers['translation'][0]))
                counter += 1
    sys.stdout.write(f"{len(sequences)} sequences to fastas...")

    block_paths = list()
    for i, (name, sequence) in enumerate(sequences):
        if i % block_size == 0:
            block_path = processing_dir.joinpath(f"KS_sequences{i}-{min(i + block_size, len(sequences))}.fasta")
            block_paths.append(block_path)
            fasta_file = open(block_path, 'w')
        fasta_file.write(f"{name}\n{sequence}\n")
    fasta_file.close()

    # Also write all sequences to a single file
    with open(output_dir.joinpath('All_KS_sequences.fasta'), 'w') as fasta_file:
        for name, sequence in sequences:
            fasta_file.write(f"{name}\n{sequence}\n")

    sys.stdout.write('done\n')

    return block_paths, output_dir.joinpath('All_KS_sequences.fasta')


def parse_and_add_missing_KSs(processing_dir, threads):
    """
    Guppy appears to miss some KSs in the output. Check for missing ones and
    Run transPACT on these individual ones.

    :param processing_dir: pathlib Path object with the path to the directory where the temporary files are stored
    :param threads: Integer with the number of threads to use
    :return: List with the names of the missing KSs and list with the names of the KSs that were not parsed properly
    """

    # Get al KS sequences
    KSs = dict()
    files = [file for file in processing_dir.glob('./*-*.fasta')]
    for file in files:
        with open(file, 'r') as file:
            fasta = file.read()
        KSs.update({seq.split('\n')[0]: seq.split('\n')[1] for seq in fasta.split('>') if len(seq) > 2})

    # Check al properly annotated KS sequences
    parsed_KSs = list()
    missing_KSs = list()
    txts = [file for file in processing_dir.glob('./*_funcAnnotate_perKS.txt')]
    for txt in txts:
        txt = open(txt, 'r').readlines()
        for line in txt:
            if line.count('|') == 1:
                parsed_KSs.append((line.split(': ')[0], line.split(': ')[1].replace('\n', '')))
            else:
                *unparsed_KSs, parsed_KS = re.split('(\|KS[0123456789]+_)', line, maxsplit=line.count('|') - 1)
                parsed_KSs.append((parsed_KS.split(': ')[0], parsed_KS.split(': ')[1].replace('\n', '')))
                unparsed_KSs = [''.join(unparsed_KSs[i:i + 2]) for i in range(0, len(unparsed_KSs), 2)]
                unparsed_KSs = [KS[0:-1] if KS[-1] == '_' else KS for KS in unparsed_KSs]
                for unparsed_KS in unparsed_KSs:
                    missing_KSs.append(unparsed_KS)
    pool = mp.Pool(threads)
    parsed_KSs.extend(
        pool.starmap(run_pplacer_on_single_KS, [(KS_name, KSs[KS_name], processing_dir) for KS_name in missing_KSs]))
    return parsed_KSs, missing_KSs


def remove_processing_files(processing_dir):
    """
    Remove all files from the processing directory
    :param processing_dir: pathlib Path object with the path to the directory where the temporary files are stored
    :return: None
    """
    files = [file for file in processing_dir.iterdir() if any(['.txt' in file, 'fasta' in file])]
    for file in files:
        os.remove(file)


def substrate_from_faa_submitter(ks_fasta_path, processing_dir):
    """
    Simple function to submit transPACT in parallel
    :param ks_fasta_path: String with the path to the fasta file with the KS sequences
    :param processing_dir: String with the path to the directory where the temporary files are stored
    :return:
    """
    subprocess.Popen(f"python substrate_from_faa.py -input {ks_fasta_path} -output {processing_dir}", shell=True,
                     stdout=subprocess.PIPE).communicate()


def write_KS_annotation_matrix_file(parsed_KSs, clade_output):
    """
    Write a txt file with the matrix with the KS annotation for each clade
    :param parsed_KSs: List with tuples of KS names and their annotation
    :param clade_output: String with the path to the output file
    :return: None
    """
    KSs_clustered = defaultdict(list)
    for KS, clade in parsed_KSs:
        KSs_clustered['_'.join(KS.split('|')[0].split('_')[0:-1])].append((KS.split('|')[1], clade))

    max_length = max([len(KSs_clustered[clade]) for clade in KSs_clustered])
    with open(clade_output, 'w') as file:
        file.write("\t" + '\t'.join([f"KS{i + 1}" for i in range(max_length)]) + '\n')
        for cluster in KSs_clustered:
            try:
                KSs_clustered[cluster] = sorted(KSs_clustered[cluster], key=lambda x: int(x[0].replace('KS', '')))
            except:
                sys.stdout.write(f"Could not sort {cluster}\n")
                continue
            # Something in the output is going wrong here with e.g. NZ_KN323025.1.region001 KS3 not having a clade assigned
            output_string = ['NA' for _ in range(max_length)]
            output_string[0:len(KSs_clustered[cluster])] = [ks_data[1] for ks_data in KSs_clustered[cluster]]
            output_string = '\t'.join(output_string)
            file.write(f"{cluster}\t{output_string}\n")
    sys.stdout.write(f"\tSaved the clading output to {clade_output}.\n")


def run_hmmalign_on_all_KSs(processing_dir, data_dir):
    """
    Run hmmalign on all KS sequences
    :param processing_dir: String with the path to the directory where the temporary files are stored
    :param data_dir: String with the path to the directory where the hmm file is stored
    :return: String with the path to the hmmalign output file and the hmmalign stdout message
    """
    all_ks_fasta = Path(processing_dir).joinpath('All_KS_sequences.fasta')
    hmm_file = Path(data_dir).joinpath('PKS_KS.hmm')
    hmmalign_output = Path(processing_dir).joinpath('All_KS_sequences.hmmalign')
    hmmalign_message = subprocess.Popen(f"hmmalign -o {hmmalign_output} {hmm_file} {all_ks_fasta}", shell=True,
                                        stdout=subprocess.PIPE).communicate()
    return hmmalign_output, hmmalign_message


def make_diamond_table(output_dir, processing_dir, threads):
    """
    Make an all vs. all diamond table with the KS sequences
    :param output_dir: String with the path to the output directory
    :param processing_dir: String with the path to the directory where the temporary files are stored
    :param threads: Integer with the number of threads to use
    :return: String with the path to the diamond output file
    """
    all_ks_fasta = Path(output_dir).joinpath('All_KS_sequences.fasta')
    all_ks_diamond = Path(processing_dir).joinpath('All_KS_sequences.dmnd')
    all_ks_dbp = Path(processing_dir).joinpath('All_KS_sequences.dbp')

    sys.stdout.write("\tMaking the diamond database...")
    subprocess.Popen(f"diamond makedb -p {threads} -d {all_ks_diamond} --in {all_ks_fasta}", shell=True,
                     stdout=subprocess.PIPE).communicate()
    sys.stdout.write("done.\n\tAligning the amino acid sequence against the database...")
    subprocess.Popen(f"diamond blastp -d  {all_ks_diamond} -p {threads} -q {all_ks_fasta} -o {all_ks_dbp}", shell=True,
                     stdout=subprocess.PIPE).communicate()
    return all_ks_dbp

def main():
    # Main function of the transPACT rerunner
    # First check if the modules are loaded
    check_loaded_modules()

    # Parse the arguments
    parser = argparse.ArgumentParser(
        description='Run the transPACT pipeline on a folder of antismash annotated Genbanks.')
    parser.add_argument('-input',
                        help='The directories containing the antiSMASH-annotated Genbank files. Multiple directories can be called simultaneously through "-input INPUT1 -input INPUT2".',
                        action='append')
    parser.add_argument('-output', help='The output directory. Default is ./output.', default=Path.cwd() / 'output',
                        type=str)
    parser.add_argument('-threads', help='The number of threads to execute the script on. Default is 1.', default=1,
                        type=int)
    parser.add_argument('-clade_output',
                        help='Filename for the file containing the clade output. Default is clade_output.txt.',
                        default='clade_output.txt', type=str)
    parser.add_argument('-tree_output', help='Filename for the tree Newick output file. Default is tree.nwk',
                        default='tree.nwk', type=str)
    parser.add_argument('-processing_dir', help='Path to save the processing files. Default is ./processing.',
                        default='./processing/', type=str)
    parser.add_argument('-test',
                        help='Run a test case on only 100 files or all files if less than 100 Genbanks are found.',
                        action='store_true')
    parser.add_argument('-build_tree',
                        help='Switch to export the tree. Default is False, due to missing proper installation of ete3 on Morgan.',
                        action='store_true')
    parser.add_argument('-tree_image', help='Filename for tree pdf output file. Default is tree.pdf', type=str,
                        default='tree.pdf')
    parser.add_argument('-remove_processing_files',
                        help='Switch to remove all the fasta and alignment files produced during the processing of the input folder.',
                        action='store_true')
    parser.add_argument('-nrps', help='Switch to perform transPACT on NRPSs instead of trans-AT PKSs.',
                        action='store_true')
    parser.add_argument('--mode', help='The mode to run the script in.', default='client')
    parser.add_argument('--host', default = '127.0.0.1')
    parser.add_argument('--port', default = 52444)
    args = parser.parse_args()

    # Resolve the full locations of the various folders
    output_path = Path(args.output).resolve()
    args.clade_output = output_path.joinpath(args.clade_output).resolve()
    args.tree_output = output_path.joinpath(args.tree_output).resolve()
    args.tree_image = output_path.joinpath(args.tree_image).resolve()
    processing_dir = Path(args.processing_dir).resolve()
    data_dir = Path.cwd().joinpath('transPACT_data')

    # Find the input files 
    files = []
    for path in args.input:
        files.extend([file for file in Path(path).expanduser().glob('./*.gb*')])
    if args.test:
        files = files[0:min(100, len(files))]

    # Extract the KS sequences from the genbanks in the folder
    pool = mp.Pool(args.threads)
    sys.stdout.write('Starting antiSMASH domain extraction...')
    if args.nrps:
        aS_domains = [(files[i], cdss) for i, cdss in
                      enumerate(pool.starmap(find_NRPS_asdomains, [(filepath, 'PKS_KS') for filepath in files]))]
        sys.stdout.write(f"\tFound {len(aS_domains)} clusters.\n")
    else:
        aS_domains = [(files[i], cdss) for i, cdss in
                      enumerate(pool.starmap(find_transAT_asdomains, [(filepath, 'PKS_KS') for filepath in files]))]
        sys.stdout.write(f"\tFound {len(aS_domains)} clusters.\n")
        aS_domains = [(file, cdss) for file, cdss in aS_domains if any(len(pks) > 1 for pks in list(cdss.values()))]
        sys.stdout.write(
            f"\tKept {len(aS_domains)} clusters after filtering for BGC with PKSs containing more than 1 KSs.\n")

    # Filter out the clusters that only have PKSs with 1 KS
    sys.stdout.write('antiSMASH domain extraction done\n')

    # Write the KS sequences to a fasta file
    block_paths, _ = save_fastas(aS_domains, processing_dir, output_path, block_size=10)

    # Run the transPACT annotation on all the KS sequences
    sys.stdout.write('done\n')
    sys.stdout.write('Starting transPACT annotation.\n')
    pool.starmap(substrate_from_faa_submitter, [(file, processing_dir) for file in block_paths])
    sys.stdout.write('done\n')
    pool.close()

    # Obtain the output of substrate_from_faa and add the missing KSs
    sys.stdout.write('Checking for missing KSs from Guppy...')
    parsed_KSs, missing_KSs = parse_and_add_missing_KSs(processing_dir, args.threads)
    sys.stdout.write(f"\n\tAdded {len(missing_KSs)} KSs.\n")

    for file in processing_dir.glob('*|*'):
        print(file)
        os.remove(file)

    write_KS_annotation_matrix_file(parsed_KSs, args.clade_output)

    # Make the hmmalign alignment
    # See supplementary information of BiGSCAPE paper for details

    # Make the Diamond table
    sys.stdout.write("Making the Diamond table...")
    diamond_table = make_diamond_table(output_path, processing_dir, args.threads)
    sys.stdout.write('Done making DIAMOND table.\n')

    # Construct the phylogenetic tree of all GenBank files 
    subprocess.Popen(
        f"python generate_dendrogram_userweights.py -blasttable {diamond_table} -annotation_matrix {args.clade_output} -treefile {args.tree_output}",
        shell=True, stdout=subprocess.PIPE).communicate()

    # For using ete3, cite
    # Jaime Huerta-Cepas, François Serra and Peer Bork. "ETE 3: Reconstruction,analysis and visualization of phylogenomic data."  Mol Biol Evol (2016) doi: 10.1093/molbev/msw046
    if args.build_tree:
        generate_tree(args.input, args.tree_output, args.tree_image, args.clade_output, generate_new_accession_dict=True)

    # Clean processing folder if desired
    if args.remove_processing_files:
        remove_processing_files(processing_dir)
        sys.stdout.write('Deleting processing files.\n')

    sys.stdout.write('===========================================')
    sys.stdout.write('                  FINISHED                 ')
    sys.stdout.write('===========================================')


if __name__ == '__main__':
    main()
