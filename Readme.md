# Updated ***trans***PACT script
## Description
This repository holds the code and data required to rerun the ***trans***PACT pipeline 
as published by Leopold-Messer ***et al***.
The pipeline is based forked from the original ***trans***PACT pipeline ([Helfrich ***et al***., doi: ```10.1038/s41467-021-21163-x```](https://www.nature.com/articles/s41467-021-21163-x),
[```https://github.com/chevrm/transPACT/```](https://github.com/chevrm/transPACT) and has been updated to run on Python 3 and some minor bugs were fixed.

In addition, the original paper provided two separate scripts: one to annotate the KS phylogeny and one to build the tree. The ```transPACT.py``` script connects these two scripts and provides an integrated pipeline that takes antiSMASH annotated GenBank files as input and produces the phylogenetic tree as output (among others).
## Installation
You can clone with ```git clone ``` or download the code from this webpage. \
The dependencies of the ete3 package, which is required for drawing the phylogenetic tree, are a bit tricky. Therefore, we highly recommend to set up a dedicated ```conda``` environment for this pipeline through the provided ```conda_environment.yml``` file. This can be done by running the following command in the root directory of this repository:
```bash conda env create -f conda_environment.yml```\
After installation has finished, you can activate the environment by running ```conda activate transPACT```. 

## Usage
The pipeline can be run by executing the ```transPACT.py``` script. Run ```python transPACT.py -h``` to see the help message and the list of command line arguments .\


> 