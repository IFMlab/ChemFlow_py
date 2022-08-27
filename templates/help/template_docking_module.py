import os
import subprocess
import shutil
import argparse
from chemflow import error
from chemflow import read_ligands
from chemflow import pdbqt_string_to_mol2

# the directory containing the template for this module
templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/module_name/'
# the name of the output file, to check the docking/rescoring completion
output = 'output.pdbqt'
# a list with the supported scoring function
scoring_functions = ("scoring_function1", "scoring_function2")
# format rquired for the receptor
format_rec = 'pdb'
# dictionary containing the data for docking, with also the default values
data = {'center': ['X', 'Y', 'Z'], 'PARAMETERS': '...',
        'command': ['PROGRAM --CONFIG_FILE config.txt']}

# FUNCTION TO CHECK THE INPUT PARAMETERS
# REQUIRED
def check(data):
    return

# PARSER FOR THE COMMAND LINE
# REQUIRED FOR COMMAND LINE
def parse(parser):
    return

# PARSE CONFIGURATION FILE PROVIDED BY THE USER
# OPTIONAL
def read_config(config_file):
    return

# FUNCTION TO CREATE A CONFIGURATION FILE FOR DOCKING IN THE PROTOCOL DIRECTORY
# REQUIRED
def dock(data, protocol):
    return

# FUNCTION TO CREATE A CONFIGURATION FILE FOR RESCORING IN THE PROTOCOL DIRECTORY
# OPTIONAL, if the program admits rescoring
def rescore(data, protocol):
    # CHANGE THE COMMAND TO RUN RESCORING
    data['command']= ['PROGRAM --CONFIG_FILE config.txt --RESCORE']
    return

# FUNCTION TO READ THE SCORE FOR EACH LIGAND POSE FROM THE OUTPUT FILE
# REQUIRED
def results():
    return 'LIST OF SCORES'

# FUNCTION TO CONVERT THE DOCKED POSE TO MOL2 STRUCTURES
# REQUIRED
def read_computed(protocol_dir, ligand):
    return 'DICTIONARY WITH LIGAND_NAME AS KEYS AND MOL2 STRUCTURES TEXTS AS VALUES'

# FUNCTION TO PREPARE THE LIGAND AND RUN DOCKING/RESCORING
# REQUIRED
def compute(data, ligand):
    # PREPARE THE LIGAND...
    # RUN DOCKING/RESCORING, EXAMPLE:
    subprocess.run(data['command'], shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

