#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import argparse
import sys
import math
import socket
import chemflow as cf
import importlib

workdir = os.getcwd()
help_dir = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/help/'


##########
# HEADER #
##########


def chemflow_header():
    with open(help_dir + 'header.txt', 'r') as f:
        print(f.read())


######################
# DOCKFLOW FUNCTIONS #
######################


def dockflow_parse(parser=None):
    chemflow_header()
    # create a parser for command line
    if not parser:
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    help_group = parser.add_argument_group('[ Help ]')
    help_group.add_argument('-h', '--help',
                            action='help',
                            help='show this help message and exit. Add -p PROGRAM to print the help for the program data')
    main = parser.add_argument_group('[ Main ]')
    main.add_argument('method',
                      metavar='dock/rescore/consensus',
                      choices=('dock', 'rescore', 'consensus'), default='dock',
                      help='Compute docking, rescoring or consensus ranking')
    required_mol = True
    if len(sys.argv) > 0 and sys.argv[1] == 'consensus':
        required_mol = False
    main.add_argument("-r", "--receptor",
                      metavar='MOL2/PDB',
                      type=cf.check_file,
                      help='Receptor file in mol2 format',
                      required=required_mol)
    main.add_argument("-l", "--ligand",
                      metavar='MOL2',
                      type=cf.check_file,
                      help='Ligands file in mol2 format',
                      required=required_mol)
    main.add_argument("-p", "--program",
                      metavar='STRING',
                      default='smina',
                      help='Docking/rescoring program')

    optional = parser.add_argument_group('[ Optional ]')
    optional.add_argument("--config",
                          metavar='FILE', type=cf.check_file, choices=('vina', 'smina', 'qvina', 'plants'),
                          help='Configuration file, valid only for programs that require just one config file')
    optional.add_argument('--crystal',
                          metavar='FILE',
                          help='File .mol2 containing a reference ligand structures in the binding site')
    optional.add_argument("--protocol",
                          metavar='STRING',
                          help='Protocol name (default: program name')
    optional.add_argument("--overwrite",
                          action='store_true',
                          help='Overwrite previous protocol')
    optional.add_argument("--write_only",
                          action='store_true',
                          help='Write the directories and configurations without running')
    optional.add_argument("--run_only",
                          action='store_true',
                          help='Compute a protocol without writing configuration')
    optional.add_argument('-v', "--verbose",
                          action='store_true',
                          help='Print on screen the computed ligands')
    optional.add_argument("--load",
                          metavar='DIR',
                          help='Directory with a computed protocol were to load program data')
    # postprocessing
    post = parser.add_argument_group('[ PostProcessing ]')
    post.add_argument("--postprocess",
                      action='store_true',
                      help='Produce a mol2 file with the computed poses')
    post.add_argument('-k', "--k_poses",
                      metavar='INT',
                      type=int, default=3,
                      help='How many poses to keep for postprocessing')
    # Parallel execution
    hpc = parser.add_argument_group('[ High Performance Computing ]')
    hpc.add_argument('-j', "--job_scheduler",
                     metavar='STRING',
                     type=str, choices=('slurm', 'pbs'),
                     help='Job scheduler to run on HPC')
    hpc.add_argument("-mj", '--molecule_job',
                     metavar='INT',
                     type=int,
                     help='Number of molecule to compute each HPC job')
    hpc.add_argument('--header',
                     metavar='File',
                     type=cf.check_file,
                     help='File containing the header text for the HPC submission')
    # consensus options
    consens = parser.add_argument_group('[ Consensus ]')
    consens.add_argument('-pl', '--protocol_list',
                         nargs='+',
                         metavar='DIR',
                         required=not required_mol,
                         help='The list of protocols to combine in consensus')
    consens.add_argument("-cm", "--consensus_method",
                         metavar='STR',
                         nargs='+',
                         default=['z_score', 'rbv'],
                         help='Consensus method to combine the docking rankings')

    progr = parser.add_argument_group('[ Program Data ]')
    if '-p' in sys.argv:
        ind = list(sys.argv).index('-p')
        method = importlib.import_module(sys.argv[ind + 1])
        method.parse(progr)
    elif '--program' in sys.argv:
        ind = list(sys.argv).index('--program')
        method = importlib.import_module(sys.argv[ind + 1])
        method.parse(progr)
    var = vars(parser.parse_args())
    return var


def resume(chem_class, input_var):
    print(f'''
DockFlow summary:
-------------------------------------------------------------------------------
[ General info ]
            HOST: {socket.gethostname()}
            USER: {os.path.expanduser('~').split('/')[-1]}
        PROTOCOL: {chem_class.protocol}
         WORKDIR: {os.getcwd()}
 
[ Docking setup ] 
   RECEPTOR FILE: {chem_class.receptor}
     LIGAND FILE: {chem_class.ligand}
  LIGANDS NUMBER: {len(chem_class.ligands)}
         PROGRAM: {chem_class.program}''')
    for variable, value in chem_class.data.items():
        print(' ' * (16 - len(variable)) + variable.upper() + f': {value}')
        # checking config file
    if input_var['config']:
        print(f'     CONFIG FILE: {input_var["config"]}')

    # checking  job scheduler
    print(f'''
[ Run options ]
   JOB SCHEDULER: {input_var['job_scheduler']}
          HEADER: {input_var['header']}
       OVERWRITE: {input_var['overwrite']}
      WRITE ONLY: {input_var["write_only"]}
        RUN ONLY: {input_var["run_only"]}

[ PostProcess]
     POSTPROCESS: {input_var['postprocess']}
      KEEP POSES: {input_var['k_poses']}''')

    print('-------------------------------------------------------------------------------\n')
