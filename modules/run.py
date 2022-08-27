#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import subprocess
import os
import sys
import argparse
from alive_progress import alive_bar
# multiprocessing over one machine
from multiprocessing import Pool as mp_Pool
# multiprocssing over nodes
from ray.util.multiprocessing import Pool as ray_Pool
import json
import importlib


# variable to check whohc multiprocessing tools to use
multinodes = False

CHEMFLOW_HOME = os.environ.get('CHEMFLOW_PY_HOME')
# this line are a precaution in case run is running on HPC
src_script = CHEMFLOW_HOME + '/modules/'
sys.path.append(src_script)

red_start = '\033[0;31m'
green_start = '\033[0;32m'
yellow_start = '\033[0;33m'
color_end = '\033[0;0;m'
purple_start = '\033[0;35m'


def compute(ligand):
    """
    It takes a ligand as input, compute docking/rescoring (using the function in the docking/rescoring module) and
    return a string with the ligand name colored in green (succeed) or red (failed)

    :param ligand: The ligand to be docked
    """
    path = protocol_dir + os.sep + ligand
    os.chdir(path)
    try:
        method.compute(data, ligand)
    except:
        return 'failed ' + ligand
    return 'computed ' + ligand


def run_ligands(protocoldir, start, end):
    """
    This function runs the ligands in the range from start to end. It creates a multiprocessing pool and run the
    'compute' function for each ligand. If there are errors during the run, the ligand name is saved in 'faliled.log',
    otherwise is saved in 'succeed.log'

    :param protocoldir: the directory where the protocol files are located
    :param start: the starting ligand index
    :param end: the last ligand index to be processed
    """
    with open(protocoldir + os.sep + 'input_ChemFlow.json', 'r') as f:
        step = json.load(f)[-1]
    global protocol_dir
    global method
    global data
    if 'data' in step:
        data = step['data']
    else:
        data = step
    protocol_dir = protocoldir
    method = importlib.import_module(step['program'])
    ligands = step['ligands']
    if not multinodes:
        pool = mp_Pool()
    else:
        pool = ray_Pool()
    os.chdir(protocol_dir)
    try:
        with alive_bar(len(ligands[start:end])) as bar:
            for mapped_results in pool.imap_unordered(compute, ligands[start:end]):
                ligand = mapped_results.split()[1]
                if 'failed' in mapped_results:
                    print(red_start + '[ Error ] ' + ligand + color_end)
                    with open(protocol_dir + '/failed.log', 'a') as f:
                        f.write(ligand + '\n')
                else:
                    if step['verbose']:
                        print(green_start + 'Computed  ' + ligand + color_end)
                    with open(protocol_dir + '/succeed.log', 'a') as f:
                        f.write(ligand + '\n')
                bar()
        pool.close()
    except KeyboardInterrupt:
        print(red_start + "Caught KeyboardInterrupt, terminating workers" + color_end)
        pool.close()
        pool.terminate()


def hpc(var, protocol_dir, molecule_job):
    """
    This function takes a dictionary with the data, the run directory and the number of moleculer to run for each job.
    It submits the jobs on the HPC. The submitted files contains the header and a line to run 'run.py' with
    parameter to be parsed

    :param var: a dictionary of variables that are used in the protocol
    :param protocol_dir: the directory where the protocol is located
    :param molecule_job: the number of molecule to run each job
    """
    runpy = os.path.abspath(__file__)
    with open(var['header'], 'r') as f:
        head = f.read()
    start = 0
    while start < len(var['ligands']):
        end = start + molecule_job
        file = protocol_dir + os.sep + f'ChemFlow_{start}_{end}.' + var['hpc_scheduler']
        body = f'python3 {runpy} -p {protocol_dir} -s {start} -e {end} --ligands {" ".join(var["ligands"])} --multinodes'
        with open(file, 'w') as f:
            f.write(head + '\n' + body)
        try:
            if var['hpc_scheduler'] == 'slurm':
                subprocess.run(f'sbatch {file}', shell=True, check=True)
            else:
                subprocess.run(f'qsub {file}', shell=True, check=True)
            print(green_start + f'job {file} submitted' + color_end)
        except:
            print(red_start + f'[ Error ] job {file} not submitted' + color_end)
        start = end


def run_parse():
    # create a parser for command line
    parser = argparse.ArgumentParser()

    # required arguments
    required = parser.add_argument_group('[ Required ]')

    required.add_argument('-p', "--protocol_dir",
                          metavar='DIR',
                          required=True,
                          help='Protocol directory')

    required.add_argument('-s', "--start",
                          metavar='INT',
                          type=int,
                          help='Start with molecule n')

    required.add_argument('-e', "--end",
                          metavar='INT',
                          type=int,
                          help='End with molecule n')

    optional = parser.add_argument_group('[ Optional ]')

    optional.add_argument("-l", "--ligands",
                          metavar='MOLS',
                          nargs='+',
                          help='molecules list to analyze')

    optional.add_argument("--multinodes",
                        action='store_true',
                        help='Run the ChemFlow over many nodes')

    args = parser.parse_args()
    var = vars(args)
    return var


if __name__ == "__main__":
    var = run_parse()
    multinodes = var['multinodes']
    run_ligands(var['protocol_dir'], var['start'], var['end'])
