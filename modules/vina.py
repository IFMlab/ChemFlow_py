#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import subprocess
import shutil
import argparse
from chemflow import error
from chemflow import read_ligands
from chemflow import pdbqt_string_to_mol2

templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/vina/'
output = 'output.pdbqt'
scoring_functions = ("vina")
format_rec = 'mol2'
config_file_name = 'config.txt'
data = {'center': ['X', 'Y', 'Z'], 'size': [15, 15, 15], 'scoring': 'vina',
        'exhaustiveness': 8, 'energy_range': 3.0, 'n_poses': 10,
        'command': 'vina --config ../config.txt'}


def check(data):
    if not shutil.which('vina'):
        error('Autodock Vina is not installed or on PATH, command "vina" not found')
        exit()
    if not shutil.which("prepare_ligand4.py"):
        error('MglTools is not installed or on PATH, command "prepare_receptor4.py" not found')
        exit()
    for coord in data['center']:
        if type(coord) != float and type(coord) != int:
            error('Center must be a list of float values')
            exit()
    for coord in data['size']:
        if type(coord) != float and type(coord) != int:
            error('Size must be a list of float values')
            exit()
    float_val = ('exhaustiveness', 'energy_range', 'n_poses')
    for val in float_val:
        if type(data[val]) != float and type(data[val]) != int:
            error(f'Value of {val} is {data[val]}, a number is required')
            exit()
    if data['scoring'] not in scoring_functions:
        error(f'Scoring Function "{data["scoring"]}" not supported, choose between: {str(scoring_functions)}')
        exit()


def parse(parser):
    parser.add_argument("-c", "--center",
                        metavar='FLOAT',
                        nargs=3, type=float,
                        help='List of the XYZ coordinates of the center of binding')
    parser.add_argument("--size",
                        metavar='FLOAT',
                        nargs=3,
                        type=float, default=[15, 15, 15],
                        help='List of the XYZ coordinates of the binding box size')
    parser.add_argument("--exhaustiveness",
                        metavar='INT',
                        type=int, default=8,
                        help='Searching exhaustiveness')
    parser.add_argument("--energy_range",
                        metavar='FLOAT',
                        type=float, default=3.0,
                        help='Searching energy range')
    parser.add_argument('-n', "--n_poses",
                        metavar='INT',
                        type=int, default=10,
                        help='Number of docking poses to generate')
    parser.add_argument("-sf", "--scoring",
                        metavar='STRING',
                        default='vina', choices=scoring_functions,
                        help='Scoring function')
    return parser


def read_config(config_file):
    print('[ Reading the config file ]')
    read_data = dict()
    read_data['center'] = list()
    read_data['size'] = list()
    #reading the config file
    with open(config_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        for key in data:
            try:
                if key in line and not line.startswith('#'):
                    param = line.split('=')[1]  # parameter are after the equal sign
                    if key == 'center' or key == 'size':
                        read_data[key].append(float(param))
                    else:
                        # reading the paramter with the same type as the default value
                        read_data[key] = type(data[key])(param)
                if 'num_modes' in line:  # the key for num_modes is n_poses in data
                    param = line.split('=')[1]
                    read_data['n_poses'] = type(data['n_poses'])(param)
            except Exception as e:
                print(e)
                error('Error in reading data, check your config file')
                exit()
                return None
    return read_data


def config_set_up(data, protocol):
    global protocol_dir
    protocol_dir = protocol
    os.chdir(protocol)
    check(data)
    if not os.path.isfile('receptor.pdbqt'):
        subprocess.run(f'python2 {shutil.which("prepare_receptor4.py")} -r {"receptor.mol2"} -o {"receptor.pdbqt"}',
                       shell=True, check=True)  # overwrite by default


def dock(data, protocol):
    config_set_up(data, protocol)
    # reading config file
    with open(templates + 'config.txt', 'r') as f:
        raw_text = f.read()
    # inserting the parameter into the config file
    text = raw_text.format(center=data['center'], size=data['size'], n_poses=data['n_poses'],
                           energy_range=data['energy_range'], exhaustiveness=data['exhaustiveness'])
    # saving the config file
    filename = protocol + os.sep + 'config.txt'
    with open(filename, 'w') as file:
        file.write(text)
    return 'config.txt'


def rescore(data, protocol):
    dock(data, protocol)
    data['command']= 'vina --score_only --config ../config.txt'
    return 'config.txt'


def results():
    energy = list()
    with open('output.pdbqt', 'r') as f:
        lines = f.readlines()
    key_word = 'RESULT:'
    for line in lines[1:]:
        if key_word in line:
            energy.append(float(line.split(key_word)[1].split()[0]))
    return energy


def read_computed(protocol_dir, ligand):
    lig_path = protocol_dir + ligand + os.sep
    pdbqt_dict = read_ligands(lig_path + 'output.pdbqt')
    mol2_dict = dict()
    for model in pdbqt_dict:
        name = ligand.split('_conf_')[0] + f'_conf_{model}'
        mol2_structure = pdbqt_string_to_mol2(lig_path + 'ligand.mol2', pdbqt_dict[model], ligand_name=name)
        mol2_dict[name] = mol2_structure
    computed_ligs = list(mol2_dict.values())
    with open(lig_path + 'computed_ligands.mol2', 'w') as f:
        f.write('\n'.join(computed_ligs))
    return mol2_dict


def compute(data, ligand):
    if not os.path.isfile('ligand.pdbqt'):
        subprocess.run(
        f'python2 {shutil.which("prepare_ligand4.py")} -l ligand.mol2 -o ligand.pdbqt -U "lps"',
        shell=True, check=True)
    subprocess.run(data['command'], shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
