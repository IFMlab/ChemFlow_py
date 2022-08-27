#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import shutil
import subprocess
from chemflow import error
from chemflow import check_file
from sys import exit

templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/'
output = "ranking.csv"
scoring_functions = ("chemplp", "plp", "plp95")
format_rec = 'mol2'
config_file_name = 'config.in'
data = {'center': ['X', 'Y', 'Z'], 'radius': 15, 'speed': 'speed1', 'ants': 20,
        'evap': 2, 'sigma': 1.0, 'cluster_rmsd': 1, 'scoring': 'chemplp',
        'n_poses': 10, 'water': None, 'water_xyzr': ('X', 'Y', 'Z', 'R'),
        'command': 'PLANTS1.2_64bit --mode screen ../config.in'}


def check(data):
    # check command
    if not shutil.which('PLANTS1.2_64bit'):
        error('PLANTS is not installed or on PATH, command "PLANTS1.2_64bit" not found')
        exit()
    # check input
    for coord in data['center']:
        if type(coord) != float and type(coord) != int:
            error('Center must be a list of float values')
            exit()
    float_val = ('radius', 'ants', 'evap', 'sigma', 'cluster_rmsd', 'n_poses')
    for val in float_val:
        if val not in data:
            continue
        if type(data[val]) != float and type(data[val]) != int:
            error(f'Value of {val} is {data[val]}, a number is required')
            exit()
    # check water file
    if 'water' in data and data['water']:
        if not os.path.isfile(data['water']):
            error(f'File with water: "{data["water"]}" not found.')
            exit()
        # check water coordinate
        for coord in data['water_xyzr']:
            if type(coord) != float and type(coord) != int:
                error('water_xyzr must be a list of float values')
                exit()
    if data['scoring'] not in scoring_functions:
        error(f'Scoring Function "{data["scoring"]}" not supported, choose between: {str(scoring_functions)}')
        exit()


def parse(parser):
    parser.add_argument("--radius",
                        metavar='FLOAT',
                        type=float, default=15,
                        help='Radius of the binding sphere')
    parser.add_argument("--speed",
                        metavar='STRING',
                        type=str, choices=('speed1', 'speed2', 'speed4'), default='speed1',
                        help='Speed of the ant during optimization')
    parser.add_argument("--ants",
                        metavar='FLOAT',
                        type=float, default=20,
                        help='Number of ants in the search algorithm')
    parser.add_argument("--evap",
                        metavar='FLOAT',
                        type=float, default=0.15,
                        help='Evaporation rate fo the pheromones during search')
    parser.add_argument("--sigma",
                        metavar='FLOAT',
                        type=float, default=1.0,
                        help='Iteration scaling factor during search')
    parser.add_argument("--cluster_rmsd",
                        metavar='FLOAT',
                        type=float, default=1.0,
                        help='RMSD cutoff for clustering')
    parser.add_argument("--water",
                        metavar='MOL2',
                        type=check_file,
                        help='Mol2 file containing a single water molecule')
    parser.add_argument("--water_xyzr",
                        metavar='FLOAT',
                        nargs=4, type=float,
                        help='The center and the radius of the sphere where the water is allowed to move')
    parser.add_argument('-n', "--n_poses",
                        metavar='INT',
                        type=int, default=10,
                        help='Number of docking poses to generate')
    parser.add_argument("-sf", "--scoring",
                        metavar='STRING',
                        default='chemplp', choices=scoring_functions,
                        help='Scoring function')
    return parser


def read_config(config_file):
    print('[ Reading the config file ]')
    read_data = dict()
    read_data['center'] = list()
    # reading the config file
    with open(config_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        for key in data:
            param = line.split()[1:]  # parameters are after a space
            try:
                if key in line and not line.startswith('#'):
                    if key == 'center':
                        read_data['center'] = [float(x) for x in param]
                    else:
                        # read the values with the same type as the default one
                        read_data[key] = type(data[key])(param[0])
                if 'cluster_structures' in line:  # the key for cluster_structures is n_poses in data
                    read_data['n_poses'] = type(data['n_poses'])(param[0])
            except Exception as e:
                print(e)
                error('Error in reading data, check your config file')
                exit()
                return None
    return read_data


def dock(data, protocol):
    global protocol_dir
    protocol_dir = protocol
    os.chdir(protocol)
    check(data)
    # reading config file
    if not data['water']:
        with open(templates + 'plants/config.in', 'r') as f:
            raw_text = f.read()
        # inserting the parameter into the config file
        text = raw_text.format(scoring=data['scoring'], speed=data['speed'], ants=data['ants'],
                               evap=data['evap'], sigma=data['sigma'], center=data['center'],
                               radius=data['radius'], n_poses=data['n_poses'],
                               cluster_rmsd=data['cluster_rmsd'])
    else:
        with open(templates + 'plants/config_water.in', 'r') as f:
            raw_text = f.read()
        # inserting the parameter into the config file
        text = raw_text.format(scoring=data['scoring'], speed=data['speed'], ants=data['ants'],
                               evap=data['evap'], sigma=data['sigma'], center=data['center'],
                               radius=data['radius'], n_poses=data['n_poses'],
                               cluster_rmsd=data['cluster_rmsd'], water=data['water'],
                               water_xyzr=data['water_xyzr'])
        # saving the water file
        shutil.copyfile(data['water'], protocol + os.sep + data['water'])
    # saving the config file
    filename = protocol + '/config.in'
    with open(filename, 'w') as file:
        file.write(text)

    return 'config.in'


def results():
    energy = list()
    with open('ranking.csv', 'r') as f:
        lines = f.readlines()
    for line in lines[1:]:
        if ',' in line:
            energy.append(float(line.split(',')[1]))
    return energy


def read_computed(protocol_dir, ligand):
    if 'docked_ligands.mol2' in os.listdir():
        os.rename('docked_ligands.mol2', 'computed_ligands.mol2')
    ligs_dir = dict()
    with open(protocol_dir + ligand + os.sep + 'computed_ligands.mol2', 'r') as f:
        computed_lig = f.read()
    computed_lig_correct = computed_lig.replace('_entry_00001', '')
    ligand_list = computed_lig_correct.split('@<TRIPOS>MOLECULE')
    ligand_list.remove('')
    for i in range(len(ligand_list)):
        lig_conf = '@<TRIPOS>MOLECULE' + ligand_list[i]
        molecule_name = ligand.split('_conf_')[0] + f'_conf_{i + 1}'
        ligs_dir[molecule_name] = lig_conf
    return ligs_dir


def rescore(data, protocol):
    global protocol_dir
    protocol_dir = protocol
    os.chdir(protocol)
    check(data)
    # reading config file
    with open(templates + 'plants/rescore.in', 'r') as f:
        raw_text = f.read()
    # inserting the parameter into the config file
    text = raw_text.format(scoring=data['scoring'], center=data['center'],
                           radius=data['radius'])
    # saving the config file
    filename = protocol + '/config.in'
    with open(filename, 'w') as file:
        file.write(text)
    # rewriting the variable command in order to do rescoring
    data['command'] = 'PLANTS1.2_64bit --mode rescore ../config.in'
    return 'config.in'


def compute(data, ligand):
    subprocess.run(data['command'], shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
