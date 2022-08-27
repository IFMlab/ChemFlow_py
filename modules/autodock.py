#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import subprocess
import shutil
from chemflow import error
from chemflow import read_ligands
from chemflow import pdbqt_string_to_mol2
import re
import codecs

templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/autodock/'
output = 'output.log'
scoring_functions = ("autodock")
format_rec = 'mol2'
config_file_name = 'config.dpf'
data = {'center': ['X', 'Y', 'Z'], 'grid_points': [60, 60, 60],
        'spacing': 0.375, 'dielectric': -0.1465, 'n_poses': 10, 'scoring': 'autodock',
        'command': 'autodock4 -p config.dpf -l output.log'}


def check(data):
    if not shutil.which('autodock4'):
        error('Autodock4 is not installed or on PATH, command "autodock4" not found')
        exit()
    if not shutil.which("prepare_ligand4.py"):
        error('MglTools is not installed or on PATH, command "prepare_receptor4.py" not found')
        exit()
    for coord in data['center']:
        if type(coord) != float and type(coord) != int:
            error('Center must be a list of float values')
            exit()
    for coord in data['grid_points']:
        if type(coord) != int:
            error('Grid_points must be a list of int values')
            exit()
    float_val = ('spacing', 'dielectric', 'n_poses')
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
    parser.add_argument("--grid_points",
                        metavar='FLOAT',
                        nargs=3,
                        type=float, default=[60, 60, 60],
                        help='List of the XYZ point for the binding grid')
    parser.add_argument("--spacing",
                        metavar='INT',
                        type=int, default=0.375,
                        help='Spacing between grid points')
    parser.add_argument("--dielectric",
                        metavar='FLOAT',
                        type=float, default=-0.1465,
                        help='Dielectric function flag')
    parser.add_argument('-n', "--n_poses",
                        metavar='INT',
                        type=int, default=10,
                        help='Number of docking poses to generate')
    parser.add_argument("-sf", "--scoring",
                        metavar='STRING',
                        default='autodock', choices=scoring_functions,
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


def config_grid(data, protocol):
    global protocol_dir
    protocol_dir = protocol
    os.chdir(protocol)
    check(data)
    config_dir = 'config_files/'
    os.mkdir(config_dir)
    if not os.path.isfile(config_dir+'receptor.pdbqt'):
        subprocess.run(
            f'python2 {shutil.which("prepare_receptor4.py")} -r {"receptor.mol2"} -o {config_dir+"receptor.pdbqt"}',
            shell=True, check=True)  # overwrite by default
    # reading grid file
    with open(templates + 'grid.gpf', 'r') as f:
        raw_grid = f.read()
    # inserting the parameter into the config file
    grid = raw_grid.format(center=data['center'], grid_points=data['grid_points'], spacing=data['spacing'],
                           dielectric=data['dielectric'])
    # saving the grid file
    filename = protocol + config_dir + 'grid.gpf'
    with open(filename, 'w') as file:
        file.write(grid)
    os.chdir(config_dir)
    # run grid parametrizion
    with open('grid.log', 'w') as f:
        subprocess.run('autogrid4 -p grid.gpf', shell=True, check=True,
                   stdout=f, stderr=f)
    os.chdir('..')
    return


def dock(data, protocol):
    config_grid(data, protocol)
    # reading config file
    with open(templates + 'config.dpf', 'r') as f:
        raw_text = f.read()
    # inserting the parameter into the config file
    text = raw_text.format(n_poses=data['n_poses'])
    # saving the config file
    filename = protocol + 'config_files/config.dpf'
    with open(filename, 'w') as file:
        file.write(text)
    return 'config.dpf'


def rescore(data, protocol):
    config_grid(data, protocol)
    # reading config file
    with open(templates + 'rescore.dpf', 'r') as f:
        raw_text = f.read()
    # inserting the parameter into the config file
    text = raw_text.format(n_poses=data['n_poses'])
    # saving the config file
    filename = protocol + 'config_files/config.dpf'
    with open(filename, 'w') as file:
        file.write(text)
    return 'config.dpf'


def results():
    energy = list()
    with open('output.log', 'r') as f:
        text = f.read()
    # if contains epdb it is a rescoring file
    if 'epdb: USER' in text:
        energy_text = text.split('Estimated Free Energy of Binding')[1]
        en_bind = float(energy_text.split()[1])
        # epdb doesn't take into account the unbound energy, so we have to add
        internal = text.split('Final Total Internal Energy')[1]
        en_int = float(internal.split()[1])
        energy.append(en_bind - en_int)
        shutil.copy('ligand.pdbqt', 'output.pdbqt')
        return energy
    # repr(text) is used to avoid special character newline '\n'
    matches = re.findall('DOCKED: MODEL(.*?)DOCKED: ENDMDL', repr(text))
    models = str()
    models_no_escape = str()
    for match in matches:
        if 'Estimated Free Energy of Binding' not in match:
            continue
        energy_line = match.split('\\n')[4]
        energy.append(float(energy_line.split()[8]))
        models = models + 'MODEL' + match.replace('DOCKED: ', '') + 'ENDMDL\n'
        models_no_escape = codecs.decode(models, 'unicode_escape')
    with open('output.pdbqt', 'w') as f:
        f.write(models_no_escape)
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
    conf_dir = '../config_files/'
    for file in os.listdir(conf_dir):
        if file in os.listdir():
            os.remove(file)
        os.symlink(conf_dir+file, file)
    subprocess.run(data['command'], shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

''' 
Other way to use autodock more automatically:
python2 [PYTHON-ENV]/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l ligand.pdbqt -r receptor.pdbqt
autogrid4 -p receptor.gpf
python2 [PYTHON-ENV]/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -l ligand.pdbqt -r receptor.pdbqt
autodock4 -p ligand_receptor.dpf

But this python scripts are not very updated'''