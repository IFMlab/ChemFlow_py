#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import shutil
from chemflow import error
from chemflow import check_file
import subprocess
import parmed as pmd
import pytraj as pt

red_start = '\033[0;31m'
color_end = '\033[0;0;m'

templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/gromacs/'
output = 'MM_BSA.dat'
scoring_functions = ("mmgbsa", 'mmpbsa')
charges = ('bcc', 'resp')
format_rec = 'pdb'
data = {'md': False, 'scoring': 'mmgbsa', "receptor_itp": None, 'command': 'gmx', 'ff_dir': 'charmm36-feb2021.ff'}


def check(data):
    if not shutil.which(data['command']):
        error(f'Gromacs is not installed or on PATH, command "{data["command"]}" not found')
        exit()
    if type(data['md']) != bool:
        error('md parameter must be a boolean value')
        exit()
    if data['scoring'] not in scoring_functions:
        error(f'Scoring Function "{data["scoring"]}" not supported, choose between {str(scoring_functions)}')
        exit()
    mm_path = os.environ['CONDA_PREFIX'] + '/bin/'
    if not shutil.which('ante-MMPBSA.py', path=mm_path):
        error(f'ante-MMPBSA.py not found in {mm_path}, check to have installed AmberTools in the conda environment')
        exit()
    if not shutil.which('MMPBSA.py', path=mm_path):
        error(f'MMPBSA.py not found in {mm_path}, check to have installed AmberTools in the conda environment')
        exit()
    data['mm_path'] = mm_path
    return


def parse(parser):
    parser.add_argument("-md",
                        action='store_true',
                        help='Run Molecular Dynamics simulation')
    parser.add_argument("--ff_dir",
                        metavar='DIR',
                        default='charmm36-feb2021.ff',
                        help='Directory inside templates/cgenff/ where are stored the default itp files')
    parser.add_argument("--receptor_itp",
                        metavar='FILE',
                        type=check_file,
                        help='ITP file for the receptor')
    parser.add_argument("-sf", "--scoring",
                        metavar='STRING',
                        default='mmgbsa', choices=scoring_functions,
                        help='Scoring function')
    return parser


def rescore(data, protocol):
    global protocol_dir
    protocol_dir = protocol
    os.chdir(protocol)
    check(data)
    gmx = data['command']
    # copy force field directoryies
    shutil.copytree(templates + 'cgenff/' + data['ff_dir'], protocol + data['ff_dir'])
    shutil.copyfile(templates + 'topol-complex.tmpl', protocol + 'topol-complex.tmpl')
    shutil.copytree(templates + 'MDP', protocol + 'MDP')
    for itp in os.listdir(templates + 'itp_files'):
        shutil.copyfile(templates + f'itp_files/{itp}', protocol + data['ff_dir'] + os.sep + itp)
    if data['receptor_itp']:
        file = data['receptor_itp']
        if not os.path.isfile(file):
            error(f'Itp receptor file "{file}" not found.')
            exit()
        shutil.copyfile(file, protocol + data['ff_dir'] + os.sep + file)
    else:
        # prepare cgenff force field for receptor
        print('Processing the receptor itp file ')
        with open('rec_itp.log', 'w') as f:
            subprocess.run(f'echo "1" | {gmx} pdb2gmx -f receptor.pdb -water none -merge all', shell=True, check=True,
                           stdout=f, stderr=f)
            with open('topol.top', 'r') as f_topol:
                topol = f_topol.read()
            itp = '[ moleculetype ]' + topol.split('[ moleculetype ]')[1]
            itp = itp.split('#endif')[0] + '#endif'
            itp_list = itp.split('\n')
            itp_fin = []
            for line in itp_list:
                if 'Protein' in line:
                    itp_fin.append('receptor             3')
                else:
                    itp_fin.append(line)
            with open(protocol + data['ff_dir'] + '/receptor.itp', 'w') as f_itp:
                f_itp.write('\n'.join(itp_fin))
            if 'posre.itp' in os.listdir(protocol):
                shutil.move(protocol + 'posre.itp', protocol + data['ff_dir'] + '/posre.itp')
            # convert the file gro to pdb and use this file for creating the complex
            subprocess.run(f'echo "0" |  {gmx} trjconv -f conf.gro -s conf.gro -o receptor_gro.pdb', shell=True,
                           check=True,
                           stdout=f, stderr=f)
    # copy mmgbpbsa input file
    mm_template = templates + data['scoring'] + '.template'
    mm_dest = protocol + os.sep + data['scoring'] + '.in'
    shutil.copyfile(mm_template, mm_dest)
    # copy gromacs commands to launch
    with open(templates + 'gmx_minimization.tmpl', 'r') as f:
        gmx_cmds = f.read().strip()
    if data['md']:
        with open(templates + 'gmx_production.tmpl', 'r') as f:
            gmx_cmds = gmx_cmds + '\n' + f.read().strip()
    with open('gmx_commands.txt', 'w') as f:
        f.write(gmx_cmds.format(gmx=data['command']))
    return


def results():
    with open(output, 'r') as f:
        text = f.read()
    text_line = text.split('DELTA TOTAL')[-1]
    energy = [float(text_line.split()[0])]
    return energy


# to improve:extract the last structure of the ligand from the trajectory
def read_computed(protocol_dir, ligand):
    shutil.copyfile('ligand.mol2', 'computed_ligands.mol2')
    ligs_dir = dict()
    with open(protocol_dir + ligand + os.sep + 'computed_ligands.mol2', 'r') as f:
        computed_lig = f.read()
    if '_conf_' in ligand:
        name = ligand
    else:
        name = ligand + '_conf_1'
    keyword = '@<TRIPOS>MOLECULE\n'
    computed_lig = computed_lig.replace(keyword + 'MOL\n', keyword + name + '\n')
    ligs_dir[name] = computed_lig
    return ligs_dir


# prepare for cgenff force field for ligand
def prepare_lig():
    # using python script to  produce MOl itp file (REQUIRE PYTHON 3.7 AND NETWORKX 2.3)
    with open('MOL.str', 'w') as f:
        subprocess.run(f'{templates}cgenff/cgenff/cgenff MOL.mol2', shell=True, check=True,
                       stdout=f, stderr=f)
    with open('MOL_cgenff.log', 'w') as f:
        subprocess.run(f'python3 {templates}cgenff/cgenff_charmm2gmx_py3_nx2.py  MOL MOL.mol2 MOL.str {data["ff_dir"]}',
                       shell=True, check=True, stdout=f, stderr=f)


def traj_conv(trajin, topol, trajout):
    tr = pt.load(trajin, topol)
    tr.strip(':SOL,SOD,CLA')
    tr.autoimage()
    tr.rmsfit()
    tr.save(trajout)
    return


def traj_MM_BSA(data, name):
    mm_path = data['mm_path']
    ante_mmpbsa_py = mm_path + 'ante-MMPBSA.py'
    mmpbsa_py = mm_path + 'MMPBSA.py'
    # convert the file to amber format, if not present
    if f'{name}.nc' not in os.listdir():
        try:
            if f'{name}.prmtop' not in os.listdir():
                # read topology with parmed
                topol_gro = pmd.load_file('topol-complex.top', xyz=f'{name}.gro')
                # save amber topology file
                topol_gro.save(f'{name}.prmtop')
            # convert trajectory
            traj_conv(f'{name}-wrap.xtc', f'{name}.prmtop', f'{name}.nc')
            ## COMMAND LINE OPTIONS
            # subprocess.run(f'parmed -i {templates}parmed_min.txt',shell=True, check=True, stdout=f, stderr=f)
            # subprocess.run(f'$(which cpptraj) -i {templates}traj_min.txt', shell=True, check=True, stdout=f, stderr=f)
        except:
            print(red_start + f'Failed conversion of {name}-wrap.gro' + color_end)
            raise SystemError  # AssertionError
    # MM_BSA
    # prepare the topologies and dry simulation
    with open('MM_BSA.log', 'w') as f:
        # if we dont' remvoe this files ante_mmpbsa doesn't run
        if 'com.top' in os.listdir(): os.remove('com.top')
        if 'rec.top' in os.listdir(): os.remove('rec.top')
        if 'lig.top' in os.listdir(): os.remove('lig.top')
        ante_cmd = f"python3 {ante_mmpbsa_py} -p {name}.prmtop -c com.top -r rec.top -l lig.top -s ':SOL,SOD,CLA' -n ':MOL' --radii=bondi"
        subprocess.run(ante_cmd, shell=True, check=True, stdout=f, stderr=f)
        # write mdin
        mdin_cmd = f'python3 {mmpbsa_py} -O -i ../{data["scoring"]}.in -cp com.top -rp rec.top -lp lig.top -y {name}.nc -o MM_BSA.dat -eo MM_BSA.csv -make-mdins'
        subprocess.run(mdin_cmd, shell=True, check=True, stdout=f, stderr=f)
        # change mdin
        with open('_MMPBSA_gb.mdin', 'r') as f2:
            mdin = f2.read()
        with open('_MMPBSA_gb.mdin', 'w') as f2:
            f2.write(mdin.replace('extdiel=80.0', 'intdiel=4, extdiel=80.0'))
        # complete mmgbsa
        mm_cmd = f'python3 {mmpbsa_py} -O -i ../{data["scoring"]}.in -cp com.top -rp rec.top -lp lig.top -y {name}.nc -o MM_BSA.dat -eo MM_BSA.csv -use-mdins'
        subprocess.run(mm_cmd, shell=True, check=True, stdout=f, stderr=f)
        return


def compute(data, ligand):
    # copy the force field directory in the current directory
    if data['ff_dir'] not in os.listdir():
        shutil.copytree('../' + data['ff_dir'], data['ff_dir'])
    if 'MDP' not in os.listdir():
        os.symlink('../MDP', 'MDP', target_is_directory=True)
        # other option if we want to copy the directory
        # shutil.copytree('../MDP', 'MDP')
    try:
        # convert the molecule name to MOL
        with open('ligand.mol2', 'r') as f:
            lig = f.readlines()
        lig[1] = 'MOL\n'
        with open('MOL.mol2', 'w') as f:
            f.writelines(lig)
        prepare_lig()
    except:
        # convert the molecule name to MOL and chek the mol2, usign obabel to solve possible error in the ligand (for PLANTS)
        subprocess.run(f'$(which obabel) -h -imol2 ligand.mol2 -omol2 -OMOL.mol2 --title MOL', shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        prepare_lig()
    if not os.path.isfile(data['ff_dir'] + '/MOL.itp'):
        shutil.move('mol.itp', data['ff_dir'] + '/MOL.itp')
        shutil.move('mol.prm', data['ff_dir'] + '/MOL.prm')

    # create the complex file, removing all symbols apart from ATOMS
    with open('../receptor_gro.pdb', 'r') as f:
        rec = f.readlines()
        rec_fin = [at for at in rec if 'ATOM' in at]
    with open('mol_ini.pdb', 'r') as f:
        mol = f.readlines()
        mol_fin = [at for at in mol if 'ATOM' in at]
    complex = rec_fin + mol_fin
    with open('complex.pdb', 'w') as f:
        f.writelines(complex)
    if 'topol-complex.top' in os.listdir():
        os.remove('topol-complex.top')
    else:
        shutil.copyfile('../topol-complex.tmpl', 'topol-complex.top')
    with open('../gmx_commands.txt', 'r') as f:
        gmx_cmds = f.readlines()
    if data['md']:
        name = 'prod'
    else:
        name = 'mini'
    # launch all gmx command and collect the output
    with open('gmx.log', 'a') as f:
        for cmd in gmx_cmds:
            subprocess.run(cmd.strip(),
                           shell=True, check=True, stdout=f, stderr=f)
    # if not md convert the results and calculate the enrgy
    traj_MM_BSA(data, name)
    return
