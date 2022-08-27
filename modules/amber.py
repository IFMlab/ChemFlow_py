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
import tempfile
import subprocess

red_start = '\033[0;31m'
color_end = '\033[0;0;m'

templates = os.environ.get('CHEMFLOW_PY_HOME') + '/templates/config/amber/'
output = 'MM_BSA.dat'
scoring_functions = ("mmgbsa", 'mmpbsa')
charges = ('bcc', 'resp')
format_rec = 'pdb'
data = {'md': False, 'water': 'explicit', 'maxcyc': 1000, 'scoring': 'mmgbsa',
        'charge': 'bcc', 'command': 'pmemd', 'score_protocol': None,
        'init': None, 'trajectory': None}


def check(data):
    if not shutil.which('sander'):
        error('AmberTools 17+  is not installed or on PATH, command "sander" not found')
        exit()
    # pmed command doesn't exist
    if shutil.which('pmemd.cuda'):
        amber_exec = "pmemd.cuda"
        # add an argument for double cuda?
        # if shutil.which('pmemd.cuda_DPFP'):
        #     amber_exec="pmemd.cuda_DPFP"
    elif shutil.which('pmemd.cuda_SPFP'):
        amber_exec = "pmemd.cuda_SPFP"
    elif shutil.which('pmemd.MPI'):
        amber_exec = "pmemd.MPI"
    elif shutil.which('pmemd'):
        amber_exec = "pmemd"
    else:
        error("Amber (pmemd) is not installed or on PATH, changing to SANDER.")
        amber_exec = "sander"
        if shutil.which('sander.MPI'):
            amber_exec = "mpirun sander.MPI"
    print(f'Amber executable: {amber_exec}')
    data['command'] = amber_exec
    if type(data['md']) != bool:
        error('md parameter must be a boolean value')
        exit()
    if type(data['maxcyc']) != int:
        error(f'Value of maxcyc is {data["maxcyc"]}, an integer is required')
        exit()
    if data['water'] not in ('explicit', 'implicit'):
        error(f'Water "{data["water"]}" not supported, choose between "explicit" or "implicit"')
        exit()
    if data['charge'] not in charges:
        error(f'Charge "{data["charge"]}" not supported, choose between {str(charges)}')
        exit()
    if data['charge'] == 'resp' and not shutil.which('g09'):
        error('Gaussian09 not installed on PATH, command "g09" not found')
        exit()
    if data['charge'] == 'bcc' and not shutil.which('antechamber'):
        error('Ambertools not installed on PATH, command "antechamber" not found')
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
    parser.add_argument("--maxcyc",
                        metavar='INT', type=int,
                        default=1000,
                        help='The maximum number of minimization steps')
    parser.add_argument("--water",
                        metavar='STRING',
                        type=str, choices=('implicit', 'explicit'), default='explicit',
                        help='Use explicit or implicit solvent')
    parser.add_argument("--charge",
                        metavar='STRING',
                        choices=('bcc', 'resp'), default='bcc',
                        help='Charge to compute for the ligand')
    parser.add_argument("-sf", "--scoring",
                        metavar='STRING',
                        default='mmgbsa', choices=scoring_functions,
                        help='Scoring function')
    return parser


def rescore(data, protocol):
    check(data)
    os.chdir(protocol)
    global protocol_dir
    protocol_dir = protocol
    tleap_config(data)
    if data['water'] == "explicit":
        init = 'ionized_solvated_SALT'
        if data['md']:
            score_protocol = ['min1', 'min2', 'min3', 'min4', 'heat_nvt', 'heat_npt', 'prod']
            trajectory = "prod.nc"
        else:
            score_protocol = ['min1', 'min2', 'min3', 'min4']
            trajectory = "min4.rst7"
    else:
        init = 'complex'
        if data['md']:
            score_protocol = ['min', 'md']
            trajectory = "md.nc"
        else:
            score_protocol = ["min"]
            trajectory = "min.rst7"
    for file in score_protocol:
        file_path = templates + data['water'] + os.sep + file
        file_dest = protocol_dir + os.sep + file
        with open(file_path + '.template', 'r') as f:
            text = f.read()
            # if data['water']=='implicit' and file=='min':
            # if {MAXCYC} is present in the text, it will put your valur, if it is not
            # present in the file it doesn't do anythin
            text = text.format(MAXCYC=data['maxcyc'])
        with open(file_dest + '.in', 'w') as f:
            f.write(text)
    data['init'] = init
    data['trajectory'] = trajectory
    data['score_protocol'] = score_protocol
    mm_template = templates + data['scoring'] + '.template'
    mm_dest = protocol_dir + os.sep + data['scoring'] + '.in'
    shutil.copyfile(mm_template, mm_dest)
    return


def tleap_explicit():
    # generate ionized_solvated.prmtop
    with open('water.log', 'w') as f:
        subprocess.run('$(which tleap) -f ../tleap_water.in', shell=True, check=True, stdout=f, stderr=f)
    # extract the number of water molecules
    with open('ionized_solvated.prmtop', 'r') as f:
        WAT = f.read().count('WAT')
    # write the number of water molecule in water.dat
    # calculate the total number of ions, Cl- and Na+
    nsalt = round(0.0187 * 0.15 * WAT)
    with open('../tleap_salt.in', 'r') as f:
        salt_in = f.read()
    # insert the number of ions in tleap_salt.in, one for each ions
    salt_in = salt_in.replace('nna', str(nsalt))
    salt_in = salt_in.replace('ncl', str(nsalt))
    with open('../tleap_salt-tot.in', 'w') as f:
        f.write(salt_in)
    # run tleap to have ionized_solvated_salt
    with open('tleap.log', 'w') as f:
        subprocess.run('$(which tleap) -f ../tleap_salt-tot.in', shell=True, check=True, stdout=f, stderr=f)


def tleap_implicit():
    with open('tleap.log', 'w') as f:
        subprocess.run('$(which tleap) -f ../tleap_implicit.in', shell=True, check=True, stdout=f, stderr=f)


###oppure usa la funzione format per string interpolation
def tleap_config(data):
    if data['water'] == 'implicit':
        # read
        with open(templates + 'tleap/implicit.template', 'r') as f:
            implicit = f.read()
        # write
        with open(protocol_dir + '/tleap_implicit.in', 'w') as f:
            f.write(implicit.replace('CHARGE', data['charge']))
    else:
        # read
        with open(templates + 'tleap/water.template', 'r') as f:
            explicit = f.read()
        # write
        with open(protocol_dir + '/tleap_water.in', 'w') as f:
            f.write(explicit.replace('CHARGE', data['charge']))
        # read
        with open(templates + 'tleap/salt.template', 'r') as f:
            salt = f.read()
        # write
        with open(protocol_dir + '/tleap_salt.in', 'w') as f:
            f.write(salt.replace('CHARGE', data['charge']))
    return


def results():
    with open(output, 'r') as f:
        text = f.read()
    text_line = text.split('DELTA TOTAL')[-1]
    energy = [float(text_line.split()[0])]
    return energy


# to improve:extract the last structure of the ligand from the trajectory
def read_computed(protocol_dir, ligand):
    files = os.listdir()
    if 'ligand_bcc.mol2' in files:
        shutil.copyfile('ligand_bcc.mol2', 'computed_ligands.mol2')
    elif 'ligand_resp.mol2' in files:
        shutil.copyfile('ligand_resp.mol2', 'computed_ligands.mol2')
    else:
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


def read_net_charge(file, dat=False):
    with open(file, 'r') as f:
        text = f.read()
    lines = text.split('@<TRIPOS>')[2].split('\n')
    charge_dat = list()
    for line in lines:
        if 'BOND' in line:
            break
        if len(line.split()) >= 8 and line.split()[8]:
            charge_dat.append(float(line.split()[8]))
    if dat:
        return [str(x) + '\n' for x in charge_dat]
    netc = round(sum(charge_dat))
    return netc


def run_charge(path, charge):
    if charge == 'bcc':
        cmd_charge = [
            '$(which antechamber) -fi mol2  -i ligand.mol2 -fo mol2  -o bcc.mol2 -c bcc -eq 2 -s 2 -at gaff2 -rn MOL -dr no -pf y -nc {net_charge}']
    else:
        cmd_charge = [
            '$(which antechamber) -fi mol2 -i ligand.mol2 -fo gcrt -o ligand.gau -ge ligand.gesp -ch ligand -eq 1 -gv 1 -gm %mem=2Gb -gn %nproc=$(nproc --all) -rn MOL -dr no -pf y',
            'g09 <ligand.gau>ligand.gout',
            '$(which antechamber) -fi gout -i ligand.gout -fo mol2 -o resp.mol2 -c resp -s 2 -at gaff2 -rn MOL -pf y -dr y -nc {net_charge}']
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        try:
            shutil.copyfile(path + '/ligand.mol2', tmpdir + '/ligand.mol2')
            for cmd in cmd_charge:
                nc = read_net_charge(tmpdir + '/ligand.mol2')
                subprocess.run(cmd.format(net_charge=nc), shell=True, check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.STDOUT)
            shutil.copyfile(tmpdir + os.sep + charge + '.mol2', path + os.sep + charge + '.mol2')
            os.chdir(path)
        except:
            erdir = path + os.sep + 'error_' + charge
            shutil.rmtree(erdir, ignore_errors=True)
            os.mkdir(erdir)
            for file in os.listdir(tmpdir):
                shutil.copyfile(tmpdir + os.sep + file, erdir + os.sep + file)
            raise ValueError


def compute(data, ligand):
    # run charge
    init = data['init']
    trajectory = data['trajectory']
    score_protocol = data['score_protocol']
    amber_exec = data['command']
    charge = data['charge']
    path = os.getcwd()
    if os.path.isfile(path + os.sep + charge + '.mol2'):
        print('\nCharge already calculated for ' + ligand + '\nDelete the directory to re-compute\n')
    else:
        run_charge(path, charge)
    # check and verify charge
    charge_path = path + os.sep + charge + '.mol2'
    lig_file = path + '/ligand.mol2'
    charge_dat = read_net_charge(charge_path, dat=True)
    charge_dat_path = path + '/charges.dat'
    with open(charge_dat_path, 'w') as f:
        f.writelines(charge_dat)
    # Prepare .mol2 with right charges
    with open('conversion.log', 'w') as f:
        subprocess.run(
            f'$(which antechamber) -i {lig_file} -o ligand_{charge}.mol2 -fi mol2 -fo mol2 -cf charges.dat -c rc -rn MOL -pf yes -dr no -at gaff2',
            shell=True, check=True, stdout=f, stderr=f)
        # Prepare the .frcmod file.
        subprocess.run(f'$(which parmchk2) -i {charge_path} -o {path}/ligand.frcmod -s 2 -f mol2',
                       shell=True, check=True, stdout=f, stderr=f)
    # run tleap
    if data['water'] == 'implicit':
        tleap_implicit()
    else:
        tleap_explicit()

    # run minimization protocol
    prev = init
    for run in score_protocol:
        input_file = f'../{run}'
        with open(f'{run}.log', 'w') as f:
            subprocess.run(
                f'{amber_exec} -O -i {input_file}.in -o {run}.mdout -e {run}.mden -r {run}.rst7 -x {run}.nc -v  {run}.mdvel -inf {run}.mdinfo -c {prev}.rst7 -p {init}.prmtop -ref {prev}.rst7',
                shell=True, check=True, stdout=f, stderr=f)
        prev = run

    # CHEKC MMPBSA AND RUN IN
    mm_path = data['mm_path']
    ante_mmpbsa_py = mm_path + 'ante-MMPBSA.py'
    mmpbsa_py = mm_path + 'MMPBSA.py'
    subprocess.run('rm -rf com.top rec.top ligand.top',
                   shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    with open('MM_BSA.log', 'w') as f:
        subprocess.run(
            f"python3 {ante_mmpbsa_py} -p {init}.prmtop -c com.top -r rec.top -l ligand.top -n :MOL -s ':WAT,Na+,Cl-' --radii=mbondi2",
            shell=True, check=True, stdout=f, stderr=f)
        if data['water'] == 'explicit':
            mm_run = f'python3 {mmpbsa_py} -O -i ../{data["scoring"]}.in -sp {init}.prmtop -cp com.top -rp rec.top -lp ligand.top -o MM_BSA.dat -eo MM_BSA.csv -y {trajectory}'
        else:
            mm_run = f'python3 {mmpbsa_py} -O -i ../{data["scoring"]}.in -cp com.top -rp rec.top -lp ligand.top -o MM_BSA.dat -eo MM_BSA.csv -y {trajectory}'
        subprocess.run(mm_run, shell=True, check=True, stdout=f, stderr=f)
    subprocess.run("rm -rf reference.frc",
                   shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
