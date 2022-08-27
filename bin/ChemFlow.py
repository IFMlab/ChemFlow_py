#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import os
import sys
import shutil
import importlib
import chemflow_parse
import chemflow as cf

CHEMFLOW_HOME = os.environ.get('CHEMFLOW_PY_HOME')
src_script = CHEMFLOW_HOME + '/modules/'
sys.path.append(src_script)


if __name__ == "__main__":
    var = chemflow_parse.dockflow_parse()
    if var['method'] == 'dock':
        docking = cf.Dock(var['receptor'], var['ligand'], program=var['program'], protocol=var['protocol'],
                          verbose=var['verbose'], header=var['header'], scheduler=var['job_scheduler'],
                          molecule_job=var['molecule_job'])
    else:
        docking = cf.Rescore(var['receptor'], var['ligand'], program=var['program'], protocol=var['protocol'],
                          verbose=var['verbose'], header=var['header'], scheduler=var['job_scheduler'],
                          molecule_job=var['molecule_job'])
    if var['load']:
        previous = cf.load(var['load'], n=0)
        docking.data = previous.data
    for key in docking.data:
        if key in var and var[key]:
            docking.data[key] = var[key]
    if var['config']:
        module = importlib.import_module(var['program'])
        docking.data = module.read_config(var['config'])
    if var['crystal']:
        docking.bounding_shape(var['crystal'], 'both')
    docking.check_input()
    chemflow_parse.resume(docking, var)
    if docking.hpc_scheduler and not docking.molecule_job:
        docking.molecule_job = int(input(f'How many docking for {docking.scheduler} job? '))
    start = input('Do you want to continue? ([y]/n) ')
    if 'n' in start.lower():
        exit()
    if not var['run_only']:
        docking.setup(overwrite=var['overwrite'])
    if var['config']:
        shutil.copyfile(var['config'], docking.protocol + module.config_file_name)
    if var['write_only']:
        exit()
    computed = False
    k_poses = var['k_poses']
    if var['postprocess']:
        computed = True
        k_poses = var['k_poses']
    docking.run(keep_poses=k_poses, computed=computed)
