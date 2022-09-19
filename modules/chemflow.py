#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

###########
# MODULES #
###########

import os
import importlib
import json
import shutil
import subprocess
import warnings
from sys import exit
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn import metrics
import matplotlib.pyplot as plt
from alive_progress import alive_bar
# multiprocessing over one machine
import multiprocessing as mp
from functools import partial
# CHEMFLOW modules
import chemflow_parse
import run
import bound_shape
import consensus_functions
import re
import codecs

#####################
# GENERAL VARIABLES #
#####################
# string to color the text output
red = '\033[0;31m'
green = '\033[0;32m'
color_end = '\033[0;0;m'
# working directory
workdir = os.getcwd()
# POSTPROCESS HEADER
header = ['PROGRAM', 'SCORING_FUNCTION', 'PROTOCOL', 'RECEPTOR', 'LIGAND', 'POSE', 'SCORE']


##################
# ERROR HANDLING #
##################

def error(sentence=""):
    """
    It prints the sentence passed to it in red, and then exits the program.

    :param sentence: The sentence to be spoken
    """
    print(red + f'[ ERROR ] {sentence}\n' + color_end)


########################
# GENERAL OS FUNCTIONS #
########################


def check_file(file, file_format=''):
    """
    This function checks if a file exists and if it is in the correct format.
    It returns the file path. If the file path is not specified, the file is inside the workdir.

    :param file: The file to check
    :param file_format: The file format you want to check for
    """
    if not os.path.isfile(file):
        error(f'"{file}", file not found')
        exit()
    if file_format not in file:
        error(f'File "{file}" is not the format required: {format}')
        exit()
    # the user is providing the path to the file
    if workdir in file:
        return file
    # the file is in workdir
    else:
        return workdir + os.sep + file


def check_path_dir(folder):
    """
    This function checks if a directory has the full path and if the last character is a slash.

    :param folder: The folder you want to check
    """
    if workdir not in folder:
        return workdir + os.sep + folder + os.sep
    if os.sep != folder[-1]:
        return folder + os.sep
    return folder


def empty_output(file):
    """
    It returns True if the file is empty, and False otherwise.

    :param file: The file you want to check
    """
    with open(file, 'r') as f:
        n_lines = len(f.readlines())
    if n_lines <= 1:
        return True
    else:
        return False


def check_same_file(file1, file2):
    """
    It checks if two files are the same.

    :param file1: The first file to compare
    :param file2: The file to compare to
    """
    text1 = Path.read_text(Path(file1))
    text2 = Path.read_text(Path(file2))
    if text1 != text2:
        error(f"The file {file1} file doesn't correspond to {file2}. ")
        exit()


def rel_symb_link(source, dest, link_dir):
    """
    It creates a relative symbolic link from the source file to the destination file in the link directory.

    :param source: The file to be linked
    :param dest: The name of the link file
    :param link_dir: The directory where the symbolic link will be created
    """
    rel_path = os.path.relpath(source, link_dir)
    os.symlink(rel_path, dest)


#############################
# USEFUL CHEMFLOW FUNCTIONS #
#############################


def pdbqt_string_to_mol2(original_mol2_file, pdbqt_string, save_as=None, ligand_name='LIGAND', keep_charge=True):
    """
    It takes a pdbqt text and the original mol2 file, and returns a mol2 text with the coordinates of the pdbqt text.

    :param original_mol2_file: The original mol2 file with the atoms' connectivity
    :param pdbqt_string: The string of the pdbqt file
    :param save_as: the name of the file to save the mol2 file as. If None, the file will not be saved
    :param ligand_name: The name of the ligand in the mol2 file, defaults to LIGAND (optional)
    :param keep_charge: If True, the charge of the original mol2 file will be kept. If False, the charge of the pdbqt
                        file will be set, defaults to True (optional)
    """
    check_mol2 = check_file(original_mol2_file, file_format='.mol2')
    with open(check_mol2, 'r') as f:
        mol2_framework = f.readlines()
    # remove empty lines
    if '\n' in mol2_framework:
        mol2_framework = list(filter(lambda a: a != '\n', mol2_framework))
    # rename the mol2
    mol2_framework[1] = ligand_name + '\n'
    # regenerate the mol2 text
    mol2_framework = ''.join(mol2_framework)
    raw_coordinates = pdbqt_string.split('\n')
    atom_list = list()
    for line in raw_coordinates:
        if "ATOM" in line or 'HETATM' in line:
            atom_list.append(line)
    atom_dict = dict()
    for line in atom_list:
        data = line.split()
        atom_index = data[2]
        # creating an index for each atom in which store information
        atom_dict[atom_index] = dict()
        atom_dict[atom_index]['x'] = data[6]
        atom_dict[atom_index]['y'] = data[7]
        atom_dict[atom_index]['z'] = data[8]
        atom_dict[atom_index]['charge'] = data[11]
    # repr(text) is used to avoid special character newline '\n'
    old_atoms_raw = re.search('@<TRIPOS>ATOM\\\\n(.*?)\\\\n@<TRIPOS>BOND', repr(mol2_framework)).group(1)
    old_atoms = codecs.decode(old_atoms_raw, 'unicode_escape').strip('\n')
    old_atoms_list = old_atoms.split('\n')
    if len(old_atoms_list) != len(atom_dict):
        error("[ Pdbqt conversion ] The number of atoms doesn't match")
        return
    new_atoms_list = list()
    for line in old_atoms_list:
        data = line.split()
        atom_index = data[1]
        if atom_index not in atom_dict:
            error(f'[ Pdbqt conversion ] Atom {atom_index} not found in pdbqt')
            return
        data[2] = atom_dict[atom_index]['x']
        data[3] = atom_dict[atom_index]['y']
        data[4] = atom_dict[atom_index]['z']
        if not keep_charge:
            data[8] = atom_dict[atom_index]['charge']
        new_line = '\t'.join(data)
        new_atoms_list.append(new_line)
    new_atoms = '\n'.join(new_atoms_list)
    final_mol2 = mol2_framework.replace(old_atoms, new_atoms)
    if save_as:
        with open(save_as, 'w') as f:
            f.write(final_mol2)
    return final_mol2


def pdbqt_file_to_mol2(original_mol2_file, pdbqt_file, save_as=None, ligand_name='LIGAND', keep_charge=True):
    """
    This function takes a pdbqt file and a mol2 file and returns a mol2 text with the coordinates of the pdbqt file.

    :param original_mol2_file: The original mol2 file with the atoms' connectivity
    :param pdbqt_file: The pdbqt file that you want to convert to mol2
    :param save_as: the name of the file to save the mol2 file as. If None, the file will not be saved
    :param ligand_name: The name of the ligand in the mol2 file, defaults to LIGAND (optional)
    :param keep_charge: If True, the charge of the original mol2 file will be kept. If False, the charge of the pdbqt
                        file will be set, defaults to True (optional)
    """
    check_pdqt = check_file(pdbqt_file, file_format='.pdbqt')
    with open(check_pdqt, 'r') as f:
        raw_coordinates = f.read()
    final_mol2 = pdbqt_string_to_mol2(original_mol2_file, raw_coordinates, save_as, ligand_name=ligand_name,
                                      keep_charge=keep_charge)
    return final_mol2


def ligand_list(file):
    """
    This function takes a file mol2 as an argument and returns a list of ligands.

    :param file: The name of the mol2 file containing the ligands
    """
    try:
        ligands = Path.read_text(Path(file))
    except OSError:
        error(f'Ligand file "{file}" not found')
        exit()
    try:
        splitter = '@<TRIPOS>MOLECULE'
        ligand_list_full = ligands.split(splitter)
        ligands = [lines.split()[0] for lines in ligand_list_full if lines]
        if not ligands:
            error('No ligand found in the ligand file')
    except Exception as e:
        print(e)
        error('Could not read the ligands file, check the format')
    if len(ligands) != len(set(ligands)):
        error('There are duplicates in the ligands file')
    return list(set(ligands))


def read_ligands(ligand_file):
    """
    This function reads in a mol2 or pdbqt ligands file and returns a dictionary, in which the keys are the ligand names
    and the values are the structure for each ligand. The ligands structure maintain the \n characters.

    :param ligand_file: The name of the mol2 or pdbqt file containing the ligands
    """
    ligands_dict = dict()
    ligand_file = check_file(ligand_file)
    ligands = Path.read_text(Path(ligand_file))
    if '.mol2' in ligand_file:
        splitter = '@<TRIPOS>MOLECULE'
    elif '.pdbqt' in ligand_file:
        splitter = 'MODEL'
    else:
        error(f'[ Read ligands ] File format not supported')
        return
    text_ligand_list = ligands.split(splitter)
    for ligand in text_ligand_list:
        # avoid empty elements in the list
        if ligand:
            ligand_name = ligand.split()[0].strip()
            ligands_dict[ligand_name] = splitter + ligand
    return ligands_dict


def load(protocol_path, n=-1, receptor=None, ligands=None):
    """
    This function take a path to a protocol directory (already docked or rescored) and return the Dock or Rescore ChemFlow
    object used to run the calculations of the protocol.

    :param protocol_path: The path to the protocol directory
    :param n: Which ChemFlow object that has been run to load. If n=-1, load the last object; if n=0, it load the first
             object.
    :param receptor: The path to the receptor file
    :param ligands: The path to the ligand file
    """
    protocol = check_path_dir(protocol_path)
    input_path = protocol_path + os.sep + 'input_ChemFlow.json'
    if not os.path.isfile(input_path):
        return error(f'input_ChemFlow.json not found in {protocol}')
    try:
        with open(input_path, 'r') as f:
            previous_data = json.load(f)[n]
            if not ligands:
                ligands = previous_data['ligand']
            if not receptor:
                receptor = previous_data['receptor']
            loaded_step = Step(receptor, ligands,
                               program=previous_data['program'], protocol=previous_data['protocol'])
            loaded_step.__dict__ = previous_data
            loaded_step.protocol = protocol
            loaded_step.name = loaded_step.protocol.split(os.sep)[-2]
    except Exception as e:
        print(e)
        error('Can not read input_ChemFlow.json')
        print(print_dict(previous_data))
        return
    print(green + 'Load completed' + color_end)
    return loaded_step


def check_job(protocol, n_ligands=None, ligands_list=None, ligands_file=None):
    """
    This function checks if a job (docking or rescoring) is finished; returns True if it is finished and False
    otherwise. It compares the number of ligands in the input with the number of ligands on the failed and succeed log
    files.

    :param protocol: the protocol to use
    :param n_ligands: number of ligands to check
    :param ligands_list: a list of ligands to check
    :param ligands_file: a file mol2 containing the ligands to check
    """
    if not n_ligands and not ligands_list and not ligands_file:
        error('To check a job you have to provide the number of ligands or the ligand file or the ligands list')
        return False
    if n_ligands:
        n_ligs = n_ligands
    elif ligands_list:
        n_ligs = len(ligands_list)
    elif ligands_file:
        check_file(ligands_file, '.mol2')
        n_ligs = len(ligand_list(ligands_file))
    protocol = check_path_dir(protocol)
    try:
        n_fail = 0
        n_suc = 0
        fail_file = protocol + 'failed.log'
        suc_file = protocol + 'succeed.log'
        if os.path.exists(fail_file):
            with open(fail_file, 'r') as fail:
                n_fail = len(fail.readlines())
        if os.path.exists(suc_file):
            with open(suc_file, 'r') as suc:
                n_suc = len(suc.readlines())
        n_computed = n_suc + n_fail
        if n_computed == n_ligs:
            return True
        elif n_computed > n_ligs:
            error(f'''The number of ligands in the log files is higher than the provided ligands.
Check:
{suc_file}
{fail_file}''')
            return True
        print(str(n_computed) + ' computed over ' + str(n_ligs))
        return False
    except OSError:
        return False


def consensus(protocol_list, methods=['ass'], folder='consensus'):
    """
    This function takes a list of docking/rescoring protocols and a list of consensus methods and returns the consensus
    rankings

    :param protocol_list: a list of protocols to be used in the consensus
    :param methods: a list of consensus methods to compute the rankings
    :param folder: the folder where the consensus files will be saved, defaults to consensus (optional)
    """
    folder = check_path_dir(folder)
    Path.mkdir(Path(folder), exist_ok=True)
    df0 = read_csv_input(protocol_list[0])[0]
    ligands_ordered = tuple(df0['LIGAND'])
    scores_dict = {}
    # prepare all docked and rescored dataframe
    for protocol in protocol_list:
        df, df_type = read_csv_input(protocol)
        score = np.array(df['SCORE'], dtype=np.float32) * -1
        scores_dict[df_type] = dict(scores_dict.get(df_type, {}), **{protocol: score})
        scores_dict[df_type] = dict(scores_dict.get(df_type, {}),
                                    **{protocol + '_index': np.array(df.index, dtype=np.float32)})
    # run  consensus
    with alive_bar(len(methods)) as bar:
        for cons_method in methods:
            scores_dict[cons_method] = {}
            consensus_df = getattr(consensus_functions, cons_method)
            print('[ Consensus ] Start ' + cons_method)
            temp_dict={}
            for method_dict in scores_dict:
                if method_dict not in cons_method:
                    # add docking and rescoring dictionaries to the temporary dictionary
                    temp_dict.update(scores_dict[method_dict])
            # use a temporary dictionary with the all the score of docking, rescoring and current consensus
            scores, name = consensus_df(protocol_list, temp_dict, separator='--')
            temp_dict.update({name: scores})
            scores_dict[cons_method].update({name: scores})
            bar()
            print(green + '[ Consensus ] computed ' + cons_method + color_end)
    # save csv with consesus results
    for method in methods:
        for name, scores in scores_dict[method].items():
            df = pd.DataFrame(ligands_ordered, columns=['LIGANDS'])
            df['SCORES'] = scores
            df.to_csv(folder + f'{method}_{name}.csv', index=False)
    return scores_dict


####################################
# PROCESSING AND METRICS FUNCTIONS #
####################################


def best(df, element, value, negative=False, ascending=True):
    """
    This function takes a dataframe, the column name with the element to sort,the column name with the values.
    If returns a dataframe with only the elements with the best value.
    If negative=True, the best values is the minimum (used to evaluate the negative scores)
    Example:
    - best(ChemFlow_df, element='LIGAND', value='SCORE', negative=True)
    it returns a dataframe with the best pose for each ligand (the lowest negative score)
    - best(AUC_df, element='TYPE', value='AUC')
    it returns a dataframe with the best docking/rescoring/consensus for each ranking method (the highest AUC value).

    :param df: The dataframe you want to search
    :param element: The column name of the element you want to find the best value for
    :param value: The value you want to compare to
    :param negative: If True, the function will return the lowest value, defaults to False (optional)
    :param ascending: Sort the dataframe in an increasing order (True) or decreasing order (False)
    """
    best_list = []
    # avoid repetition in ligands
    ligs_set = set(df[element])
    for ligand in ligs_set:
        # select all element with the same name of the ligand
        all_elem = df.loc[df[element] == ligand]
        # select the maximum/minimum values
        if negative:
            best_value = min(all_elem[value])
        else:
            best_value = max(all_elem[value])
        # select the poses that have the best value
        best_poses = all_elem.loc[all_elem[value] == best_value]
        # store the best pose, the first line of the best_poses dataframe
        best_list.append(list(best_poses.iloc[0]))
    best_df = pd.DataFrame(best_list, columns=list(df.columns))
    best_df.sort_values(value, ascending=ascending, inplace=True)
    return best_df


def read_csv_input(protocol):
    """
    This function returns the dataframe ChemFlow_best of the protocol and the protocol type
    (the types are 'dock' or 'score').

    :param protocol: The protocol you want to use
    """
    protocol = check_path_dir(protocol)
    df = pd.read_csv(protocol + 'ChemFlow_best.csv')
    df.sort_values(by=['LIGAND'], inplace=True)
    with open(protocol + 'input_ChemFlow.json', 'r') as f:
        df_type = json.load(f)[-1]['type']
    return df, df_type


def classify(ligands, decoys):
    """
    It takes a list of ligands and a set of decoys. If a ligand is a decoy, it is classified as 0, otherwise as 1.
    It returns the classification list.

    :param ligands: a list of ligands
    :param decoys: set of decoys
    """
    # check that the decoys is a set
    if set != type(decoys):
        decoys = set(decoys)
    classification_list = []
    for ligand in ligands:
        if ligand in decoys:
            classification_list.append(0)
        else:
            classification_list.append(1)
    return classification_list


def dict_metrics(metric_dict, name, df_type, scores):
    """
    It takes the metrics dictionary, the method name, the type name and the list of scores. It calculates the AUC score
    and the enrichment factor (EF) for the dataframe and add the name, type and metrics to the metrics dictionary.

    :param metric_dict: a dictionary that will be used to store the AUC values for each fold
    :param name: the name of the method
    :param df_type: the type of the method (docking, rescoring, consensus protocols)
    :param scores: a list of scores
    """
    auc = metrics.roc_auc_score(classification, scores)
    ef_01 = consensus_functions.enrichment_factor(classification, scores, 0.1)
    ef_1 = consensus_functions.enrichment_factor(classification, scores, 1)
    ef_5 = consensus_functions.enrichment_factor(classification, scores, 5)
    ef_10 = consensus_functions.enrichment_factor(classification, scores, 10)
    metric_dict['PROTOCOL'].append(name)
    metric_dict['AUC'].append(auc)
    metric_dict['EF_0.1'].append(ef_01)
    metric_dict['EF_1'].append(ef_1)
    metric_dict['EF_5'].append(ef_5)
    metric_dict['EF_10'].append(ef_10)
    metric_dict['TYPE'].append(df_type)
    return metric_dict


## FIX DOCUMENTATION
def methods_metrics(auc_dict, metrics_dict, folder='analysis/'):
    """
    It takes a dictionary of AUC values and a folder name. Then creates a dataframe of the AUC values for each type of
    ranking (docking, rescoring, consensus, ..). In the end, save the histogram distribution and the csv for each
    ranking type in the specified folder. It automatically overwrites previous existing graphs/csv.

    :param auc_dict: a dictionary of AUC values for each method
    :param folder: the folder where the results are stored
    """
    df = pd.DataFrame.from_dict(auc_dict)
    type_set = {method for method in df['TYPE'] if method not in metrics_dict['METHOD']}
    # avoid annoying warning message from pandas
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method_class in type_set:
            df_type = df[df['TYPE'] == method_class]
            fig, ax = plt.subplots()
            df_type.hist(ax=ax)
            fig.savefig(folder + method_class + '.png')
            plt.close(fig)
            df_metr = df_type.describe()
            for metric in list(df_metr.keys()):
                metrics_dict['METHOD'].append(method_class)
                metrics_dict['metric'].append(metric)
                for index, metric_name in enumerate(list(df_metr.index)):
                    metrics_dict[metric_name].append(df_metr[metric][index])


def analysis(protocol_list, decoy_list, consensus_deep=None, methods=('ass', 'rbr'), verbose=False,
             repetition=False, folder='analysis/'):
    """
    This function takes a list of protocols and a list of decoys and returns a dictionary containing the
    results of the analysis (AUC scores). It uses function 'methods_metrics' to save the metrics and
    results for each scoring method.

    :param protocol_list: a list of the protocols you want to analyze
    :param decoy_list: a list of decoys
    :param consensus_deep: the maximum number of protocol to combine in consensus
    :param methods: a list of consensus method to compute
    :param verbose: if True, prints out the computed consensus method
    :param repetition: if True, the consensus combinations contain repetition in the scoring function,
                        defaults to False (optional)
    :param folder: the folder where the analysis files will be saved, defaults to analysis/ (optional)
    """
    folder = check_path_dir(folder)
    Path.mkdir(Path(folder), exist_ok=True)
    auc_dict = {'PROTOCOL': [], 'AUC': [], 'TYPE': [], 'EF_0.1': [], 'EF_1': [], 'EF_5': [], 'EF_10': []}
    metrics_dict = {'METHOD': [], 'metric': [], 'count': [], 'mean': [], 'std': [],
                    'min': [], '25%': [], '50%': [], '75%': [], 'max': []}
    scores_dict = {}
    df0 = read_csv_input(protocol_list[0])[0]
    ligands_ordered = tuple(df0['LIGAND'])
    global classification
    # using dtype=np.float32 to save space and because scores have not many digits, no need for big storage/preicision
    classification = np.array(classify(ligands_ordered, set(decoy_list)))
    df0['CLASS'] = classification
    df0.loc[:, ['LIGAND', 'CLASS']].to_csv(folder + 'LIGANDS_order.csv', index=False)
    # prepare all docked and rescored dataframe
    for protocol in protocol_list:
        df, df_type = read_csv_input(protocol)
        score = np.array(df['SCORE'], dtype=np.float32) * -1
        scores_dict[df_type] = dict(scores_dict.get(df_type, {}), **{protocol: score})
        scores_dict[df_type] = dict(scores_dict.get(df_type, {}),
                                    **{protocol + '_index': np.array(df.index, dtype=np.float32)})
        dict_metrics(auc_dict, protocol, df_type, score)
    methods_metrics(auc_dict, metrics_dict, folder)
    # calculate the complete number of combinations to run for each method
    print('[ Consensus ] Calculating combinations ')
    combos_dict, total_len = consensus_functions.combos(protocol_list, consensus_deep, repetition=repetition)
    # run  consensus
    for cons_method in methods:
        # import the consensus function
        consensus_df = getattr(consensus_functions, cons_method)
        scores_dict[cons_method] = dict()
        temp_dict = {}
        for method in scores_dict:
            if method not in cons_method:
                # add docking and rescoring dictionaries to the temporary dictionary
                temp_dict.update(scores_dict[method])
        print('[ Consensus ] Start ' + cons_method)
        with alive_bar(total_len) as bar:
            for deep in combos_dict:
                combos = combos_dict[deep]
                df_type = cons_method + '_' + str(deep)
                # tried to do multiprocessing here
                for combo in combos:
                    # use a temporary dictionary with the all the score of docking, rescoring and current consensus
                    scores, name = consensus_df(combo, temp_dict, separator='--')
                    temp_dict.update({name: scores})
                    scores_dict[cons_method].update({name: scores})
                    dict_metrics(auc_dict, name, df_type, scores)
                    bar()
                methods_metrics(auc_dict, metrics_dict, folder)
                if verbose:
                    print(green + 'computed ' + df_type + color_end)
        print(green + '[ Consensus ] computed ' + cons_method + color_end)
    res = pd.DataFrame.from_dict(auc_dict)
    met = pd.DataFrame.from_dict(metrics_dict)
    auc_sorted = res.sort_values(by=['AUC'], ascending=False)
    auc_best = best(auc_sorted, 'TYPE', 'EF_1', negative=False, ascending=False)
    print('[ Analysis ] saving metrics.csv')
    met.to_csv(folder + 'metrics.csv', index=False)
    print('[ Analysis ] saving AUC_sorted.csv and AUC_best.csv')
    auc_sorted.to_csv(folder + 'AUC_sorted.csv', index=False)
    auc_best.to_csv(folder + 'AUC_best.csv', index=False)
    size_fig = 5 + 5 * len(methods) * consensus_deep // 2
    # avoid annoying warning message from pandas
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig, ax = plt.subplots(figsize=(size_fig, size_fig))
        res.hist(by='TYPE', ax=ax)
        fig.savefig(folder + 'histograms.png')
        plt.close(fig)
    return auc_dict


def roc(protocol_list, decoy_list, consensus_tuples=[()], save_as='ROC_curve'):
    """
    It takes a list of protocols, a list of decoys, and a list of consensus tuples, and returns
    a ROC curves graph.

    :param protocol_list: a list of the names of the protocols you want to compare
    :param decoy_list: a list of decoy names
    :param consensus_tuples: a list of tuples of the form [(consensus method, [protocol_list to combine])]
    :param save_as: the name of the file to save the ROC curve as, defaults to 'ROC_curve' (optional)
    """
    auc_dict = {'PROTOCOL': [], 'AUC': [], 'EF_1%': []}
    scores_dict = {}
    df0 = read_csv_input(protocol_list[0])[0]
    ligands_ordered = tuple(df0['LIGAND'])
    classification = np.array(classify(ligands_ordered, set(decoy_list)))
    # prepare all dataframe and insert them in a list 
    for protocol in protocol_list:
        df, df_type = read_csv_input(protocol)
        score = np.array(df['SCORE'], dtype=np.float32) * -1
        scores_dict[protocol] = score
    # run the consensus analysis
    for elem in consensus_tuples:
        cons_method = elem[0]
        combo = elem[1]
        consensus_df = getattr(consensus_functions, cons_method)
        scores, name = consensus_df(combo, scores_dict, separator='--')
        scores_dict[cons_method + '_' + name] = scores
    plt.figure()
    lw = 6
    opacity = 0.7
    for name in scores_dict:
        auc = metrics.roc_auc_score(classification, scores_dict[name])
        ef_1 = consensus_functions.enrichment_factor(classification, scores_dict[name], 1)
        auc_dict['PROTOCOL'].append(name)
        auc_dict['AUC'].append(auc)
        auc_dict['EF_1%'].append(ef_1)
        fpr, tpr, thresholds = metrics.roc_curve(classification, scores_dict[name], pos_label=1)
        plt.plot(fpr, tpr, label=f"{name} (AUC = {auc:.2f}), (EF1% = {ef_1:.2f})", lw=lw, alpha=opacity)
    plt.title("Receiver operating characteristic (ROC)\n", fontweight="bold", fontsize=15)
    plt.plot([0, 1], [0, 1], color="navy", lw=3, linestyle="--")
    plt.tick_params(axis='both', labelsize=15)
    plt.xlabel("\nFalse Positive Rate", fontsize=20)
    plt.ylabel("True Positive Rate\n", fontsize=20)
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if save_as:
        plt.savefig(save_as + ".png", dpi=300, facecolor='white', transparent=False, bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
    plt.show()
    plt.close()
    plt.clf()
    scores_dict["LIGAND"] = ligands_ordered
    scores_dict['CLASS'] = classification
    return scores_dict


def print_dict(param, ind=0, hide=False):
    """
    This function prints a dictionary in a readable format.

    :param param: the dictionary to be printed
    :param ind: indentation level, defaults to 0 (optional)
    :param hide: if True, the function will not print empty values or strings that starts with an underscore
                defaults to False (optional)
    """
    txt = str()
    if type(param) != dict:
        line = str(param)
        if len(line) > 255:
            line = line[:255] + ' ...]'
        return line + '\n'
    if type(param) == dict and ind != 0:
        txt = '\n'
    for key in param:
        if hide and key[0] == '_':
            continue
        if hide and not param[key]:
            continue
        txt = txt + ind * ' ' + key + (13 - len(key)) * ' ' + ': ' + print_dict(param[key], ind + 5)
    return txt


class Step:
    """
    A base class used to run docking or rescoring

    ...

    Attributes
    ----------
  program : str
        the name of the program module to import to run docking/rescoring.
    data : dict
        a dictionary containing the program variable. It is imported from the porgram module.
    protocol : str
        the path of the protocol directory
    receptor : str
        the path to the receptor file
    ligand : str
        the path to the ligand file
    ligands_dir : str
        the path to the directory containing the mol2 ligands extracted from the ligand file. If the ligands_dir is
        equal to the protocol, each ligand file has been stored in its own directory in the protocol folder
    receptor_dir : str
        the path to the directory containing the receptor file.
    ligands : list
        the list of ligands name to compute
    hpc_scheduler : str (slurm or pbs)
        the name of the job scheduler to use for HPC cluster
    header : str
        the path to the file containing the header for HPC job
    molecule_job : int
        number of ligands to compute in a HPC job
    type : str (dock or rescore)
        the type of analysis to run, docking or rescoring
    verbose : bool
        print all the computed ligand on screen


    Methods
    -------
    check_job(self):
            It checks the status of the job and returns the job status (True: finished, False: not finished)

    prepare_directories(self, inplace=True, overwrite=False, ignore=False):
            This function creates the protocol directory, the ligands and receptor directory.

    prepare_receptor(self):
            This function takes a receptor file and prepares it for docking.

    update_ligands(self):
            This function updates the ligands in the system, checking if the ligands files are already computed.
            It will compute only the ligands that are not previously computed (no score detected in the output file)

    prepare_ligands(self, inplace=True, overwrite=False):
            This function takes the ligand file and split all the ligands in single files and directories.

    write(self):
            It writes the file configuration files according to the docking/rescoring method.

    setup(self, inplace=True, overwrite=False, ignore=False):
            This function set up all the files for docking/rescoring: it prepares the directories, receptor file, ligand
            files and configuration files.

    run(self, setup=False, postprocess=True, header=None, scheduler=None, molecule_job=None, computed=False,
                keep_poses=10):
            The function run docking/rescoring both locally or on cluster. It saves the protocol data before running.

    save(self):
            It saves the object parameters inside "input_ChemFlow.json". The file is organized as a list: when the function
            is called, it saves the parameters' dictionary as the last element of the list.

    load(self, n=-1):
            It loads the object parameters stored in input_ChemFlow.json into the current object.

    load_ligand(self, protocol):
            This function loads the computed ligands contained in a protocol.

    bounding_shape(self, crystal, shape='both', padding=15, spacing=0.375):
            This function returns the bounding box of the reference ligand, which is the smallest box that contains the
            crystal structure. It automatically adds all the dimensions to the method' data.

    read_energy_convert_ligand(self, ligand, computed=False, r_name=None):
            This function reads the score of a ligand and eventually convert the docked pose to a mol2 structure.

    postprocess(self, keep_poses=3, computed=True):
            It reads the ligand scores, converts the ligands to mol2 format, and returns csv files with the ligands ranking.

    roc_curve(self, decoy_list, save=True):
            This function takes the list of decoys and calculates the Receiver Operating Characteristic (ROC) curve and the
            Area Under the Curve (AUC)
    """

    def __init__(self, receptor, ligand, protocol=None, program='qvina', header=None, scheduler=None, verbose=False,
                 molecule_job=None):
        global method
        global workdir
        workdir = os.getcwd() + os.sep
        # PROGRAM
        self.program = program
        method = importlib.import_module(program)
        self.data = method.data
        # PROTOCOL
        if not protocol:
            self.protocol = check_path_dir(program)
        else:
            self.protocol = check_path_dir(protocol)
        self.name = self.protocol.split(os.sep)[-2]
        # RECEPTOR
        self.receptor = check_file(receptor)
        # LIGAND
        self.ligand = check_file(ligand, '.mol2')
        # OTHER
        self.ligands_dir = None
        self.receptor_dir = None
        self.ligands = ligand_list(self.ligand)
        self.hpc_scheduler = scheduler
        self.header = None
        if header:
            self.header = check_file(header)
        self.molecule_job = molecule_job
        self.type = 'dock'
        self.verbose = verbose

    def __repr__(self):
        return print_dict(self.__dict__)

    def __str__(self):
        return print_dict(self.__dict__, hide=True)

    def reload_program(self, program):
        """
        It reloads the docking or rescoring program.

        :param program: The program to reload
        """
        global method
        self.program = program
        method = importlib.import_module(program)

    def check_input(self):
        """
        It checks the input data of the docking/rescoring program. It calls the check function of the program.
        """
        method.check(self.data)
        print(green + f'Check completed for {self.program} inputs, no error detected' + color_end)

    def check_job(self):
        """
        It checks the status of the job and returns the job status (True: finished, False: not finished)
        """
        not_comp = 0
        tot = len(self.ligands)
        for ligand in self.ligands:
            lig_path = self.protocol + ligand + os.sep
            files = os.listdir(lig_path)
            # docking/rescoring is not completed in these cases
            if method.output not in files or empty_output(lig_path + method.output):
                not_comp += 1
                if self.verbose:
                    print(ligand)
                continue
            # full check of completion
            energy_list, _ = self.read_energy_convert_ligand(ligand)
            energy_found = len(energy_list)
            # at least one energy found
            if energy_found == 0:
                not_comp += 1
                if self.verbose:
                    print(ligand)
            elif energy_found != int(self.data['n_poses']):
                if self.verbose:
                    print(f'[ Check job ] Found only {str(energy_found)} poses for {ligand}')
        print(f'[ {self.name} ]: computed {tot - not_comp} over {tot}')
        os.chdir(workdir)
        if not_comp == 0:
            return True
        return False

    def prepare_directories(self, inplace=True, overwrite=False, ignore=False):
        """
        This function creates the protocol directory, the ligands and receptor directory.

        :param inplace: If True, the files will be stored in the protocol directory. If False, the files will be
                        copied to the receptor_dir and ligands_dir directories. If these are not defined, the files will
                        be copied to a new directory ('ChemFlow_files'). Then the files are linked with relative symbolic
                        links. Defaults to True (optional)
        :param overwrite: If True, overwrite existing files, defaults to False (optional)
        :param ignore: If True, don't raise an error if the directory already exists, defaults to False (optional)
        """
        # check protocol in case it is changed after initialization of the object
        self.protocol = check_path_dir(self.protocol)
        if os.path.exists(self.protocol) and not ignore:
            if overwrite:
                shutil.rmtree(self.protocol)
            else:
                error(f'Protocol {self.protocol} already exists, use overwrite=True')
                exit()
        print(f'[ Prepare settings ] {self.name}')
        # PREPARE DIRECTORIES
        if inplace:
            self.receptor_dir = self.protocol
            self.ligands_dir = self.protocol
        else:
            # folders are not defined, putting al files in ChemFlow_files
            if not self.receptor_dir:
                self.receptor_dir = os.sep.join(self.protocol.split(os.sep)[:-2]) + 'ChemFlow_files' + os.sep
            else:  # checking the path
                self.receptor_dir = check_path_dir(self.receptor_dir)
            if not self.ligands_dir:
                self.ligands_dir = os.sep.join(self.protocol.split(os.sep)[:-2]) + 'ChemFlow_files' + os.sep
            else:  # checking the path
                self.ligands_dir = check_path_dir(self.ligands_dir)
        dirs = [self.receptor_dir, self.ligands_dir, self.protocol]
        # check if the directories exist and create them
        for fold in dirs:
            Path.mkdir(Path(fold), exist_ok=True)

    def prepare_receptor(self):
        """
        This function takes a receptor file and prepares it for docking. If the receptor file format is different from
        the one required, it tries to convert it with openbabel.
        """
        rec_dir = self.receptor_dir
        rec_format = self.receptor.split('.')[-1]
        # convert receptor in case of different rec_format
        if rec_format != method.format_rec:
            rec_name = self.receptor.split('/')[-1].split('.')[0]
            try:
                print(f'Converting receptor from {rec_format} to {method.format_rec}')
                with open('conversion.log', 'w') as file:
                    subprocess.run(
                        f'$(which obabel) -i{rec_format} {rec_name}.{rec_format} -o{method.format_rec} -O{rec_name}.{method.format_rec} --title receptor',
                        shell=True, check=True, stdout=file, stderr=file)
                self.receptor = workdir + os.sep + f'{rec_name}.{method.format_rec}'
                rec_format = method.format_rec
            except Exception as e:
                print(e)
                error(f'{self.program} requires receptor in {method.format_rec} format')
                exit()
        if rec_format == 'mol2':
            rec_string = 'receptor.mol2'
            rec_mol2 = rec_dir + rec_string
            # copy the receptor.mol2 file if it is not present
            if not os.path.isfile(rec_mol2):
                shutil.copyfile(self.receptor, rec_mol2)  # overwrite by default
            else:  # check that the receptor is the same
                check_same_file(self.receptor, rec_mol2)
            if rec_string not in os.listdir(self.protocol):
                rel_symb_link(rec_mol2, self.protocol + rec_string, self.protocol)

        # prepare the pdb file for mmgbsa
        elif rec_format == 'pdb':
            name_pdb = 'receptor.pdb'
            rec_pdb = rec_dir + os.sep + name_pdb
            if not os.path.isfile(rec_pdb):
                # copy the receptor.pdb file if it is not present
                shutil.copyfile(self.receptor, rec_pdb)  # overwrite by default
            else:  # check that the receptor is the same
                check_same_file(self.receptor, rec_pdb)
            if name_pdb not in os.listdir(self.protocol):
                rel_symb_link(rec_pdb, self.protocol + name_pdb, self.protocol)

    def update_ligands(self):
        """
        This function updates the ligands in the system, checking if the ligands files are already computed.
        It will compute only the ligands that are not previously computed (no score detected in the output file)
        """
        docked_name = {f.name for f in os.scandir(self.protocol) if f.is_dir()}
        common = set(self.ligands).intersection(docked_name)
        # remove common ligands
        ligand_set_new = set(self.ligands) - common
        # add common ligand without output
        for name in common:
            name_path = self.protocol + name
            files = os.listdir(name_path)
            # docking is not completed in these cases
            if method.output not in files or empty_output(name_path + os.sep + method.output):
                ligand_set_new.add(name)
                continue
            # full check of completion
            energy_list, _ = self.read_energy_convert_ligand(name)
            # at least one energy found
            if len(energy_list) == 0:
                ligand_set_new.add(name)
        self.ligands = list(ligand_set_new)

    def prepare_ligands(self, inplace=True, overwrite=False):
        """
        This function takes the ligand file and split all the ligands in single files and directories.

        :param inplace: If True, the ligand files will be written in each ligand directory inside the protocol folder.
                        If False, the ligand files will be written to the ligands_dir folder if specified (written in
                        the "ChemFlow_files" folder otherwise) and linked with relative symbolic links.
                        Defaults to True (optional)
        :param overwrite: If True, will overwrite the existing ligand files, defaults to False (optional)
        """
        self.update_ligands()
        # PREPARE LIGANDS
        ligs_dir = self.ligands_dir
        # reading the ligands and create a dictionary with inchi key as key
        ligands = read_ligands(self.ligand)
        # prepare file that are selected in ligand set
        with alive_bar(len(self.ligands)) as bar:
            for inchi_key in self.ligands:
                inchi_dir = self.protocol + inchi_key
                inchi_mol2 = ligs_dir + inchi_key + '.mol2'
                # check the ligand directory in the protocol
                Path.mkdir(Path(inchi_dir), exist_ok=True)
                if inplace or ligs_dir == self.protocol:  # writing the file in the ligand directory inside protocol
                    inchi_mol2 = inchi_dir + '/ligand.mol2'
                # copy the mol2 file
                if not os.path.isfile(inchi_mol2) or overwrite:
                    with open(inchi_mol2, 'w') as ligand_file:
                        ligand_file.write(ligands[inchi_key])
                # check if the link exist
                if 'ligand.mol2' not in os.listdir(inchi_dir):
                    rel_symb_link(inchi_mol2, inchi_dir + os.sep + 'ligand.mol2', inchi_dir)
                bar()
        os.chdir(workdir)

    def write(self):
        """
        It writes the file configuration files according to the docking/rescoring method.
        """
        if self.type not in dir(method):
            error(f"The method '{self.program}' doesn't contain a function to {self.type}")
            return
        else:
            config = getattr(method, self.type)
            config(self.data, self.protocol)
        os.chdir(workdir)

    def setup(self, inplace=True, overwrite=False, ignore=False):
        """
        This function set up all the files for docking/rescoring: it prepares the directories, receptor file, ligand
        files and configuration files.

        :param inplace: If True, the files will be stored in the protocol directory. If False, the files will be
                        copied to the receptor_dir and ligands_dir directories. If these are not defined, the files will
                        be copied to a new directory ('ChemFlow_files'). Then the files are linked with relative symbolic
                        links. Defaults to True (optional)
        :param overwrite: If True, overwrite existing files, defaults to False (optional)
        :param ignore: If True, ignore the file if it already exists, defaults to False (optional)
        """
        self.prepare_directories(inplace=inplace, overwrite=overwrite, ignore=ignore)
        self.prepare_receptor()
        self.prepare_ligands(inplace=inplace, overwrite=overwrite)
        try:
            self.write()
        except Exception as e:
            print(e)
            print('Error in writing configuration file')
            if not ignore:
                exit()
        print(green + f'[ Prepare settings ] Completed {self.name}' + color_end)
        os.chdir(workdir)

    def run(self, setup=False, postprocess=True, header=None, scheduler=None, molecule_job=None, computed=False,
            keep_poses=10):
        """
        The function run docking/rescoring both locally or on cluster. It saves the protocol data before running.

        :param setup: if True, the setup function will be called with overwrite option, defaults to False (optional)
        :param postprocess: if True, will run the postprocess function, defaults to True (optional)
        :param header: the header file for the cluster
        :param scheduler: The HPC scheduler to use for the cluster. If not specified, it will run locally
        :param molecule_job: The number of molecules to run per job on the cluster
        :param computed: if True, the postprocess function will produce the computed.mol2 file, defaults to False (optional)
        :param keep_poses: number of poses to keep for each ligand for postprocessing, defaults to 10 (optional)
        """
        if setup:
            self.setup(overwrite=True)
        print(f'[ Running ] {self.program} for {self.name}')
        if scheduler:
            self.hpc_scheduler = scheduler
        if self.hpc_scheduler:
            if not molecule_job:
                molecule_job = 1
            self.molecule_job = molecule_job
            if not self.hpc_scheduler or self.hpc_scheduler not in ('pbs', 'slurm'):
                error('HPC scheduler not allowed, chose between "pbs" or "slurm"')
                exit()
            if not self.header:
                self.header = check_file(header)
            self.save()
            run.hpc(self.__dict__, self.protocol, molecule_job)
        else:
            self.save()
            run.run_ligands(self.protocol, 0, len(self.ligands))
        print(green + f'[ Running ] Completed {self.program}' + color_end)
        os.chdir(workdir)
        if not self.hpc_scheduler and postprocess:
            self.postprocess(keep_poses=keep_poses, computed=computed)

    def save(self):
        """
        It saves the object parameters inside "input_ChemFlow.json". The file is organized as a list: when the function
        is called, it saves the parameters' dictionary as the last element of the list.
        """
        input_path = self.protocol + os.sep + 'input_ChemFlow.json'
        if not os.path.isfile(input_path):
            with open(input_path, 'w') as f:
                f.write(json.dumps([self.__dict__]))
        else:
            with open(input_path, 'r+') as f:
                previous_data = json.load(f)
                previous_data.append(self.__dict__)
                # reset the pointer to rewrite the file
                f.seek(0)
                f.write(json.dumps(previous_data))

    def load(self, n=-1):
        """
        It loads the object parameters stored in input_ChemFlow.json into the current object.

        :param n: The index of the input_ChemFlow list to load. 0 means to load the first parameters saved, -1 means
                to load the last parameters saved. Default to -1
        """
        input_path = self.protocol + os.sep + 'input_ChemFlow.json'
        if not os.path.isfile(input_path):
            return error(f'input_ChemFlow.json not found in {self.protocol}')
        try:
            with open(input_path, 'r') as f:
                previus_data = json.load(f)[n]
        except OSError:
            error('Enable to read input_ChemFlow.json')
        self.__dict__ = previus_data
        print(green + 'Load completed' + color_end)

    def load_ligand(self, protocol):
        """
        This function loads the computed ligands contained in a protocol.

        :param protocol: The protocol folder containing the computed ligand file
        """
        file = protocol + os.sep + 'computed_ligands.mol2'
        self.ligand = file
        self.ligands = ligand_list(self.ligand)
        self.molecule_job = len(self.ligands)
        print(green + 'Ligand loading completed' + color_end)

    def bounding_shape(self, crystal, shape='both', padding=15, spacing=0.375):
        """
        This function returns the bounding box of the reference ligand, which is the smallest box that contains the
        crystal structure. It automatically adds all the dimensions to the method' data.

        :param crystal: reference ligand in the binding site
        :param shape: Shape of the grid, 'box', 'sphere' or 'both',  defaults to both (optional)
        :param padding: the number of Angstrom to add to the bounding box, defaults to 15 (optional)
        :param spacing: the spacing between two nodes in the grid (autodock4 parameter), default to 0.375
        """
        if shape not in ('box', 'sphere', 'both'):
            return error(f'Shape {shape} not supported')
        if shape == 'box' or shape == 'both':
            center, size = bound_shape.get_shape(crystal, 'box', padding)
            self.data['size'] = size
            self.data['spacing'] = spacing
            # calculate grid point form the size
            grid = []
            for coord in size:
                # integer division (//) cut the decimals, we want to round the bigger number (+1)
                points = coord // spacing + 1
                # the grid points must be even numbers,
                if points % 2 != 0:
                    points += 1
                grid.append(int(points))
            self.data['grid_points'] = grid
        if shape == 'sphere' or shape == 'both':
            center, radius = bound_shape.get_shape(crystal, 'sphere', padding)
            self.data['radius'] = radius
        self.data['center'] = center

    def read_energy_convert_ligand(self, ligand, computed=False, r_name=None):
        """
        This function reads the score of a ligand and eventually convert the docked pose to a mol2 structure.

        :param ligand: the ligand name
        :param computed: if True, it converts the docked pose to a mol2 structure, defaults to False (optional)
        :param r_name: the name of the receptor
        """
        lig_path = self.protocol + ligand
        os.chdir(lig_path)
        results = list()
        computed_ligands = dict()
        try:
            lig_energy = method.results()
        except OSError:
            error('Error in reading output file for ' + ligand)
            return results, dict()
        if not lig_energy:
            error('Energy not found for ' + ligand)
            return results, dict()
        # creating a list with the csv text-line
        for pose, E in enumerate(lig_energy):
            name = ligand.split('_conf_')[0]
            if self.type == 'rescore' and len(ligand.split('_conf_'))>1:
                pose = int(ligand.split('_conf_')[1]) - 1
            results.append(
                [self.program, self.data['scoring'], self.name, r_name, name, pose + 1, E])
        # read computed ligands
        if computed:
            computed_ligands = method.read_computed(self.protocol, ligand)
        if self.verbose:
            print(green + ligand + color_end + ' ' * 25)
        os.chdir(workdir)
        return results, computed_ligands

    # to implement first poses
    def postprocess(self, keep_poses=3, computed=True):
        """
        It reads the ligand scores, converts the ligands to mol2 format, and returns csv files with the ligands ranking.

        :param keep_poses: The number of poses to keep for each ligand, defaults to 3 (optional)
        :param computed: if True, it will create a file with the mol2 computed structures, defaults to True (optional)
        :return: The final_csv, final_best, final_sorted are being returned.
        """
        rec_name = self.receptor.split(os.sep)[-1]
        print(f'[ PostProcessing ] Postprocessing {self.name} directory')
        results = list()
        computed_ligands = dict()
        pool = mp.Pool()
        worker = partial(self.read_energy_convert_ligand, computed=computed, r_name=rec_name)
        with alive_bar(len(self.ligands)) as bar:
            for mapped_results in pool.imap_unordered(worker, self.ligands):
                results = results + mapped_results[0]
                computed_ligands.update(mapped_results[1])
                bar()
        pool.close()
        # obtain the results dataframe
        final_csv = pd.DataFrame(results, columns=header)
        # keep configurations
        final_csv = final_csv.loc[final_csv['POSE'] <= keep_poses]
        final_csv.sort_values(by=['LIGAND', 'POSE'], inplace=True)
        # sort the other list with the data
        final_sorted = final_csv.sort_values(by=['SCORE'])
        final_best = best(final_sorted, 'LIGAND', 'SCORE', negative=True)
        # save the file
        path_prot = self.protocol + os.sep
        final_csv.to_csv(path_prot + 'ChemFlow.csv', index=False)
        final_sorted.to_csv(path_prot + 'ChemFlow_sorted.csv', index=False)
        final_best.to_csv(path_prot + 'ChemFlow_best.csv', index=False)
        computed_ligands_file = str()
        # create computed ligands file with the poses we select
        if computed:
            for ligand in computed_ligands:
                if int(ligand.split('_conf_')[1]) <= keep_poses:
                    computed_ligands_file = computed_ligands_file + computed_ligands[ligand] + '\n'
            with open(self.protocol + 'computed_ligands.mol2', 'w') as f:
                f.writelines(computed_ligands_file)
        print(green + f'[ PostProcessing ] Completed {self.name}' + color_end)
        os.chdir(workdir)
        return final_csv, final_best, final_sorted

    def roc_curve(self, decoy_list, save=True):
        """
        This function takes the list of decoys and calculates the Receiver Operating Characteristic (ROC) curve and the
        Area Under the Curve (AUC)

        :param decoy_list: the list of the decoys
        :param save: If True, the ROC curve will be saved as a ROC_curve.png file, defaults to True (optional)
        """
        df = pd.read_csv(self.protocol + 'ChemFlow_best.csv')
        df['SCORE'] = df['SCORE'] * -1
        classification = list()
        for ligand in df['LIGAND']:
            if ligand in decoy_list:
                classification.append(0)
            else:
                classification.append(1)
        df['classification'] = classification
        plt.figure(figsize=(8, 6), dpi=100)
        lw = 6
        opacity = 0.7
        auc = metrics.roc_auc_score(df['classification'], df['SCORE'])
        fpr, tpr, thresholds = metrics.roc_curve(df['classification'], df['SCORE'], pos_label=1)
        plt.plot(fpr, tpr, label=f"{self.program} (AUC = %0.2f)" % auc, lw=lw, alpha=opacity)
        plt.title("Receiver operating characteristic (ROC) - Combined\n", fontweight="bold", fontsize=15)
        plt.plot([0, 1], [0, 1], color="navy", lw=3, linestyle="--")
        plt.tick_params(axis='both', labelsize=15)
        plt.xlabel("\nFalse Positive Rate", fontsize=20)
        plt.ylabel("True Positive Rate\n", fontsize=20)
        lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
        if save:
            plt.savefig(self.protocol + "ROC_curve.png", dpi=300, facecolor='white', transparent=False, bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
        return auc


class Dock(Step):
    """
    A class used to run docking.
    Basically is the Step class setting the object.type equal to 'dock'.
    """
    def __init__(self, receptor, ligand, protocol=None, program='qvina', header=None, scheduler=None, verbose=False,
                 molecule_job=None):

        super().__init__(receptor, ligand, protocol=protocol, program=program, header=header, scheduler=scheduler,
                         verbose=verbose, molecule_job=molecule_job)
        self.type = 'dock'


class Rescore(Step):
    """
    A class used to run rescoring.
    Basically is the Step class setting the object.type equal to 'rescore'.
    """
    def __init__(self, receptor, ligand, protocol=None, program='qvina', header=None, scheduler=None, verbose=False,
                 molecule_job=None):

        super().__init__(receptor, ligand, protocol=protocol, program=program, header=header, scheduler=scheduler,
                         verbose=verbose, molecule_job=molecule_job)
        self.type = 'rescore'

        
if __name__ == "__main__":
    var = chemflow_parse.dockflow_parse()
    if var['method'] == 'dock':
        docking = Dock(var['receptor'], var['ligand'], program=var['program'], protocol=var['protocol'],
                          verbose=var['verbose'], header=var['header'], scheduler=var['job_scheduler'],
                          molecule_job=var['molecule_job'])
    elif var['method'] == 'rescore':
        docking = Rescore(var['receptor'], var['ligand'], program=var['program'], protocol=var['protocol'],
                          verbose=var['verbose'], header=var['header'], scheduler=var['job_scheduler'],
                          molecule_job=var['molecule_job'])
    elif var['method'] == 'consensus':
        consensus(var['consensus_list'], var['consensus_method'])
        exit()
    elif var['method'] == 'postprocess':
        for protocol in var['postprocess_list']:
            pp_protocol = load(protocol)
            pp_protocol.postprocess(keep_poses=var['k_poses'])
        exit()
    if var['load']:
        previous = load(var['load'], n=0)
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
    if not var['no_postprocess']:
        docking.run(keep_poses=var['k_poses'], computed=True)

