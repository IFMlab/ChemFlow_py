#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari", "Katia Galentino", "Marco Cecchini"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import numpy as np
from itertools import combinations
import json

n_type = np.float32

##########################
# CONSENSUS COMBINATIONS #
##########################


def combos(protocol_list, consensus_deep, repetition=True):
    combos_dict = {}
    total_len = 0
    # read which scoring function is used in each protocol and store the results
    if not repetition:
        scoring_dict = dict()
        for protocol in protocol_list:
            with open(protocol + '/input_ChemFlow.json', 'r') as f:
                scoring = json.load(f)[-1]['data']['scoring']
                scoring_dict[protocol] = scoring
    for deep in range(2, consensus_deep + 1):
        # don't care about repetition of scoring function
        if repetition:
            combos_dict[deep] = list(combinations(protocol_list, deep))
            total_len += len(combos_dict[deep])
        # check the results and avoid repetitions of scoring functions
        else:
            combos_rep = list(combinations(protocol_list, deep))
            combos_no_rep = list()
            for combo in combos_rep:
                scoring_set = {scoring_dict[protocol] for protocol in combo}
                if len(scoring_set) != len(combo):
                    continue
                combos_no_rep.append(combo)
            combos_dict[deep] = combos_no_rep
            total_len += len(combos_dict[deep])
    return combos_dict, total_len


# all next methods are summarized here:
# https://doi.org/10.1038/s41598-019-41594-3
#####################
# ENRICHMENT FACTOR #
#####################

# enrichment factor
# percentage of hits inside the treshold over the percentage of the molecule inside the trehsold
def enrichment_factor(classification, scores, threshold):
    hits_tot = sum(classification)
    n_tot = len(classification)
    n_mol = int(n_tot//100*threshold)
    # calculate the index from the score, the highest positive score has the first index:
    indexes = scores.argsort()[::-1]
    # order the classification
    clas_ord = classification[indexes]
    n_hits = sum(clas_ord[:n_mol])
    ef = n_hits / n_mol * n_tot / hits_tot
    return ef


#############################
# CONSENSUS BASED ON SCORES #
#############################


# Auto-Scaled Scores
# Normalizes each score between 0 and 1 for each ranking, sum the rankings to get the score
def ass(protocol_list, scores_dict, separator='--'):
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        scores = scores_dict[protocol_list[-1]]
        final_scores = previous_scores + (scores - scores.min()) / (scores.max() - scores.min())
        return final_scores, name
    # normal calculation
    for protocol in protocol_list:
        # read and auto_scale the score
        scores = scores_dict[protocol]
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores),  dtype=n_type)
        final_scores = final_scores + (scores - scores.min()) / (scores.max() - scores.min())
    # average by the number of programs
    return final_scores, name

'''    # check if the previous results is already calculated to speed up the calculation
    previous = '--'.join(protocol_list[:-1])
    if previous in scores_dict:
        prev_scores = scores_dict[previous]
        scores = scores_dict[protocol_list[-1]]
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores))
        final_scores = prev_scores + (scores - scores.min()) / (scores.max() - scores.min())
        return final_scores'''


# Auto-Normalized Score
# Similar to ASS. Divide the scores by the best values for each ranking (best scores=1), sum the rankings
def ans(protocol_list, scores_dict, separator='--'):
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        scores = scores_dict[protocol_list[-1]]
        final_scores = previous_scores + scores / scores.max()
        return final_scores, name
    # normal calculation
    for protocol in protocol_list:
        # read and auto_scale the score
        scores = scores_dict[protocol]
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores),  dtype=n_type)
        final_scores = final_scores + scores / scores.max()
    # average by the number of programs
    return final_scores, name


# Rank by number
# Sum the scores of a ligand for each scoring function and average
def rbn(protocol_list, scores_dict, separator='--'):
    deep = len(protocol_list)
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        scores = scores_dict[protocol_list[-1]]
        final_scores = previous_scores*(deep-1) + scores
        return final_scores / deep, name
    for protocol in protocol_list:
        # read and sum the score
        scores = scores_dict[protocol]
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores),  dtype=n_type)
        final_scores = final_scores + scores
    # average by the number of programs
    return final_scores / deep, name


# Z-score
# The score is scaled using the average (μ) and standard deviation (σ) of the  docking program ranking.
# The final score is the average of the scaled-score among all the scoring functions
def z_score(protocol_list, scores_dict, separator='--'):
    deep = len(protocol_list)
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        scores = scores_dict[protocol_list[-1]]
        final_scores = previous_scores*(deep-1) + (scores-scores.mean()) / scores.std()
        return final_scores / deep, name
    for protocol in protocol_list:
        scores = scores_dict[protocol]
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores),  dtype=n_type)
        final_scores = final_scores + (scores-scores.mean()) / scores.std()
    # average by the number of programs
    return final_scores / deep, name


##############################
# CONSENSUS BASED ON RANKING #
##############################
# These methods require the ranking, the ranking can be calculate with the argosrt, but this way is way more slow
# than using directly the index for each initial ranking. There forse is suggested to add to the score_dict the
# array of the indices with the suffix '_index' as dictionary key.

def check_compute_rank(protocol, scores_dict, suffix='_index'):
    if protocol + suffix in scores_dict:
        rank = scores_dict[protocol + suffix]
    else:
        scores_raw = scores_dict[protocol]
        score_index = scores_raw.argsort()[::-1]
        rank = np.array([int(np.where(score_index == i)[0])+1 for i in range(len(scores_raw))], dtype=n_type)
    return rank


# Exponential Consensus Ranking
# https://doi.org/10.1038/s41598-019-41594-3
# results are combined using an exponential distribution for each individual rank
# assuming all docking programs have same index
def ecr(protocol_list, scores_dict, sigma_percent=5, separator='--', suffix='_index'):
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        rank = check_compute_rank(protocol_list[-1], scores_dict, suffix=suffix)
        sigma = len(rank) / 100 * sigma_percent
        scores = np.exp(-rank / sigma)
        final_scores = previous_scores + scores / sigma
        return final_scores, name
    for protocol in protocol_list:
        rank = check_compute_rank(protocol, scores_dict, suffix=suffix)
        sigma = len(rank) / 100 * sigma_percent
        scores = np.exp(-rank/sigma)
        if type(final_scores) != bool:
            final_scores = np.zeros(len(scores),  dtype=n_type)
        final_scores = final_scores + scores
    return final_scores / sigma, name


# Rank by Rank
# Compute the average ranking position for each target. The score is considered as 1/rank
def rbr(protocol_list, scores_dict, separator='--', suffix='_index'):
    deep = len(protocol_list)
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        rank = check_compute_rank(protocol_list[-1], scores_dict, suffix=suffix)
        # we are considering the negative of the scores, and che negative score has already changes sign
        final_scores = previous_scores*(deep-1) - rank
        return final_scores/deep, name
    for protocol in protocol_list:
        rank = check_compute_rank(protocol, scores_dict, suffix=suffix)
        if type(final_scores) != bool:
            final_scores = np.zeros(len(rank),  dtype=n_type)
        final_scores = final_scores + rank
    # the SCORE is calculated as minor rank. This is done because the best ranked compounds have the lowest index, while
    # they must have highest score when computing AUC. With minus rank, the lowest index is the highest score.
    return - final_scores/deep, name


# Rank by Vote
# Each molecule receives a vote (+1) if it is ranked in the top x% of the ranking
# The final score for each molecule is given by the sum of votes obtained from all the programs.
def rbv(protocol_list, scores_dict, top_ranking=5, separator='--', suffix='_index'):
    final_scores = False
    name = separator.join(protocol_list)
    previous = separator.join(protocol_list[:-1])
    # check if the previous calculation has already been done
    if previous in scores_dict:
        previous_scores = scores_dict[previous]
        rank = check_compute_rank(protocol_list[-1], scores_dict, suffix=suffix)
        top = len(rank)/100*top_ranking
        # vote the results
        votes = np.array([1 if ind < top else 0 for ind in rank], dtype=n_type)
        # we are considering the negative of the scores, and che negative score has already changes sign
        final_scores = previous_scores + votes
        return final_scores, name
    for protocol in protocol_list:
        rank = check_compute_rank(protocol, scores_dict, suffix=suffix)
        top = len(rank)/100*top_ranking
        # vote the results
        votes = np.array([1 if ind < top else 0 for ind in rank],  dtype=n_type)
        if type(final_scores) != bool:
            final_scores = np.zeros(len(votes),  dtype=n_type)
        final_scores = final_scores + votes
    return final_scores, name



