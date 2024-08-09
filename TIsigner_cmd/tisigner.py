#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified TIsigner script for integration with Google Colab notebook
"""

import os
import re
import time
import warnings
from multiprocessing import Pool
warnings.filterwarnings("ignore", category=FutureWarning)
import numpy as np
from libs import data, functions

class SubstitutionException(Exception):
    '''Exception when codon substitution range greater then sequence.'''
    pass

def valid_input_seq(seq):
    '''check if given sequence is valid.'''
    seq = re.sub('\s+', '', seq.upper()).rstrip()
    pattern = re.compile('^[ATGCU]*$')
    cod = functions.Optimiser.splitter(seq)
    if list(set(cod[1:-1]) & set(data.STOP_CODONS)):
        raise ValueError('Premature stop codons.')
    elif len(seq) < 75:
        raise ValueError('Sequence too short. Min length = 75 nucleotides.')
    elif cod[0] != 'ATG':
        raise ValueError('No start codon.')
    elif not pattern.match(seq):
        raise ValueError('Unknown nucleotides.')
    elif len(seq)%3 != 0:
        raise ValueError('Sequence is not divisible by 3.')
    return seq

def get_threshold(n):
    '''return threshold opening energy(accessibility) from given score'''
    try:
        if float(n) < 5 or float(n) > 30:
            raise ValueError('Please input opening energy in range 5 to 30 only.')
        threshold = float(n)
        if threshold < 1:
            threshold = 1
        if threshold > 30:
            threshold = 30
    except ValueError:
        raise ValueError('Please input opening energy in range 5 to 30 only.')
    return threshold

def tisigner_optimize(sequence, output='result', codons=9, utr=data.pET21_UTR, 
                      host='ecoli', niter=50, result=20, target_opening_energy=None, 
                      filter_sites=None, termcheck=False, seed=0):
    '''Main function to run TIsigner optimization'''
    
    # Validate input sequence
    s = valid_input_seq(sequence)
    
    # Set up optimization parameters
    try:
        c = int(codons) + 1
        if c*3 >= len(s):
            raise SubstitutionException("Substitution out of range.")
    except (ValueError, SubstitutionException):
        c = int((len(s) - len(s)%3)/3) - 1  # strip stop codon
    
    if host == 'yeast':
        plfold_args = data.RNAPLFOLD_YEAST
    elif host == 'ecoli':
        plfold_args = data.RNAPLFOLD_ECOLI
    else:  # default to E coli
        host = 'ecoli'
        plfold_args = data.RNAPLFOLD_ECOLI
    
    if result > 50:
        result = 50
    
    # Set threshold
    if host == 'ecoli' and utr == data.pET21_UTR:
        threshold = target_opening_energy
    else:
        threshold = None
    
    # Instantiate optimization with given parameters
    seeds = list(range(seed, seed + result))
    rand_states = [np.random.RandomState(i) for i in seeds]
    new_opt = functions.Optimiser(seq=s, host=host, ncodons=c, utr=utr,
                                  niter=niter, threshold=threshold,
                                  rms_sites=filter_sites)
    
    # Run optimizer (multiprocess)
    with Pool(result) as pool:
        results = pool.map(new_opt.simulated_anneal, rand_states)
    
    # Format results
    result_df = functions.sort_results(functions.sa_results_parse(results,
                                       threshold=threshold, termcheck=termcheck),
                                       direction=new_opt.direction, termcheck=termcheck)
    
    # Prepare output
    cols = ['Type', 'Sequence', 'Accessibility', 'pExpressed',
            'Hits', 'E_val', 'Mismatches']
    
    if 'Hits' not in result_df.columns:
        cols.remove('Hits')
        cols.remove('E_val')
    if 'pExpressed' not in result_df.columns:
        cols.remove('pExpressed')
    tmp_df = result_df[cols].copy()
    
    columns_rename = {'pExpressed':'Score',
                      'Accessibility':'Opening Energy',
                      'Sequence':'Sequence', 'Hits':'Term. Hits'}
    tmp_df.rename(columns=columns_rename, inplace=True)
    export_df = tmp_df.reindex(np.roll(tmp_df.index, shift=1)).reset_index(drop=True)
    
    return export_df

# If you need to run TIsigner as a standalone script, you can add this block
if __name__ == '__main__':
    # Add argument parsing and main function call here if needed
    pass
