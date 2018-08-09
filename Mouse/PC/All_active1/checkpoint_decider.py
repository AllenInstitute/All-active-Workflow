#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 15:44:22 2018

@author: anin
"""

import pickle
import glob
import numpy as np 
import logging

logger = logging.getLogger(__name__)


def best_seed(cp_dir):
    checkpoint_dir  = cp_dir + '/*.pkl'
    file_list = glob.glob(checkpoint_dir)

    cp_min = float('inf')
    for filename in file_list:
        cp = pickle.load(open(filename,'r'))
        cp_log = cp['logbook']
        cp_log_min =  np.min(np.array(cp_log.select('min')))
        if cp_log_min < cp_min:
            cp_min = cp_log_min
            cp_best = filename
             
    logger.debug('Best checkpoint file is %s with min. objective %s'\
                     %(cp_best.split('/')[-1],cp_min))
    return cp_best