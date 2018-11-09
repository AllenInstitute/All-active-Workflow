#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 11:54:39 2018

@author: anin
"""

import os
import json
import pandas as pd


def get_sec(time_str):
    h, m, s = time_str.split(':')
    return float(h) * 3600 + float(m) * 60 + float(s)

def Main():
    
    path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
    with open(path_to_cell_metadata,'r') as metadata:
            cell_metadata = json.load(metadata) 
            
    cell_id = cell_metadata['Cell_id']
    layer = cell_metadata['Layer']
    area = cell_metadata['Area']
    species = cell_metadata['Species']
    cre_line = cell_metadata['Cre_line']
    dendrite_type = cell_metadata['Dendrite_type']
    machine = cell_metadata['Machine']
    
    total = 0
    counter = 0
    time_vec = []
    
    with open('time_info.txt', 'r') as inp, open('total_time.txt', 'w') as outp: 
       for line in inp:
           try:
               num = get_sec(line.strip())
               total += num
               counter += 1
               time_vec.append(num)
           except ValueError:
               print('{} is not a number!'.format(line))
               
       outp.write('#Generations %s: %s seconds'%(counter,total))
    
    print('{} Generations took: {} seconds'.format(counter,total))
    
    time_metrics = pd.DataFrame({'species' : [species for i in range(len(time_vec))],
                        'cell_id' : [cell_id for i in range(len(time_vec))],
                        'layer' : [layer for i in range(len(time_vec))],
                        'area' : [area for i in range(len(time_vec))],
                        'cre_line' : [cre_line for i in range(len(time_vec))],
                        'dendrite_type' : [dendrite_type for i in range(len(time_vec))],
                        'machine' : [machine for i in range(len(time_vec))],
                        'time' : time_vec
                        })
    
    time_metrics_filename = 'time_metrics_'+cell_id+'.csv' 
    time_metrics.to_csv(time_metrics_filename)  
    
if __name__ == '__main__':
    Main()