#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:54:54 2018

@author: anin
"""

import json
import os
import glob

parent_dir = os.path.abspath(os.path.join('.', os.pardir))
path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0]
       
with open(path_to_cell_metadata,'r') as metadata:
    cell_metadata = json.load(metadata)


sup=0.1
inf=10

if cell_metadata['Area'] == 'DG':
    passive_and_Ih_params = { 
                                'cm' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                                        'bounds':{'soma':[0.1,5], 'apic':[0.1,5], 'dend' : [0.1,5],
                                        'axon':[0.1,5]},
                                        },
                                'Ra' : {'section' : ['all'],
                                        'bounds' : {'all':[100, 1000]}
                                        },
                                'g_pas' : {'section' : ['all'],
                                        'bounds' : {'all':[sup*7.77399164824e-08, inf*0.000571761326433]}
                                        },
                                'e_pas' : {'section' : ['all'],
                                        'bounds' : {'all':[-120, -70]}
                                        },
                                'gbar_HCN' : {'section' : ['apic', 'dend'],
                                        'mechanism': 'HCN',
                                        'bounds' : {'apic':[sup*1.41549809896e-08, inf*0.000159936662528],
                                                    'dend' : [sup*4.57696391733e-08,inf*0.000961887464337]}
                                        },
                            }
                                
else:

    passive_and_Ih_params = {   
                                'cm' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                                          'bounds':{'soma':[1e-1,10], 'apic':[1e-1,10], 'dend' : [1e-1,10],
                                                    'axon':[1e-1,10]},
                                        },
                                'Ra' : {'section' : ['all'],
                                        'bounds' : {'all':[50, 200]}
                                        },
                                'g_pas' : {'section' : ['all'],
                                        'bounds' : {'all':[1e-7, 1e-2]}
                                        },
                                'e_pas' : {'section' : ['all'],
                                        'bounds' : {'all':[-110, -60]}
                                        },
                                'gbar_Ih' : {'section' : ['soma', 'apic', 'dend'],
                                        'mechanism': 'Ih',
                                        'bounds' : {'soma':[1e-7,1e-4], 'apic':[1e-7, 1e-4], 'dend' : [1e-7,1e-4]}
                                        },
                                }



def main():
    with open('passive_and_Ih_bounds.json','w') as bound_file:
        json.dump(passive_and_Ih_params,bound_file,indent=4)   
        
        
if __name__ == '__main__':
    main()