#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:54:54 2018

@author: anin
"""

from ateamopt.utils import utility


passive_params = {   
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
                            }
                    }



def main():
    path = 'param_bounds_stage0.json'
    utility.save_json(path,passive_params)   
        
        
if __name__ == '__main__':
    main()