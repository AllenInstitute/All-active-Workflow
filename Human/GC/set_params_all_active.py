#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:54:54 2018

@author: anin
"""

# Setting the parameters for the Granule Cell 

import json

all_params = {   
                'cm' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'bounds':{'soma':[1e-1,5], 'apic':[1e-1,5], 'dend' : [1e-1,5],
                                    'axon':[1e-1,5]},
                        },
                'Ra' : {'section' : ['all'],
                        'bounds' : {'all':[50, 1000]}
                        },
                'g_pas' : {'section' : ['all'],
                        'bounds' : {'all':[1e-7, 1e-2]}
                        },
                'e_pas' : {'section' : ['all'],
                        'bounds' : {'all':[-110, -70]}
                        },
                'gbar_HCN' : {'section' : ['apic', 'dend'],
                        'mechanism': 'HCN',
                        'bounds' : {'apic':[1e-7, 1e-4], 'dend' : [1e-7,1e-4]}
                        },
                'gkbar_Kir21': {'section' : ['axon', 'apic', 'dend'],
                            'mechanism': 'Kir21',
                            'bounds':{'axon':[1e-7,1e-2], 'apic':[1e-7,1e-2], 'dend' : [1e-7,1e-2]}
                              }, 
                'gkbar_Kv11': {'section' : ['axon'],
                            'mechanism': 'Kv11',
                            'bounds':{'axon':[1e-7,1e-2]}
                              },
                'gkbar_Kv14': {'section' : ['axon'],
                            'mechanism': 'Kv14',
                            'bounds':{'axon':[1e-7,1e-2]}
                              },               
                'gkbar_Kv34': {'section' : ['axon'],
                              'mechanism' : 'Kv34',
                              'bounds': {'axon':[1e-7,1e-2]}
                              }, 
                'gkbar_Kv723': {'section' : ['axon'],
                              'mechanism' : 'Kv723',
                              'bounds': {'axon':[1e-7,1e-2]}
                              },        
#                'gbar_Kv2like': {'section' : ['axon'],
#                                'mechanism' : 'Kv2like'
#                                },
                'gbar_Cav12' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'mechanism' : 'Cav12',
                          'bounds':{'soma':[1e-7,1e-2], 'apic':[1e-7,1e-2], 'dend' :[1e-7,1e-2],
                                    'axon':[1e-7,1e-2]}
                              },
                'gbar_Cav13' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'mechanism' : 'Cav13',
                          'bounds':{'soma':[1e-7,1e-2], 'apic':[1e-7,1e-2], 'dend' :[1e-7,1e-2],
                                    'axon':[1e-7,1e-2]}
                              },
                'gkbar_SK2' : {'section' : ['apic', 'axon', 'dend'],
                        'mechanism' : 'SK2',
                        'bounds':{'apic':[1e-7,1e-2], 'dend' :[1e-7,1e-2],'axon':[1e-7,1e-2]},
                                },
                'gkbar_Kv21' : {'section' : ['soma'],
                        'mechanism' : 'Kv21',
                        'bounds':{'soma':[1e-7,1e-2]},
                        },
                'gakbar_BK_gc' : {'section' : ['soma','axon'],
                        'mechanism' : 'Kv21',
                        'bounds':{'soma':[1e-7,1e-2],'axon':[1e-7,1e-2]},
                        },
                'gbar_na8st': {'section' : ['soma', 'axon'],
                       'mechanism': 'na8st',
                       'bounds':{'soma':[1e-4,20], 'axon':[1e-4,20]},
                        },
                'gkbar_Kv42':{'section' : ['apic', 'dend'],
                        'mechanism' : 'Kv42',
                        'bounds':{'apic':[1e-7,1e-2], 'dend':[1e-7,1e-2]},
                        },
                 'gbar_Cav22':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cav22',
                         'bounds':{'soma':[1e-7,1e-2], 'apic':[1e-7,1e-2], 'dend' :[1e-7,1e-2],
                                    'axon':[1e-7,1e-2]}
                         },
                 'gbar_Cav32':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cav32',
                         'bounds':{'soma':[1e-7,1e-2], 'apic':[1e-7,1e-2], 'dend' :[1e-7,1e-2],
                                    'axon':[1e-7,1e-2]}
                         },
                 'tau_Cabuffer':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cabuffer',
                         'bounds':{'soma':[20,1000], 'apic':[20,1000], 'dend' :[20,1000],
                                    'axon':[20,1000]}
                         },
                 'brat_Cabuffer':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cabuffer',
                         'bounds':{'soma':[100,1000], 'apic':[100,1000], 'dend' :[100,1000],
                                    'axon':[100,1000]}
                         },
                }



def main():
    with open('all_param_bounds.json','w') as bound_file:
        json.dump(all_params,bound_file,indent=4)   
        
        
if __name__ == '__main__':
    main()