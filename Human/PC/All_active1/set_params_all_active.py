#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:54:54 2018

@author: anin
"""

import json

all_params = {   
                'cm' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'bounds':{'soma':[1e-1,10], 'apic':[1e-1,10], 'dend' : [1e-1,10],
                                    'axon':[1e-1,10]},
                        },
                'Ra' : {'section' : ['all'],
                        'bounds' : {'all':[50, 1000]}
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
                'gbar_NaTs2_t': {'section' : ['soma', 'apic', 'dend'],
                            'mechanism': 'NaTs2_t',
                            'bounds':{'soma':[1e-7,1], 'apic':[1e-7, 1e-1], 'dend' : [1e-7,1e-1]}
                              }, 
                'gbar_NaTa_t': {'section' : ['axon'],
                            'mechanism': 'NaTa_t',
                            'bounds':{'axon':[1e-7,4]}
                              },
                'gbar_Nap_Et2': {'section' : ['axon','soma', 'apic', 'dend'],
                            'mechanism': 'Nap_Et2',
                            'bounds':{'axon':[1e-7,4],'soma':[1e-7,1], 'apic':[1e-7, 1], 
                                      'dend' : [1e-7,1]}
                              },               
                'gbar_K_Tst': {'section' : ['axon','soma', 'apic', 'dend'],
                              'mechanism' : 'K_Tst',
                              'bounds': {'axon':[1e-7 ,1e-1],'soma':[1e-7,1e-1],'apic':[1e-7, 1e-1], 
                                      'dend' : [1e-7,1e-1]}
                              }, 
                'gbar_K_Pst': {'section' : ['axon','soma', 'apic', 'dend'],
                              'mechanism' : 'K_Pst',
                              'bounds': {'axon':[1e-7 ,1],'soma':[1e-7,1],'apic':[1e-7, 1], 
                                      'dend' : [1e-7,1]}
                              },        
#                'gbar_Kv2like': {'section' : ['axon'],
#                                'mechanism' : 'Kv2like'
#                                },
                'gbar_Kv3_1' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'mechanism' : 'Kv3_1',
                          'bounds':{'soma':[1e-7,1], 'apic':[1e-7, 1e-1], 'dend' : [1e-7,1e-1],
                                    'axon':[1e-7,2]}
                              },
                'gbar_SK' : {'section' : ['soma', 'axon'],
                        'mechanism' : 'SK',
                        'bounds':{'soma':[1e-7,1e-1], 'axon':[1e-7, 1e-1]},
                            },
                'gbar_Ca_HVA' : {'section' : ['soma', 'axon'],
                        'mechanism' : 'Ca_HVA',
                        'bounds':{'soma':[1e-7,1e-3], 'axon':[1e-7, 1e-3]},
                                },
                'gbar_Ca_LVA' : {'section' : ['soma', 'axon'],
                        'mechanism' : 'Ca_LVA',
                        'bounds':{'soma':[1e-7,1e-2], 'axon':[1e-7, 1e-2]},
                        },
                'gamma_CaDynamics': {'section' : ['soma', 'axon'],
                       'mechanism': 'CaDynamics',
                       'bounds':{'soma':[5e-4,5e-2], 'axon':[5e-4,5e-2]},
                        },
                'decay_CaDynamics':{'section' : ['soma', 'axon'],
                        'mechanism' : 'CaDynamics',
                        'bounds':{'soma':[20,1000], 'axon':[20,1000]},
                        },
                 'gbar_Im':{'section' : ['apic', 'dend'],
                         'mechanism' : 'Im',
                         'bounds':{'apic':[1e-7,1e-3], 'dend':[1e-7,1e-3]},
                         }
                }



def main():
    with open('all_param_bounds.json','w') as bound_file:
        json.dump(all_params,bound_file,indent=4)   
        
        
if __name__ == '__main__':
    main()