from ateamopt.utils import utility

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
                                }
                            }

def main():
    path = 'param_bounds_stage1.json'
    utility.save_json(path,passive_and_Ih_params)   
        
        
if __name__ == '__main__':
    main()