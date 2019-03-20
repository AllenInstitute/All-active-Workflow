from ateamopt.utils import utility


sup=0.1
inf=10

all_params = {   
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
                'gkbar_Kir21': {'section' : ['axon', 'apic', 'dend'],
                            'mechanism': 'Kir21',
                            'bounds':{'axon':[sup*1.35603352584e-05,inf*0.0303186971763], 'apic':[sup*8.18108406065e-06,inf*7.29395027445e-05], 'dend' : [sup*1.00460816768e-05,inf*0.000486630511843]}
                              }, 
                'gkbar_Kv11': {'section' : ['axon'],
                            'mechanism': 'Kv11',
                            'bounds':{'axon':[sup*0.000123418742201,inf*0.0400476045528]}
                              },
                'gkbar_Kv14': {'section' : ['axon'],
                            'mechanism': 'Kv14',
                            'bounds':{'axon':[sup*0.00128661048642,inf*0.0878442241784]}
                              },               
                'gkbar_Kv34': {'section' : ['axon'],
                              'mechanism' : 'Kv34',
                              'bounds': {'axon':[sup*0.00258730371942,inf*0.178448481113]}
                              }, 
                'gkbar_Kv723': {'section' : ['axon'],
                              'mechanism' : 'Kv723',
                              'bounds': {'axon':[sup*0.000225323920443,inf*0.0679642899908]}
                              },        
                'gbar_Cav12' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'mechanism' : 'Cav12',
                          'bounds':{'soma':[sup*1.7376675793e-05,inf*0.00353917485432], 'apic':[sup*3.04177409039e-06,inf*0.000226752881512], 'dend' :[sup*6.8874637183e-07,inf*0.00263032002298],
                                    'axon':[1.12821431072e-06,inf*0.000842698527358]}
                              },
                'gbar_Cav13' : {'section' : ['soma', 'apic', 'dend', 'axon'],
                          'mechanism' : 'Cav13',
                          'bounds':{'soma':[sup*1.5400937259e-07,inf*0.000635941821159], 'apic':[sup*1.94161629728e-06,inf*0.000172822983153], 'dend' :[sup*3.06423928474e-06,inf*0.000236480590298],
                                    'axon':[sup*1.45387100194e-05,inf*0.000629066689443]}
                              },
                'gkbar_SK2' : {'section' : ['apic', 'axon', 'dend'],
                        'mechanism' : 'SK2',
                        'bounds':{'apic':[sup*1.52260080013e-06,inf*0.000192536416871], 'dend' :[sup*2.85244439822e-07,inf*6.68072472856e-05],'axon':[sup*1.02328012851e-05,inf*0.00168688933165]},
                                },
                'gkbar_Kv21' : {'section' : ['soma'],
                        'mechanism' : 'Kv21',
                        'bounds':{'soma':[sup*9.6025404782e-05,inf*0.0280205630894]},
                        },
                'gakbar_BK_gc' : {'section' : ['soma','axon'],
                        'mechanism' : 'BK_gc',
                        'bounds':{'soma':[sup*0.00105842731675,inf*0.17517529657],'axon':[sup*0.00343693386871,inf*0.396954855292]},
                        },
                'gbar_na8st': {'section' : ['soma', 'axon'],
                       'mechanism': 'na8st',
                       'bounds':{'soma':[sup*0.160913935191,inf*0.832102859566], 'axon':[sup*0.00908143072529,inf*32.2956646203]},
                        },
                'gkbar_Kv42':{'section' : ['apic', 'dend'],
                        'mechanism' : 'Kv42',
                        'bounds':{'apic':[sup*0.00327741227103,inf*0.228060170782], 'dend':[sup*0.0021937383449,inf*0.325187429003]},
                        },
                 'gbar_Cav22':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cav22',
                         'bounds':{'soma':[sup*5.32650855577e-05,inf*0.00236306205943], 'apic':[sup*7.56682447764e-06,inf*0.00395642353009], 'dend' :[sup*5.72817131992e-06,inf*0.00446163756026],
                                    'axon':[sup*1.0345852628e-05,inf*0.00945514284236]}
                         },
                 'gbar_Cav32':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cav32',
                         'bounds':{'soma':[sup*5.85163506309e-06,inf*0.00215020604104], 'apic':[sup*2.02311633564e-06,inf*0.000310549834406], 'dend' :[sup*1.17589732856e-05,inf*0.00133638093067],
                                    'axon':[sup*7.73391424472e-06,inf*0.000654675620729]}
                         },
                 'tau_Cabuffer':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cabuffer',
                         'bounds':{'soma':[10,1000], 'apic':[10,1000], 'dend' :[10,1000],
                                    'axon':[10,1000]}
                         },
                 'brat_Cabuffer':{'section' : ['soma', 'apic', 'dend', 'axon'],
                         'mechanism' : 'Cabuffer',
                         'bounds':{'soma':[10,1000], 'apic':[10,1000], 'dend' :[10,1000],
                                    'axon':[10,1000]}
                         },
                }



def main():
    path = 'param_bounds_stage2.json'
    utility.save_json(path,all_params)   
        
        
if __name__ == '__main__':
    main()