from ateamopt.utils import utility

feature_set = {
        "features": [
                    "voltage_base", 
                    "steady_state_voltage", 
                    "voltage_deflection_vb_ssse", 
                    "decay_time_constant_after_stim",
                	    "sag_amplitude",
                    "Spikecount"
                    ]
            }


def main():
    path = 'feature_set_stage1.json'
    utility.save_json(path,feature_set)
        
        
if __name__ == '__main__':
    main()