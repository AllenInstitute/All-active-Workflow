from ateamopt.utils import utility

def main():
    feature_set = {
        'features' : [
                      'voltage_base',
                      'steady_state_voltage',
                      'mean_frequency',
                      'time_to_first_spike',
                      'AP_amplitude_from_voltagebase',
                      'ISI_CV',
                      'AP_width',
                      'adaptation_index2',
                      'AHP_depth',
                      'depol_block',
                      'Spikecount'
                      ]
            }
   
    path = 'feature_set_stage2.json'
    utility.save_json(path, feature_set)


if __name__ == '__main__':
    main()
