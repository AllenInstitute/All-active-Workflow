import os
from allensdk.core.nwb_data_set import NwbDataSet
import numpy as np
import json
from collections import defaultdict
import efel
import math
from ateamopt.utils import utility
import logging
logger = logging.getLogger(__name__)

class NWB_Extractor(object):
    
    def __init__(self, cell_id, junc_potential=-14,temp=34,**kwargs):
        
        self.cell_id = cell_id
        self.junction_potential = junc_potential
        self.temperature = temp
        nwb_dir =  utility.get_filepath_for_exten(exten = '.nwb')
        self.nwb_path =[str_path for str_path in nwb_dir if 'cell_types' in str_path][0] 
        
        
    @staticmethod
    def calc_stimparams(time, stimulus_trace,trace_name):
        """Calculate stimuls start, stop and amplitude from trace"""
    
        nonzero_indices = np.where(stimulus_trace != 0)[0]
    
        # make sure if the stimulus is ok if there was no input
        # if the stimulus is zero
        if not nonzero_indices.any():   # if the list is empty
            # arbitrary values for the no-stimulus response
            stim_start = time[20000]    # after 100ms
            stim_stop = time[-1]        # until the end
            stim_amp_start = 0
            stim_amp_end = 0
            hold_curr = 0
        else:
            # if the stimulus is not zero
            stim_start = time[nonzero_indices[0]]
            stim_stop = time[nonzero_indices[-1]]
            if 'DC' in trace_name:
                hold_curr = np.mean(stimulus_trace[nonzero_indices[-1]+1000:\
                                                   nonzero_indices[-1] + 20000])*1e12
            stim_amp_start = stimulus_trace[nonzero_indices[0]] * 1e12 - hold_curr
            stim_amp_end = stimulus_trace[nonzero_indices[-1]] * 1e12 - hold_curr
            
        tot_duration = time[-1]    
        return stim_start, stim_stop, stim_amp_start, stim_amp_end, tot_duration, hold_curr
    
    @staticmethod
    def calc_stimparams_nonstandard(time, stimulus_trace,trace_name):
        """Calculate stimuls start, stop and amplitude from trace for nonstandard nwb"""
    
        # if the stimulus is not empty
        # find the max/min of the noisy signal
        gradient_thresh = 10  # arbitrary
        gradient_f = np.gradient(stimulus_trace)*1e12
        gradient_f[abs(gradient_f) <= gradient_thresh] = 0
        
        nonzero_indices = np.where(gradient_f != 0)[0]
    
        if not nonzero_indices.any():
            
            stim_start = time[20000]    # after 100ms (arbitrary)
            stim_stop = time[40000]     # after 200ms (arbitrary)
            stim_amp_start = 0.0
            stim_amp_end = 0.0
            hold_curr =  np.mean(stimulus_trace[-20000:])*1e12
       
        else:
            
            signal_max = max(gradient_f)
            signal_min = min(gradient_f)
        
            # find the max/min of the gradient
            first_ind = np.where(gradient_f == signal_max)[0][0]
            second_ind = np.where(gradient_f == signal_min)[0][0]
        
            # check for the first and second indexes
            if first_ind > second_ind:
                start_ind = second_ind
                end_ind = first_ind
            elif first_ind < second_ind:
                start_ind = first_ind
                end_ind = second_ind
        
            stim_start = time[start_ind]
            stim_stop = time[end_ind]
        
            # check for the middle part of the signal
        
            # approximate the amp, it is the mean between the start and end
            if 'DC' in trace_name:
                hold_curr = np.mean(stimulus_trace[end_ind+1000:end_ind + 20000])*1e12
            else:
                hold_curr = 0
            stim_amp = np.mean(stimulus_trace[start_ind:end_ind] ) * 1e12 - hold_curr
            stim_amp_start=stim_amp
            stim_amp_end=stim_amp
        tot_duration = time[-1]
    
        return stim_start, stim_stop, stim_amp_start, stim_amp_end, tot_duration, hold_curr


    @staticmethod
    def write_stimmap_csv(stim_map, output_dir, stim_sweep_map):
        """Write StimMap.csv"""
    
        stim_reps_sweep_map = {}
    
        stimmapreps_csv_content = "DistinctID, StimType, HoldingCurrent, "\
            "Amplitude_Start, Amplitude_End, Stim_Start, Stim_End, Duration, DataPath\n"
    
        reps = defaultdict(lambda: defaultdict(list))
        for stim_type in stim_map:
            for trace_params in stim_map[stim_type]:
    
    
                amplitude = str(trace_params[3])+'&'+ str(trace_params[6])
                reps[stim_type][amplitude].append(trace_params)
    
        for stim_type in reps:
            for amplitude in reps[stim_type]:
    
                cumul_params = reps[stim_type][amplitude][0]
    
                trace_name = cumul_params[0]
    
                cumul_params[2] = np.mean(
                    [rep_params[2] for rep_params in reps
                     [stim_type][amplitude]])
    
                cumul_params[8] = "|".join(
                    rep_params[8] for rep_params in reps[stim_type][amplitude])
    
                rep_names = [rep_params[0]
                             for rep_params in reps[stim_type][amplitude]]
                rep_sweeps = [stim_sweep_map[rep_name] for rep_name in rep_names]
                stim_reps_sweep_map[trace_name] = rep_sweeps
    
                tstart_set = set(['.1f' % rep_params[5]
                                  for rep_params in reps[stim_type][amplitude]])
                if len(tstart_set) != 1:
                    raise Exception(
                        "Stim type %s Amplitude %s don't have equal start "
                        "times: %s" %
                        (stim_type, amplitude.split('&')[0], str(tstart_set)))
    
                tstop_set = set(['.1f' % rep_params[6]
                                 for rep_params in reps[stim_type][amplitude]])
                if len(tstop_set) != 1:
                    raise Exception(
                        "Stim type %s Amplitude %s don't have equal stop "
                        "times: %s" %
                        (stim_type, amplitude.split('&')[0], str(tstop_set)))
    
                stimmapreps_csv_content += ",".join([str(x) for x in cumul_params])
                stimmapreps_csv_content += '\n'
    
        stimmap_filename = 'StimMapReps.csv'
        stimmapreps_csv_filename = os.path.join(output_dir,stimmap_filename)
    
        with open(stimmapreps_csv_filename, 'w') as stimmapreps_csv_file:
            stimmapreps_csv_file.write(stimmapreps_csv_content)
    
        return stim_reps_sweep_map,stimmap_filename
    
    
    @staticmethod
    def calculate_md5hash(filename):
        """Calculate the md5hash of a file"""
    
        import hashlib
        with open(filename, 'rb') as file_h:
            md5hash = hashlib.md5(file_h.read()).hexdigest()
    
        return md5hash
    
    
    def write_provenance(self,
            output_dir,
            nwb_filename,
            stim_sweep_map,
            stim_reps_sweep_map):
        """Writing provenance file"""
    
        provenance_filename = os.path.join(output_dir, 'provenance.json')
    
        nwb_md5hash = self.calculate_md5hash(nwb_filename)
    
        provenance = {
            'nwb_filename': os.path.abspath(nwb_filename),
            'nwb_md5hash': nwb_md5hash,
            'temperature': self.temperature,
            'junction_potential': self.junction_potential,
            'stim_sweep_map': stim_sweep_map,
            'stim_reps_sweep_map': stim_reps_sweep_map}
    
        with open(provenance_filename, 'w') as provenance_file:
            json.dump(
                provenance,
                provenance_file,
                sort_keys=True,
                indent=4,
                separators=(
                    ',',
                    ': '))
    
           
    def save_cell_data(self,acceptable_stimtypes,non_standard_nwb = False):

        bpopt_stimtype_map = utility.bpopt_stimtype_map
        distinct_id_map = utility.aibs_stimname_map
        nwb_file = NwbDataSet(self.nwb_path)
        
        stim_map = defaultdict(list)
        stim_sweep_map = {}
        output_dir = os.getcwd() +'/preprocessed'
        utility.create_dirpath(output_dir)
        
        for sweep_number in nwb_file.get_sweep_numbers():
            sweep_data = nwb_file.get_sweep_metadata(sweep_number)
            stim_type = sweep_data['aibs_stimulus_name']
            
            try:
                stim_type = stim_type.decode('UTF-8')
            except:
                pass
                
            if stim_type in acceptable_stimtypes:
                sweep = nwb_file.get_sweep(sweep_number)
    
                start_idx, stop_idx = sweep['index_range']
    
                stimulus_trace = sweep['stimulus'][start_idx:stop_idx]
                response_trace = sweep['response'][start_idx:stop_idx]
    
                sampling_rate = sweep['sampling_rate']
    
                time = np.arange(0, len(stimulus_trace)) / sampling_rate
                trace_name = '%s_%d' % (
                distinct_id_map[stim_type], sweep_number)
                
                if non_standard_nwb:
                    calc_stimparams_func = self.calc_stimparams_nonstandard
                else:
                    calc_stimparams_func = self.calc_stimparams
                
                stim_start, stim_stop, stim_amp_start, stim_amp_end, tot_duration,hold_curr = calc_stimparams_func(
                    time, stimulus_trace,trace_name)
             
    
                response_trace_short_filename = '%s.%s' % (trace_name, 'txt')
    
                response_trace_filename = os.path.join(
                    output_dir, response_trace_short_filename)
                
    
                
                time *= 1e3 # in ms
                response_trace *= 1e3 # in mV 
                response_trace = utility.correct_junction_potential(response_trace,
                                                            self.junction_potential)
                
                # downsampling
                time,response_trace = utility.downsample_ephys_data(time,response_trace)
                    
                with open(response_trace_filename, 'wb') as response_trace_file:
                    np.savetxt(response_trace_file,
                                  np.transpose([time, response_trace]))
                
                holding_current = hold_curr  # sweep['bias_current']
    
                stim_map[distinct_id_map[stim_type]].append([
                    trace_name,
                    bpopt_stimtype_map[stim_type],
                    holding_current/1e12,
                    stim_amp_start /1e12,
                    stim_amp_end/1e12,
                    stim_start * 1e3,
                    stim_stop * 1e3,
                    tot_duration * 1e3,
                    response_trace_short_filename])
    
                stim_sweep_map[trace_name] = sweep_number
                
        logger.debug('Writing stimmap.csv ...')
    
        stim_reps_sweep_map,stimmap_filename = self.write_stimmap_csv(stim_map, output_dir, stim_sweep_map)
        
        self.write_provenance(
            output_dir,
            self.nwb_path,
            stim_sweep_map,
            stim_reps_sweep_map)
            
        return output_dir,stimmap_filename    
    
    
    @staticmethod
    def get_stim_map(stim_map_filename, dend_recording = None, locations = None):
        """Get stim map"""
    
        stim_map = defaultdict(dict)
    
        with open(stim_map_filename, 'r') as stim_map_file:
            stim_map_content = stim_map_file.read()
    
        for line in stim_map_content.split('\n')[1:-1]:
            if line != '':
                stim_name, stim_type, holding_current, amplitude_start, amplitude_end, \
                    stim_start, stim_end, duration, sweeps = line.split(',')
                iter_dict1, iter_dict2 = dict(), dict()
                iter_dict1['type'] = stim_type.strip()
                iter_dict1['amp'] = 1e9 * float(amplitude_start)
                iter_dict1['amp_end'] = 1e9 * float(amplitude_end)
                iter_dict1['delay'] = float(stim_start)
                iter_dict1['duration'] = float(stim_end) - float(stim_start)
                iter_dict1['stim_end'] = float(stim_end)
                iter_dict1['totduration'] = float(duration)
                iter_dict1['sweep_filenames'] = [
                    x.strip() for x in sweeps.split('|')]
                
                if 'Ramp' in stim_name:
                    holding_current = 0
                iter_dict2['type'] = 'SquarePulse'
                iter_dict2['amp'] = 1e9 * float(holding_current)
                iter_dict2['amp_end'] = 1e9 * float(holding_current)
                iter_dict2['delay'] = 0
                iter_dict2['duration'] = float(duration)
                iter_dict2['stim_end'] = float(duration)
                iter_dict2['totduration'] = float(duration)
                
                if float(holding_current) != 0.0:
                    iter_list = [iter_dict1, iter_dict2]
                else:
                    iter_list = [iter_dict1]
                
                stim_map[stim_name]['stimuli'] = iter_list
                if dend_recording:
                    record_list = list()
                    for i, loc in enumerate(locations):
                        record_dict = dict()
                        record_dict['var'] = 'v'
                        record_dict['somadistance'] = loc
                        record_dict["seclist_name"] = "apical"
                        record_dict['name'] = 'dend'+ str(i+1)
                        record_dict['type'] = 'somadistance'
                        record_list.append(record_dict)
                    stim_map[stim_name]['extra_recordings'] = record_list
        return stim_map
    
    
    def get_ephys_features(self,feature_set_filename,ephys_data_path,stimmap_filename, 
                           filter_rule_func,*args,
                           dend_recording = None, record_locations = None):
        
        cell_name = self.cell_id
        features_write_path = 'config/'+ cell_name +'/features.json'
        untrained_features_write_path = 'config/'+ cell_name +'/untrained_features.json'
        protocols_write_path = 'config/'+cell_name+'/protocols.json'
        all_protocols_write_path = 'config/'+cell_name+'/all_protocols.json'
        
        utility.create_filepath(all_protocols_write_path)
        feature_file = utility.locate_template_file(feature_set_filename)
        feature_map = utility.load_json(feature_file)
        stim_features = feature_map['features'] # Features to extract
        
        features_meanstd = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(dict)))


        stim_map = self.get_stim_map(os.path.join(ephys_data_path,stimmap_filename),
                            dend_recording = dend_recording, locations = record_locations)

        cell_stim_map= stim_map.copy()
        training_stim_map = dict()
        
        for stim_name in stim_map.keys():
            
            logger.debug("\n### Getting features from %s of cell %s ###\n" \
                % (stim_name, cell_name))

            sweeps = []
            for sweep_filename in stim_map[stim_name]['stimuli'][0]['sweep_filenames']:
                sweep_fullpath = os.path.join(
                    ephys_data_path,
                    sweep_filename)

                data = np.loadtxt(sweep_fullpath)
                time = data[:, 0]
                voltage = data[:, 1]


                # Prepare sweep for eFEL
                sweep = {}
                sweep['T'] = time
                sweep['V'] = voltage
                sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
                sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
                sweep['T;location_AIS'] = time
                sweep['V;location_AIS'] = voltage
                sweep['stim_start;location_AIS'] = [stim_map[stim_name]['stimuli'][0]['delay']]
                sweep['stim_end;location_AIS'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
                sweeps.append(sweep)

            # Do the actual feature extraction
            feature_results = efel.getFeatureValues(
                sweeps, stim_features)
            
            for feature_name in stim_features:
                # For one feature, a list with values for every sweep
                feature_values = [np.mean(trace_dict[feature_name])
                                  for trace_dict in feature_results
                                  if trace_dict[feature_name] is not None]
                if len(feature_values) == 0:
                    continue
                elif len(feature_values) == 1:
                    mean = feature_values[0]
                    std = 0.05 * abs(mean)
                elif len(feature_values) > 1:
                    mean = np.mean(feature_values)
                    std = np.std(feature_values)
                
                if std== 0 and len(feature_values) != 1:
                    std = 0.05 * abs(mean)/math.sqrt(len(feature_values)) 
                
                if math.isnan(mean) or math.isnan(std):
                    continue
                if mean == 0:
                    std = 0.05
                
                if feature_name in ['voltage_base', 'steady_state_voltage'] \
                        and len(feature_values) == 1:
                    std = 0
                     
                
                features_meanstd[stim_name]['soma'][
                    feature_name] = [mean , std]
            if stim_name in features_meanstd.keys():
                training_stim_map[stim_name] = cell_stim_map[stim_name]
        

        features_meanstd_filtered,untrained_features_dict,training_stim_map_filtered,\
                all_stim_filtered = filter_rule_func(features_meanstd,training_stim_map,cell_stim_map,*args)    
        utility.save_json(features_write_path,features_meanstd_filtered)
        utility.save_json(untrained_features_write_path,untrained_features_dict)
        utility.save_json(protocols_write_path,training_stim_map_filtered)
        utility.save_json(all_protocols_write_path,all_stim_filtered)
        
        return features_write_path,untrained_features_write_path,\
                                protocols_write_path,all_protocols_write_path
    
    
