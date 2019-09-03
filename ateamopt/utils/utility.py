import os,errno
import glob
import json
import numpy as np
import pkg_resources
import ateamopt.template as templ
import pickle
import allensdk.core.swc as swc
import ateamopt.scripts as pyscripts
import logging

logger = logging.getLogger(__name__)

bpopt_section_map = {
                      'soma':'somatic',
                      'apic':'apical',
                      'dend':'basal',
                      'axon':'axonal',
                      'all' : 'all'
                    }

bpopt_section_map_inv = {
                    'somatic':'soma', 
                     'axonal':'axon', 
                     'apical':'apic',
                     'basal':'dend', 
                     'all':'all'
                 }

bpopt_stimtype_map = {
                        'Long Square': 'SquarePulse',
                        'Ramp': 'RampPulse',
                        'Square - 2s Suprathreshold': 'SquarePulse',
                        'Ramp to Rheobase': 'RampPulse',
                        'Short Square - Triple' : 'TriBlip',
                        'Noise 1': 'Noise',
                        'Noise 2': 'Noise'
                     }

aibs_stimname_map = {
                        'Long Square': 'LongDC',
                        'Ramp': 'Ramp',
                        'Square - 2s Suprathreshold': 'LongDCSupra',
                        'Short Square - Triple' : 'Short_Square_Triple',
                        'Ramp to Rheobase': 'RampRheo',
                        'Noise 1': 'Noise_1',
                        'Noise 2': 'Noise_2',
                    }

aibs_stimname_map_inv = {
                        'LongDC':'Long Square',
                        'Ramp':'Ramp',
                        'LongDCSupra': 'Square - 2s Suprathreshold',
                        'Short_Square_Triple':'Short Square - Triple',
                        'Ramp to Rheobase': 'RampRheo',
                        'Noise_1':'Noise 1',
                        'Noise_2':'Noise 2'
                    }

bpopt_current_play_stimtypes = ['Short Square - Triple', 'Noise 1', 'Noise 2']

rev_potential = {'ena' : 53, 'ek' : -107}

#passive_params = ['cm', 'Ra', 'g_pas', 'e_pas']

def correct_junction_potential(data,junction_potential):
    return data+junction_potential # Assumes junction potential negative (AIBS specific)

def reverse_junction_potential_correction(data,junction_potential): # For release to web
    return data-junction_potential # Assumes junction potential negative (AIBS specific)

def create_dirpath(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise   

def create_filepath(path):
    if not os.path.exists(path):
        try:
            os.makedirs(os.path.dirname(path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise                   
                
def get_filepath_for_exten(exten, topdir = '.'):
        dir_list = list()
        
        def step(ext, dirname, names):
            ext = ext.lower()
            for name in names:
                if name.lower().endswith(ext):
                    dir_list.append(os.path.join(dirname, name)) 
        
        try:            
            os.path.walk(topdir, step, exten)
        except:
            dir_list = glob.glob(topdir+'/**/*'+exten,recursive=True)
        return dir_list


def save_json(path, content):
    with open(path,'w') as json_write:
        json.dump(content,json_write,indent=4)

def save_file(path,content):
    with open(path, 'a') as handle:
        handle.write(content)

def save_pickle(path, content):
    with open(path,'wb') as pickle_write:
        pickle.dump(content,pickle_write)
        
        
def load_json(path):
    with open(path,'r') as json_read:
        json_data=json.load(json_read)
        
    return json_data

def load_pickle(path):
    try:
        with open(path,'rb') as pickle_read:
            pickle_data=pickle.load(pickle_read)
    except:
        with open(path, 'rb') as f:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            pickle_data = u.load()
        logger.debug('solving python 2-3 pickling incompatibility')
    return pickle_data


def downsample_ephys_data(time,stim,response,downsample_interval=5):
    
    time_end = time[-1]
    stim_end = stim[-1]
    response_end = response[-1]
                
    time = time[::downsample_interval]
    stim = stim[::downsample_interval]
    response = response[::downsample_interval]
    if time_end != time[-1]:
        time = np.append(time,time_end)
        stim = np.append(stim,stim_end)
        response = np.append(response,response_end)
        
    return time,stim,response

def check_swc_for_apical(morph_path):
    morphology = swc.read_swc(morph_path)
    no_apical = True
    for n in morphology.compartment_list:
        if n['type']==4 :
            no_apical=False
            break
    return no_apical

def remove_entries_dict(dict_, entries):
    for k in entries:
        dict_.pop(k, None)
    return dict_

def locate_template_file(rel_file_path):
    try:
        file_path = pkg_resources.resource_filename(templ.__name__, 
                                                    rel_file_path)
    except:
        file_path = None
    return file_path

def locate_script_file(rel_file_path):
    try:
        file_path = pkg_resources.resource_filename(pyscripts.__name__,
                                                    rel_file_path)
    except:
        file_path = None
    return file_path


