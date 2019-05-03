import matplotlib
matplotlib.use('Agg')
from allensdk.core.cell_types_cache import CellTypesCache

import pickle
import matplotlib.pyplot as plt
import os
import numpy as np
from functools import partial
from ipyparallel import Client
import pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def plot_fi_data(data, title = 'fi_curve'):
    fig,ax = plt.subplots(1, dpi = 80, figsize = (6,5))
    plt.style.use('ggplot')
    spike_onset_list = []
    for data_ in data:
        data_stim = data_['stim_list']
        data_spikes = data_['spike_list']
        data_spike_ordered = [x for _,x in sorted(zip(data_stim,data_spikes))]
        data_stim_ordered = sorted(data_stim)
        spike_onset_list.append(get_fi_intercept(data_stim_ordered,data_spike_ordered))
        ax.plot(data_stim_ordered, data_spike_ordered, alpha = 0.7, color = 'k')
    
    ax.set_xlabel('Stim amplitude (in pA)')
    ax.set_ylabel('Spike Rate')

    ax.set_xlim(-10, 360)
    ax.set_ylim(-5, 25)
    
    fig_title = title + ' (Avg. Spike onset stim = %s pA)'%np.mean(spike_onset_list)
    ax.set_title(fig_title)
    fig.savefig('fi_curve.pdf', bbox_inches = 'tight')
    plt.close(fig)
        
def get_fi_intercept(stim_fi, spike_fi):
    spiking_stim_idx = [ix for ix,rate in enumerate(spike_fi) if rate >0]
    spiking_stim = [stim for idx,stim in enumerate(stim_fi) if idx in spiking_stim_idx \
                    and stim > 0]
    intercept = spiking_stim[0]
    return intercept      

def get_fi_data(ctc, cell_id):
    
    from collections import defaultdict
    from allensdk.ephys.extract_cell_features import extract_cell_features
    import shutil

    data_set = ctc.get_ephys_data(cell_id)
    sweeps = ctc.get_ephys_sweeps(cell_id)

    # group the sweeps by stimulus 
    sweep_numbers = defaultdict(list)
    for sweep in sweeps:
        sweep_numbers[sweep['stimulus_name']].append(sweep['sweep_number'])
    
    try:
        cell_features = extract_cell_features(data_set,
                                                      sweep_numbers['Ramp'],
                                                      sweep_numbers['Short Square'],
                                                      sweep_numbers['Long Square'])
                
        all_stim_sweeps = cell_features['long_squares']['sweeps']
        stim_list = [stim_sweep['stim_amp'] for stim_sweep in all_stim_sweeps]
        spike_list = [stim_sweep['avg_rate'] for stim_sweep in all_stim_sweeps]
        fi_dict = {'id' : cell_id,
                         'stim_list' : stim_list,
                         'spike_list' : spike_list
                         }
    except:
        fi_dict = None
    shutil.rmtree("cell_types/specimen_%s"%cell_id)    
    return fi_dict

def Main():
    filter_obj = {'dendrite_type':'spiny', 'structure_layer_name' : '5', 
                   'structure_area_abbrev' : 'VISp'}
    
    ctc = CellTypesCache()
    cells = ctc.get_cells(species = ['Mus musculus'])
    
    cells_df = pd.DataFrame(cells)
    for filt_key,filt_val in filter_obj.items():
        cells_df = cells_df.loc[cells_df[filt_key] == filt_val,:]
    
    cell_ids = list(cells_df['id'].values)
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
    logger.debug('Using ipyparallel with %d engines', len(rc))
    lview = rc.load_balanced_view()
    
    func = partial(get_fi_data, ctc)
    filter_fi_data = lview.map_sync(func, cell_ids)
    filter_fi_data = [data for data in filter_fi_data if data is not None]
    file_name = 'fi_data.pkl'
    with open(file_name, 'wb') as fh:
        pickle.dump(filter_fi_data, fh)
    plot_fi_data(filter_fi_data)
    
    rc.shutdown(hub=True)

if __name__ == '__main__':
    Main()