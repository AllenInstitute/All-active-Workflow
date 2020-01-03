from scipy import interpolate
import numpy as np
import scipy.signal as signal
import pandas as pd
import os

def get_spike_shape(time,voltage,spike_times,
                    AP_shape_time, AP_shape_voltage):
    
    prefix_pad = 10.5 #should cover 2ms (default)
    postfix_pad = 15.5 #should cover 5ms
    
    for i,spike_time in enumerate(spike_times):
        min_index = np.argmax(time >= spike_time - prefix_pad) 
        max_index = np.argmax(time >= spike_time + postfix_pad) 
        AP_shape = voltage[min_index:max_index]
        t_shape  = time[min_index:max_index] - spike_time
        f_shape = interpolate.interp1d(t_shape, AP_shape)
        AP_shape_voltage += f_shape(AP_shape_time)- \
                f_shape(AP_shape_time[0])
    return AP_shape_voltage


def calculate_spike_time_metrics(expt_trains, model_train, total_length, dt, sigma):
    """
    Calculate explained variance (%) in terms of  the experiment and model spike train.
    
    Parameters
    ----------
    expt_trains : list
        list of nparrays, where each array contains the time indices of spike times (from all experimental trials)
    model_train :  np.array
        is a single nparray which contains the time indices of spike times from the model
    total_length : int
        the total number of time bins
    dt : float 
        Time step between each bin
    sigma : float
        Length of the Gaussian convolution window (same units as dt)
        
    Returns
    -------
    float
        Explained variance in percentage
        
    """
    trialtotrial_expvar_result = []
    for s in sigma:
        sigma_points = s / dt
        window_length = 10 * sigma_points
        gauss_func = signal.gaussian(window_length, sigma_points)
        binary_expt_trains = []
        for et in expt_trains:
            bt = np.zeros(total_length)
            if len(et) > 0:
                bt[et] = 1
            binary_expt_trains.append(bt)
        convolved_expt_trains = [signal.fftconvolve(gauss_func, bt_) for bt_ in binary_expt_trains]
        avg_convolved_expt_train = np.mean(convolved_expt_trains, axis=0)
        binary_model_train = np.zeros(total_length)
        if len(model_train) > 0:
            binary_model_train[model_train] = 1
        convolved_model_train = signal.fftconvolve(gauss_func, binary_model_train)

        f = np.array(convolved_expt_trains)
        trialtotrial_expvar_result.append(trial_expvar(f, convolved_model_train) / trial_expvar(f, avg_convolved_expt_train))

    return np.array(trialtotrial_expvar_result)


def trial_expvar(f, m):
    if f.shape[1] != len(m):
        print("Comparison trace has a different number of points than trial traces")
        return np.nan

    fbar = np.mean(f)
    mbar = np.mean(m)
    var_sum = np.sum((f - fbar) ** 2) + f.shape[0] * np.sum((m - mbar) ** 2)
    ev = (var_sum - (np.sum(((f - fbar) - (m - mbar)) ** 2))) / var_sum
    return ev


def save_optimization_time(time_by_gen_filename,time_metrics_filename,
                           cell_metadata):
    
    def get_sec(time_str):
        h, m, s = time_str.split(':')
        return float(h) * 3600 + float(m) * 60 + float(s)
    
    total = 0
    counter = 0
    time_vec = []
    with open(time_by_gen_filename, 'r') as inp, open('total_time.txt', 'w') as outp: 
       for line in inp:
           try:
               num = get_sec(line.strip())
               total += num
               counter += 1
               time_vec.append(num)
           except ValueError:
               print('{} is not a number!'.format(line))
               
       outp.write('#Generations %s: %s seconds'%(counter,total))
    print('{} Generations took: {} seconds'.format(counter,total))
    
    cell_id = cell_metadata['cell_id']
    time_metrics = pd.DataFrame({
                        'cell_id' : [cell_id for i in range(len(time_vec))],
                        'time' : time_vec
                        })
    time_metrics.to_csv(time_metrics_filename) 
    
    
def save_compute_statistics(opt_logbook,output_csv_filename):
    """
    Saving generation wise compute time from the optimization logbook
    
    
    Parameters
    ----------
    opt_logbook : str
        location of the optimization logbook 
    output_csv_filename : str
        location of the output csv file
    
    """
    
    if os.path.exists(opt_logbook):
        colnames = ['gen','nevals','avg','std','min','max','cp_loc']
        col_dtypes = {'gen':int,'nevals':int,'avg': float, 'std': float, 'min': float,
                      'max': float, 'cp_loc':str}
        logbook = pd.read_fwf(opt_logbook,names=colnames,
                              header=0)
        logbook_filtered = logbook[logbook['gen'] != 'gen']
        nan_seed_indices = logbook_filtered.index[logbook_filtered.cp_loc.isnull()].tolist()
        actual_seed_indices = [idx+1 for idx in nan_seed_indices]
        nan_cp_locations = logbook_filtered['cp_loc'][actual_seed_indices].tolist()
        logbook_filtered.loc[logbook_filtered.cp_loc.isnull(),'cp_loc']=nan_cp_locations
        logbook_filtered = logbook_filtered.astype(dtype=col_dtypes).reset_index()
        
        logbook_filtered['seed_index'] = logbook_filtered['cp_loc'].apply(lambda x:x.split('seed')[-1])
        logbook_filtered.to_csv(output_csv_filename)
    