from scipy import interpolate
import numpy as np


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