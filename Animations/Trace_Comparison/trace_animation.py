#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 13:25:30 2018

@author: anin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
import pickle

matplotlib.rcParams.update({'font.size': 16})
plt.style.use('ggplot')

optimized_resp = pickle.load(open('resp_opt.txt', 'r'))
released_resp = pickle.load(open('resp_release.txt','r'))

select_key = 'LongDC_47'
soma_loc = '.soma.v'
AIS_loc = '.AIS.v'

exp_response_file = 'preprocessed/' + select_key +'.txt'

exp_response = np.loadtxt(exp_response_file)
exp_resp_time = exp_response[:,0]
exp_resp_voltage = exp_response[:,1]

stim_key_somaloc = select_key + soma_loc
stim_key_AISloc = select_key + AIS_loc

select_opt_resp_soma = optimized_resp[0][stim_key_somaloc]
time_soma_opt = select_opt_resp_soma['time']
voltage_soma_opt = select_opt_resp_soma['voltage']

select_opt_resp_AIS = optimized_resp[0][stim_key_AISloc]
time_AIS_opt = select_opt_resp_AIS['time']
voltage_AIS_opt = select_opt_resp_AIS['voltage']

select_released_resp_soma = released_resp[stim_key_somaloc]
time_soma_released = select_released_resp_soma['time']
voltage_soma_released = select_released_resp_soma['voltage']


time_min = 200
time_max = 1320
num_points = 1000
time_points = np.linspace(time_min, time_max, num_points)

fig, ax = plt.subplots(figsize = (8,6),dpi=80)

ax.set_ylim([-100, 50])
color_list = ['k', 'b']
label_list = ['Experiment', 'Model']
lw_list = [1.5,1.5]
lines = []

for i,index in enumerate(range(2)):
    line, = ax.plot([], [], color = color_list[i],lw = lw_list[i],
                    label = label_list[i], alpha = .8)
    lines.append(line)
fig.legend(handles = (lines[0],lines[1]), loc = 'upper center', 
           ncol=2, fontsize = 16)
ax.set_xlabel('Time (ms)', fontsize = 16)
ax.set_ylabel('Voltage (mV)',fontsize = 16)

def init():
    for line in lines:
        line.set_data([],[])
    return lines


def animate(i):
    anim_time = time_points[i]
    data_list_time = [exp_resp_time, time_soma_opt] 
    data_list_voltage = [exp_resp_voltage, voltage_soma_opt]
    
    for i,line in enumerate(lines):
        time = data_list_time[i]
        voltage = data_list_voltage[i]
        xdata = time[time <= anim_time]
        ydata = voltage[time <= anim_time]
        line.set_data(xdata, ydata)
    ax.set_xlim([time_min, anim_time +5])
    return lines


ani = animation.FuncAnimation(fig, animate, init_func=init, frames = num_points, 
                              interval=25, blit=True)
ani.save('test.mp4')
plt.show()


######## Animation to show the AIS initiation of AP ########

num_points = 1500
time_points = np.linspace(time_min, time_max, num_points)
fig, ax = plt.subplots(figsize = (8,6),dpi=80)

ax.set_ylim([-100, 50])
#color_list = ['k', 'b', 'r']
label_list = ['Soma', 'AIS']
lw_list = [1.5,1.5]
lines = []

for i,index in enumerate(range(2)):
    line, = ax.plot([], [], lw = lw_list[i],
                    label = label_list[i], alpha = .8)
    lines.append(line)
fig.legend(handles = (lines[0],lines[1]), loc = 'upper center', 
           ncol=2, fontsize = 16)
ax.set_xlabel('Time (ms)', fontsize = 16)
ax.set_ylabel('Voltage (mV)', fontsize = 16)

def init():
    for line in lines:
        line.set_data([],[])
    return lines


def animate(i):
    anim_time = time_points[i]
    data_list_time = [time_soma_opt, time_AIS_opt] 
    data_list_voltage = [voltage_soma_opt, voltage_AIS_opt]
    
    for i,line in enumerate(lines):
        time = data_list_time[i]
        voltage = data_list_voltage[i]
        xdata = time[time <= anim_time]
        ydata = voltage[time <= anim_time]
        line.set_data(xdata, ydata)
    ax.set_xlim([time_min, anim_time +5])
    return lines


ani = animation.FuncAnimation(fig, animate, init_func=init, frames = num_points, 
                              interval=25, blit=True)
ani.save('test_AIS.mp4')
plt.show()
