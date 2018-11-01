#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 15:25:03 2018

@author: anin
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from sklearn.decomposition import PCA
import numpy.linalg as la
from allensdk.core.cell_types_cache import CellTypesCache
import allensdk.core.swc as swc
from pyquaternion import Quaternion
import glob
import re
import math

import animate3d

#topdir = '.'
#dir_list = list()
#def step(ext, dirname, names):
#    ext = ext.lower()
#    for name in names:
#        if name.lower().endswith(ext):
#            dir_list.append(os.path.join(dirname, name))
            
def get_morph_path(exten = '.swc'):
    morph_path = glob.glob('./morphology/*'+exten)
    return morph_path[0]           



def shift_origin(coord_tuple, morphology):
    
    (x,y,z) = coord_tuple
    x -= morphology.soma["x"]
    y -= morphology.soma["y"]
    z -= morphology.soma["z"]
    return (x,y,z)
   
def get_cell_morphXYZ(morph_path):
    
    morphology = swc.read_swc(morph_path) 
    
    x, x_apical =[], []
    y, y_apical =[], []
    z, z_apical =[], []

    for n in morphology.compartment_list:
        coord_tuple = (n['x'],n['y'],n['z'])
        x_coord, y_coord, z_coord =  shift_origin(coord_tuple,morphology)
        x.append(x_coord)
        y.append(y_coord)
        z.append(z_coord)
        
        if n['type']==2 :
            x_apical.append(x_coord)
            y_apical.append(y_coord)
            z_apical.append(z_coord)

    morph_data = np.array(np.column_stack((x,y,z))) 
    morph_apical = np.array(np.column_stack((x_apical,y_apical,z_apical))) 
    morph_soma = [0,0,0]
    
    return morph_data,morph_soma,morph_apical



def cal_rotation_angle(morph_apical):
    
    pca = PCA(n_components=2)
    pca.fit(morph_apical)
    v1 = pca.components_[0]  # the first principal component 
    z_axis = np.array([0, 0, 1])
    v1_unit = v1/la.norm(v1)
    theta = np.arccos(np.clip(np.dot(z_axis, v1_unit), -1.0, 1.0))
    axis_of_rot = np.cross(z_axis,v1_unit)
    return theta, axis_of_rot    

def rotate3D_point(point,theta,axis_of_rot):
    
    # use one of these two commands to rotate the morphology upside down
    
    point_rotated = Quaternion(axis=axis_of_rot,angle=-theta ).rotate(point)
#    point_rotated = Quaternion(axis=axis_of_rot,angle=-theta + math.pi).rotate(point)

    return point_rotated





def Main():

    morph_path = get_morph_path()
    morph_data,morph_soma,morph_apical = get_cell_morphXYZ(morph_path)
    theta, axis_of_rot = cal_rotation_angle(morph_apical)    
    color_dict = {4:'darkred',3:'orange',2:'royalblue',1:'black'}   
    label_dict = {4:'apical dendrite',3:'basal dendrite',2:'axon',1:'soma'}
    
    plt.style.use('classic')
    fig = plt.figure(figsize=(10, 8), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    
    morphology = swc.read_swc(morph_path) 
    
    all_x, all_y, all_z = [], [], []
    
    for n in morphology.compartment_list:
        for c in morphology.children_of(n):
            nx,ny,nz = shift_origin((n['x'],n['y'],n['z']),morphology)
            [nx_rot,ny_rot,nz_rot] = rotate3D_point([nx,ny,nz],theta,axis_of_rot)
            
            cx,cy,cz = shift_origin((c['x'],c['y'],c['z']),morphology)
            [cx_rot,cy_rot,cz_rot] = rotate3D_point([cx,cy,cz],theta,axis_of_rot)
            
            ax.plot([nx_rot, cx_rot], [ny_rot,cy_rot], [nz_rot,cz_rot], 
                    color=color_dict[n['type']],lw =1,alpha =.5,label = label_dict[n['type']])
            
            all_x.extend((nx_rot, cx_rot))
            all_y.extend((ny_rot, cy_rot))
            all_z.extend((nz_rot, cz_rot))
            
    ax.scatter(morph_soma[0], morph_soma[1],morph_soma[2],marker = 'o',
                           s=150,color= 'k', alpha =.6)
#    ax.legend()
    all_x_min, all_x_max = min(all_x), max(all_x)
    all_y_min, all_y_max = min(all_y), max(all_y)
    all_z_min, all_z_max = min(all_z), max(all_z)


    ax.set_xlim([all_x_min-10, all_x_max+10])
    ax.set_ylim([all_y_min-10, all_y_max+10])
    ax.set_zlim([all_z_min-10, all_z_max+10])
#    
#
    ax.axis('off')
    plt.show()
    
    angles = np.linspace(0,360,21)[:-1] # Take 20 angles between 0 and 360
 
#     create an animated gif (20ms between frames)
    cell_name = re.search('morphology/(.*).swc',morph_path).group(1)
    movie_name = cell_name+'_movie.gif'
    animate3d.rotanimate(ax, angles,movie_name,delay=20) 
    
if __name__== '__main__':
    Main()