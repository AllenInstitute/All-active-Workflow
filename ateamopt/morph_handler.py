import numpy as np
import matplotlib.pyplot as plt
import allensdk.core.swc as swc
import numpy.linalg as la
from sklearn.pipeline import Pipeline
import math
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation
from mpl_toolkits.mplot3d import Axes3D
import os
import neurom as nm
from neurom import morphmath as mm
from neurom.core.types import tree_type_checker, NEURITES
from ateamopt.utils import utility
import seaborn as sns
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection

swc_dict = {4: 'apical dendrite', 3: 'basal dendrite',
                      2: 'axon', 1: 'soma'}

dark_grey = (128/255,128/255,128/255)

class MorphHandler(object):
    def __init__(self, morph_file, cell_id=None):
        self.morph_file = morph_file
        self.morph = swc.read_swc(morph_file)
        self.soma_coord = np.array([self.morph.soma["x"],
                                    self.morph.soma["y"],
                                    self.morph.soma["z"]
                                    ])
        self.soma_rad = self.morph.soma['radius']
        self.cell_id = cell_id

    def save_morph_data(self, morph_stats_filename):

        if not os.path.exists(morph_stats_filename):

            #  load a neuron from an SWC file
            nrn = nm.load_neuron(self.morph_file)

            morph_stats = {}

            morph_stats['cell_id'] = self.cell_id
            morph_stats['soma_suface'] = nm.get('soma_surface_areas', nrn)[0]
            morph_stats['soma_radius'] = np.mean(nm.get('soma_radii', nrn))

            # Morph stats
            for nrn_type_ in NEURITES:

                morph_stats['length' + '.' + str(nrn_type_).split('.')[1]] = np.sum(nm.get('segment_lengths',
                                                                                           nrn, neurite_type=nrn_type_))
                morph_stats['area' + '.' + str(nrn_type_).split('.')[1]] = sum(mm.segment_area(s)
                                                                               for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(nrn_type_)))
                morph_stats['volume' + '.' + str(nrn_type_).split('.')[1]] = sum(mm.segment_volume(s)
                                                                                 for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(nrn_type_)))
                morph_stats['taper_rate' + '.' + str(nrn_type_).split('.')[1]] = \
                    np.mean([mm.segment_taper_rate(s) for s in nm.iter_segments(
                        nrn, neurite_filter=tree_type_checker(nrn_type_))])

            utility.save_json(morph_stats_filename, morph_stats)

    def get_morph_coords(self, reject_axon=True):
        x, x_apical, x_axon = [], [], []
        y, y_apical, y_axon = [], [], []
        z, z_apical, z_axon = [], [], []
        
        morph_dist_arr = []
        
        for comp_ in self.morph.compartment_list:
            x_coord, y_coord, z_coord = comp_['x'], comp_['y'], comp_['z']
            if comp_['type'] == 2:
                x_axon.append(x_coord)
                y_axon.append(y_coord)
                z_axon.append(z_coord)

            if reject_axon and comp_['type'] == 2:
                continue

            x.append(x_coord)
            y.append(y_coord)
            z.append(z_coord)
            
            dist = np.linalg.norm(self.soma_coord-np.array([x_coord,y_coord,z_coord]))
            morph_dist_arr.append(dist)
            
            if comp_['type'] == 4:
                x_apical.append(x_coord)
                y_apical.append(y_coord)
                z_apical.append(z_coord)

        morph_data = self.shift_origin(np.array(np.column_stack((x, y, z))))
        morph_apical = self.shift_origin(np.array(np.column_stack((x_apical,
                                                                   y_apical, z_apical))))
        morph_axon = self.shift_origin(np.array(np.column_stack((x_axon,
                                                                 y_axon, z_axon))))

        return morph_data, morph_apical, morph_axon,morph_dist_arr

    def shift_origin(self, coord_arr):
        shifted_coord = coord_arr - self.soma_coord
        return shifted_coord

    def calc_rotation_angle(self, morph_data, morph_apical=None):
        pca_ = PCA(n_components=2)
        pca_pipeline = Pipeline([('pca', pca_)])

        pca_pipeline.fit(morph_data)
        v1 = pca_.components_[0]  # the first principal component

        z_axis = np.array([0, 0, 1])  # target rotation direction
        v1_unit = v1/la.norm(v1)
        v1_sign = np.sign(v1_unit[2])
        v1_unit *= v1_sign
        theta = np.arccos(np.clip(np.dot(z_axis, v1_unit), -1.0, 1.0))
        axis_of_rot = np.cross(z_axis, v1_unit)
        axis_of_rot = axis_of_rot/np.linalg.norm(axis_of_rot)
        try:
            proj_dir = np.sign(np.mean(morph_apical.dot(v1_unit)))
        except:
            proj_dir = np.sign(np.mean(morph_data.dot(v1_unit)))
        if proj_dir == 1:
            theta = 2*math.pi-theta
        elif proj_dir == -1:
            theta = -theta + math.pi

        return theta, axis_of_rot

    def calc_euler_angle(self, morph_data, morph_apical=None):
        pca_ = PCA(n_components=2)
        pca_pipeline = Pipeline([('pca', pca_)])

        pca_pipeline.fit(morph_data)
        v1 = pca_.components_[0]  # the first principal component

        y_axis = np.array([0, 1, 0])  # target rotation direction
        v1_unit = v1/la.norm(v1)
        v1_sign = np.sign(v1_unit[1])
        v1_unit *= v1_sign
        theta = np.arccos(np.clip(np.dot(y_axis, v1_unit), -1.0, 1.0))
        axis_of_rot = np.cross(y_axis, v1_unit)
        axis_of_rot = axis_of_rot/np.linalg.norm(axis_of_rot)
        try:
            proj_dir = np.sign(np.mean(morph_apical.dot(v1_unit)))
        except:
            proj_dir = np.sign(np.mean(morph_data.dot(v1_unit)))
        if proj_dir == 1:
            theta = 2*math.pi-theta
        elif proj_dir == -1:
            theta = -theta + math.pi

        r = Rotation.from_rotvec(theta*axis_of_rot)
        r_euler = r.as_euler('xyz')
        return r_euler

    @staticmethod
    def rotate3D_point(point, theta, axis_of_rot):
        r = Rotation.from_rotvec(theta*axis_of_rot)
        point_rotated = r.apply(point)
        return point_rotated

    def draw_sphere(self, center_tuple):
        rad = self.soma_rad
        xCenter, yCenter, zCenter = center_tuple

        # draw sphere
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones(np.size(u)), np.cos(v))
        # shift and scale sphere
        x = rad*x + xCenter
        y = rad*y + yCenter
        z = rad*z + zCenter
        return (x, y, z)


    def add_synapses(self,section_points,n_syn,theta,axis_of_rot,ax,**kwargs):
        color = 'lime' if not kwargs.get('color') else kwargs.get('color')
        np.random.seed(0)
        select_syn_loc = section_points[np.random.choice(section_points.shape[0],n_syn,
                                                         replace=False),:]
        select_syn_loc = np.apply_along_axis(lambda x:self.rotate3D_point(x,
                                               theta, axis_of_rot),1,select_syn_loc)
        ax.scatter(select_syn_loc[:,0], select_syn_loc[:,1],
                       select_syn_loc[:,2],
                       marker='o', alpha=1,linewidth=0,
                       s=10,color=color)
        return ax


    def draw_morphology_2D(self, theta, axis_of_rot, reject_axon=True,
                           soma_loc= np.array([0,0]),**kwargs):
        color_dict = kwargs.get('color_dict')
        if not color_dict:
            color_dict = {4: 'orange', 3: 'darkred', 2: 'royalblue', 1: dark_grey}


        morph_dist_arr = kwargs.get('morph_dist_arr')
        if morph_dist_arr:
            _,max_dist = np.min(morph_dist_arr),np.max(morph_dist_arr)
            lw_min = kwargs.get('lw_min') or .2
        lw_max = kwargs.get('lw') or 1
        ax = kwargs.get('ax')
        alpha = kwargs.get('alpha') or 1
        
        if not ax:
            sns.set(style='whitegrid')
            fig,ax = plt.subplots() 
        
        all_lines =[]
        linewidths,colors = [],[]
        for comp_ in self.morph.compartment_list:
            if reject_axon and comp_['type'] == 2:
                continue
            nx, ny, nz = self.shift_origin(np.array([comp_['x'], comp_['y'], comp_['z']]))
            [nx_rot, ny_rot, nz_rot] = self.rotate3D_point([nx, ny, nz],
                                                           theta, axis_of_rot)
            nx_rot += soma_loc[0]
            nz_rot += soma_loc[1]
            
            for c in self.morph.children_of(comp_):
                cx, cy, cz = self.shift_origin(
                    np.array([c['x'], c['y'], c['z']]))
                [cx_rot, cy_rot, cz_rot] = self.rotate3D_point([cx, cy, cz],
                                                               theta, axis_of_rot)
                # Make neurites get thinner with distance
                if morph_dist_arr:
                    dist_ = np.linalg.norm(np.array([cx, cy, cz]))
                    lw = lw_min+(lw_max-lw_min)*(max_dist-dist_)/max_dist
                else:
                    lw=lw_max
                
                
                # Shift the origin to the desired soma location
                cx_rot += soma_loc[0]
                cz_rot += soma_loc[1]

                
                linewidths.append(lw)
                all_lines_x,all_lines_y = (nx_rot, nz_rot),(cx_rot, cz_rot)
                all_lines.append([all_lines_x,all_lines_y])
                colors.append(color_dict[comp_['type']])

        lc = mc.LineCollection(all_lines, colors=colors, linewidths=linewidths,alpha=alpha)
        ax.add_collection(lc)            
        shifted_soma = self.shift_origin(self.soma_coord)
        
        
        soma_rad = kwargs.get('soma_rad') or 40
        shifted_soma[0] += soma_loc[0]
        shifted_soma[2] += soma_loc[1]
        ax.scatter(shifted_soma[0], shifted_soma[2],
               marker='o', alpha=alpha,edgecolor=color_dict[1],linewidth=lw_max,
               s=soma_rad,color=color_dict[1])
        
        if kwargs.get('axis_off'):
            ax.axis('off')
        return ax


    def draw_morphology(self, theta, axis_of_rot, reject_axon=True,soma_loc= np.array([0,0,0]),
                        **kwargs):
        color_dict = kwargs.get('color_dict')
        if not color_dict:
            color_dict = {4: 'orange', 3: 'darkred', 2: 'royalblue', 1: 'dimgrey'}
        
        morph_dist_arr = kwargs.get('morph_dist_arr')
        if morph_dist_arr:
            _,max_dist = np.min(morph_dist_arr),np.max(morph_dist_arr)
            lw_min = kwargs.get('lw_min') or .2
        lw_max = kwargs.get('lw') or 1
        alpha = kwargs.get('alpha') or 1
        ax = kwargs.get('ax')
        
        if not ax:
            sns.set(style='whitegrid')
            fig = plt.figure(figsize=(3, 8), dpi=100)        
            ax = fig.add_subplot(111, projection='3d')
        ax.axis('equal')

        all_x, all_y, all_z = [], [], []
        all_lines =[]
        linewidths,colors = [],[]
        
        for comp_ in self.morph.compartment_list:
            if reject_axon and comp_['type'] == 2:
                continue
            nx, ny, nz = self.shift_origin(np.array([comp_['x'], comp_['y'], comp_['z']]))
            [nx_rot, ny_rot, nz_rot] = self.rotate3D_point([nx, ny, nz],
                                                           theta, axis_of_rot)
            nx_rot += soma_loc[0]
            ny_rot += soma_loc[1]
            nz_rot += soma_loc[2]
            for c in self.morph.children_of(comp_):
                cx, cy, cz = self.shift_origin(
                    np.array([c['x'], c['y'], c['z']]))
                [cx_rot, cy_rot, cz_rot] = self.rotate3D_point([cx, cy, cz],
                                                               theta, axis_of_rot)
                # Make neurites thinner with distance
                if morph_dist_arr:
                    dist_ = np.linalg.norm(np.array([cx, cy, cz]))
                    lw = lw_min+(lw_max-lw_min)*(max_dist-dist_)/max_dist
                else:
                    lw=lw_max
                    
                
                # Shift the origin to the desired soma location
                cx_rot += soma_loc[0]
                cy_rot += soma_loc[1]
                cz_rot += soma_loc[2]
                
                # For animation: axis limits are automatically determined here
                ax.plot([nx_rot, cx_rot], [ny_rot, cy_rot],[nz_rot, cz_rot],
                        color=color_dict[comp_['type']], lw=lw, alpha=alpha)
                
                linewidths.append(lw)
                point1,point2 = (nx_rot,ny_rot,nz_rot),(cx_rot,cy_rot,cz_rot)
                all_lines.append([point1,point2])
                all_x.extend((nx_rot, cx_rot))
                all_y.extend((ny_rot, cy_rot))
                all_z.extend((nz_rot, cz_rot))
                colors.append(color_dict[comp_['type']])

        # For efficiently plotting a large no. of lines
#        lc = Line3DCollection(all_lines, colors=colors, linewidths=linewidths,alpha=1)
#        ax.add_collection(lc) 
        
        shifted_soma = self.shift_origin(self.soma_coord)
        shifted_soma_coord = (shifted_soma[0]+soma_loc[0], shifted_soma[1]+soma_loc[1],
                              shifted_soma[2]+soma_loc[2])

        soma_rad = kwargs.get('soma_rad') or 40
        if kwargs.pop('draw_sphere', None):
            (xs, ys, zs) = self.draw_sphere(shifted_soma_coord)
            ax.plot_surface(xs, ys, zs, rstride=1, cstride=1,color=color_dict[1],
                            linewidth=0,alpha=alpha)
            
        else:
            ax.scatter(*shifted_soma_coord,marker='o', alpha=alpha,edgecolor=color_dict[1],
                       linewidth=lw_max,s=soma_rad,color=color_dict[1])

        all_x_min, all_x_max = min(all_x), max(all_x)
        all_y_min, all_y_max = min(all_y), max(all_y)
        all_z_min, all_z_max = min(all_z), max(all_z)
        
        # Set axis limits for Line3DCollection
#        ax.set_xlim([all_x_min,all_x_max])
#        ax.set_ylim([all_y_min,all_y_max])
#        ax.set_zlim([all_z_min,all_z_max])
        
        z_mid = (all_z_min+all_z_max)*0.5
        view_base = np.max([all_x_min, all_y_min, all_x_max, all_y_max])
        elev_angle = np.arctan(z_mid/view_base)
        
        if kwargs.get('axis_off'):
            ax.axis('off')
        ax.view_init(elev=elev_angle, azim=90)

        return ax,elev_angle

