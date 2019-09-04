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

        for n in self.morph.compartment_list:
            x_coord, y_coord, z_coord = n['x'], n['y'], n['z']
            if n['type'] == 2:
                x_axon.append(x_coord)
                y_axon.append(y_coord)
                z_axon.append(z_coord)

            if reject_axon and n['type'] == 2:
                continue

            x.append(x_coord)
            y.append(y_coord)
            z.append(z_coord)

            if n['type'] == 4:
                x_apical.append(x_coord)
                y_apical.append(y_coord)
                z_apical.append(z_coord)

        morph_data = self.shift_origin(np.array(np.column_stack((x, y, z))))
        morph_apical = self.shift_origin(np.array(np.column_stack((x_apical,
                                                                   y_apical, z_apical))))
        morph_axon = self.shift_origin(np.array(np.column_stack((x_axon,
                                                                 y_axon, z_axon))))

        return morph_data, morph_apical, morph_axon

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
        u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        # shift and scale sphere
        x = rad*x + xCenter
        y = rad*y + yCenter
        z = rad*z + zCenter
        return (x, y, z)

    def draw_morphology(self, theta, axis_of_rot, reject_axon=True, **kwargs):
        color_dict = {4: 'orange', 3: 'darkred', 2: 'royalblue', 1: 'grey'}
        label_dict = {4: 'apical dendrite', 3: 'basal dendrite',
                      2: 'axon', 1: 'soma'}

        fig = plt.figure(figsize=(8, 8), dpi=100)
        ax = fig.add_subplot(111, projection='3d')

        all_x, all_y, all_z = [], [], []

        for n in self.morph.compartment_list:
            if reject_axon and n['type'] == 2:
                continue
            nx, ny, nz = self.shift_origin(np.array([n['x'], n['y'], n['z']]))
            [nx_rot, ny_rot, nz_rot] = self.rotate3D_point([nx, ny, nz],
                                                           theta, axis_of_rot)
            for c in self.morph.children_of(n):
                cx, cy, cz = self.shift_origin(
                    np.array([c['x'], c['y'], c['z']]))
                [cx_rot, cy_rot, cz_rot] = self.rotate3D_point([cx, cy, cz],
                                                               theta, axis_of_rot)

                ax.plot([nx_rot, cx_rot], [ny_rot, cy_rot],
                        [nz_rot, cz_rot],
                        color=color_dict[n['type']], lw=.8, alpha=.5,
                        label=label_dict[n['type']])

                all_x.extend((nx_rot, cx_rot))
                all_y.extend((ny_rot, cy_rot))
                all_z.extend((nz_rot, cz_rot))

        shifted_soma = self.shift_origin(self.soma_coord)
        shifted_soma_coord = (shifted_soma[0], shifted_soma[1],
                              shifted_soma[2])

        if kwargs.pop('draw_sphere', None):
            (xs, ys, zs) = self.draw_sphere(shifted_soma_coord)
            ax.plot_surface(xs, ys, zs, color=color_dict[1])
        else:
            ax.scatter(shifted_soma[0], shifted_soma[1],
                       shifted_soma[2],
                       marker='o', alpha=0.8,
                       #               s=math.pi*self.soma_rad,
                       s=100,
                       color=color_dict[1])

        all_x_min, all_x_max = min(all_x), max(all_x)
        all_y_min, all_y_max = min(all_y), max(all_y)
        all_z_min, all_z_max = min(all_z), max(all_z)
        z_mid = (all_z_min+all_z_max)*0.5
        view_base = np.max([all_x_min, all_y_min, all_x_max, all_y_max])
        elev_angle = np.arctan(z_mid/view_base)

        ax.axis('off')
        ax.view_init(elev=elev_angle, azim=90)
        if kwargs.get('figname'):
            figname = kwargs['figname']
        else:
            figname = 'Morph_figure.pdf'
        fig.savefig(figname, bbox_inches='tight')
        plt.close(fig)
