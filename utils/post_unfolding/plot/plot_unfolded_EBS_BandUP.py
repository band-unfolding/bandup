#! /usr/bin/env python
# Copyright (C) 2013-2015 Paulo V. C. Medeiros, Jonas Bjork
# This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
#
# BandUP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
#  BandUP is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import division
from __future__ import print_function
import argparse
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from os.path import join, abspath, dirname, splitext, basename
import sys
import json
from fractions import Fraction
from subprocess import Popen, PIPE

class BandUpPlotOptions(argparse.ArgumentParser):
    def __init__(self):
        super(BandUpPlotOptions, self).__init__()

        self.indent = 4 * ' '
        self.add_argument('input_file', default='unfolded_EBS_symmetry-averaged.dat', nargs='?', 
                          help='Name of the input file.')
        self.add_argument('output_file', nargs='?', default=None, 
                           help='Optional: Name of the output file. If not given, it will be based on ' 
                                'the name of the input file.')
        self.add_argument('-kpts', '--kpoints_file', default='KPOINTS_prim_cell.in', 
                          help='Name of the file containing information'
                               'about the primitive cell k-points. Default: KPOINTS_prim_cell.in')
        self.add_argument('-pc_file', '--prim_cell_file', default='prim_cell_lattice.in', 
                          help='Name of the file containing information'
                               'about the primitive cell lattice vectors. Default: prim_cell_lattice.in')
        self.add_argument('-efile', '--energy_info_file', default='energy_info.in', 
                          help='Name of the file containing information about the energy grid and '
                               'Fermi energy to be used. Default: energy_info.in. '
                                'This file is optional for the plotting tool.')
        # Colormaps
        self.cmap_names, self.cmaps = self.get_available_cmaps() 
        # Change self.default_cmap if you want another colormap to be the default one. 
        # I used 'gist_ncar' in my paper. 
        self.default_cmap = 'gist_ncar' 
        if('gist_ncar' not in self.cmap_names):
            self.default_cmap = self.cmap_names[0]
        self.possible_cmap_choices = self.add_mutually_exclusive_group()
        self.possible_cmap_choices.add_argument('-cmap', '--colormap', default=None, choices=self.cmap_names, 
                                                help='Choice of colormap for the plots (name of the colormap).'
                                                     'You might want to try a few.', metavar='colormap_name')
        self.possible_cmap_choices.add_argument('-icmap', '--icolormap', type=int, default=None, choices=range(len(self.cmap_names)), 
                                                help='Choice of colormap for the plots (integer number).', metavar='0~'+str(len(self.cmap_names) - 1))
        # E-Fermi and high-symm lines
        mpl_colors = sorted(mpl.colors.cnames.keys())
        self.add_argument('--high_symm_linecolor', default=None, choices=mpl_colors, 
                          help='Color of the lines marking the positions of the high-symmetry points.', 
                          metavar='some_matplotlib_color_name')
        self.add_argument('--e_fermi_linecolor', default=None, choices=mpl_colors, 
                          help='Color of the Fermi energy line.', metavar='some_matplotlib_color_name')
        self.add_argument('--line_width_E_f', type=float, default=None)
        self.add_argument('--line_width_high_symm_points', type=float, default=None)
        linestyles = ['solid', 'dashed', 'dashdot', 'dotted', '-', '--', '-.', ':']
        self.add_argument('--line_style_high_symm_points', default=None, choices=linestyles)
        self.add_argument('--line_style_E_f', default=None, choices=linestyles)
        self.add_argument('--tick_marks_size', type=float, default=8.0)

        self.add_argument('--save', action='store_true', default=False, 
                          help='Saves the figue to a file. Default: False')
        self.add_argument('--show', action='store_true', default=False, 
                          help='Shows the figue. Default: False if --save is selected, True otherwise.')
        self.add_argument('--saveshow', action='store_true', default=False, 
                          help='Saves the figue to a file and opens the saved file instead of a pyplot window. Default: False')
        self.add_argument('-res', '--fig_resolution', default='m', choices=['l','m','h'], 
                          help='Resolution of the figure: l = 100 dpi, m = 300 dpi, h = 600 dpi.'
                               'Default: m = 300 dpi')
        self.add_argument('--disable_auto_round_vmin_and_vmax', action='store_true', 
                          help='Disable normalization of the vmin and vmax of the color scale to ' 
                               'integer numbers.')

        # File format of the output figure
        self.default_fig_format = self.set_default_fig_format('tiff') # Change this line if you want another format to be the default.
        self.add_argument('-fmt', '--file_format', default=self.default_fig_format,
                          help='File format of the figure. Default: ' + self.default_fig_format)

        self.possible_fig_orientations = self.add_mutually_exclusive_group()
        self.possible_fig_orientations.add_argument('--landscape', action='store_true', default=False)
        self.possible_fig_orientations.add_argument('--portrait', action='store_true', default=False)

        self.add_argument('-ar', '--aspect_ratio', type=type(Fraction('3/4')), default=Fraction('3/4'), 
                          help='Aspect ratio of the generated plot. Default: 3/4')
        self.add_argument('--no_ef', action='store_true',help='Hides the E-Fermi line.')

        self.add_argument('--no_symm_lines', action='store_true',help='Hides the high-symmetry lines.')
        self.add_argument('--no_symm_labels', action='store_true',help='Hides the high-symmetry labels.')
        self.add_argument('--no_cb', action='store_true',help='Hides the colorbar.')


        self.plot_spin_projections = self.add_mutually_exclusive_group()
        self.plot_spin_projections.add_argument('--plot_sigma_x', action='store_true', default=False)
        self.plot_spin_projections.add_argument('--plot_sigma_y', action='store_true', default=False)
        self.plot_spin_projections.add_argument('--plot_sigma_z', action='store_true', default=False)
        self.plot_spin_projections.add_argument('--plot_spin_perp', action='store_true', default=False)
        self.plot_spin_projections.add_argument('--plot_spin_para', action='store_true', default=False)
        self.add_argument('--clip_spin', type=float, default=None)
        self.add_argument('--spin_marker', default='_', choices=['_', 'o'])


        self.add_argument('-nlev', '--n_levels', type=int, default=101, 
                          help='Number of different levels used in the contour plot. Default: 101')
        self.add_argument('-kmin', type=float, default=None)
        self.add_argument('-kmax', type=float, default=None)
        self.add_argument('-shift_k', '--shift_kpts_coords', type=float, default=0.0, 
                          help='Shift in the k-points. Default: 0.0')

        self.add_argument('-emin', type=float, default=None)
        self.add_argument('-emax', type=float, default=None)
        self.add_argument('-dE', type=float, default=None)
        self.add_argument('-shift_e', '--shift_energy', type=float, default=0.0,
                          help='Shift in the energy grid. Default: 0.0')

        self.add_argument('-interp', '--interpolation', default=None, choices=['nearest', 'linear', 'cubic'], 
                          help='Interpolation scheme used. Default: nearest')

        self.add_argument('-vmin', '--minval_for_colorbar', type=float, 
                          help='Value to which the first color of the colormap will be normalized.')
        self.add_argument('-vmax', '--maxval_for_colorbar', type=float, 
                          help='Value to which the last color of the colormap will be normalized.')
        self.add_argument('--round_cb', type=int, default=1, 
                          help='Number of decimal digits displayed in the colobar ticks.')

        self.show_cb_label_options = self.add_mutually_exclusive_group()
        self.show_cb_label_options.add_argument('--cb_label', action='store_true', 
                                                help='Show the colorbar label.')
        self.show_cb_label_options.add_argument('--cb_label_full', action='store_true', 
                                                help='Show the colorbar label (full).')

        
        self.add_argument('--skip_grid_freq_clip', action='store_true', help=argparse.SUPPRESS)
        # Flag to know if this script has been called by the GUI
        self.add_argument('--running_from_GUI', action='store_true', help=argparse.SUPPRESS)

        # Parsing command-line arguments and flags
        self.args = self.parse_args() 

        if(self.args.icolormap is not None):
            self.args.colormap = self.cmap_names[self.args.icolormap]
        self.using_default_cmap = self.args.colormap is None
        if(self.using_default_cmap):
            self.args.colormap = self.default_cmap

        if(self.args.saveshow):
            self.args.save = True

        if(self.args.output_file is not None):
            self.args.output_file = abspath(self.output_file_with_supported_extension(self.args.output_file))
            self.args.save = True

        if(not self.args.save and not self.args.show):
            self.args.save = False
            self.args.show = True

        self.args.aspect_ratio = float(self.args.aspect_ratio)
        if(self.args.aspect_ratio > 1.0):
            self.args.aspect_ratio = 1.0 / self.args.aspect_ratio


    def get_available_cmaps(self):
        colormap_names = sorted(plt.cm.datad.keys(), key=lambda s: s.lower())
        colormaps = dict([[cmap_name, plt.get_cmap(cmap_name)] for cmap_name in colormap_names])
        # Custom colormaps - I'm still testing
        try:
            custom_cm_folder = os.path.join(os.environ['BANDUPPLOTPATH'], 'custom_colormaps')
            os_listdir_full_custom_cm_folder = [os.path.abspath(os.path.join(custom_cm_folder, cmap_file)) for cmap_file in os.listdir(custom_cm_folder)]
            custom_cmap_files = [cmap_file for cmap_file in os_listdir_full_custom_cm_folder if cmap_file.endswith('.cmap')]
            custom_cmaps = []
            for cmap_file in custom_cmap_files:
                cmap_name = os.path.splitext(os.path.basename(cmap_file))[0]
                color_dict = json.load(open(cmap_file))
                colormap_names.append(cmap_name)
                colormaps[cmap_name] = plt.cm.colors.LinearSegmentedColormap(cmap_name, color_dict, 2048)
        except:
            pass

        return sorted(colormap_names), colormaps


    def allowed_fig_formats(self):
        temp_fig = plt.figure()
        all_allowed_filetypes = [f.lower() for f in temp_fig.canvas.get_supported_filetypes().keys()]
        plt.close(temp_fig)
        preferred_filetypes = ['tiff', 'png', 'bmp', 'jpg', 'pdf', 'eps']
        allowed_filetypes = [fmt for fmt in preferred_filetypes if fmt in all_allowed_filetypes]
        allowed_filetypes += [fmt for fmt in all_allowed_filetypes if fmt not in preferred_filetypes]
        return allowed_filetypes


    def set_default_fig_format(self, default_fig_format='tiff'):
        allowed_filetypes = self.allowed_fig_formats()
        if default_fig_format not in allowed_filetypes:
            default_fig_format = allowed_filetypes[0]
        return default_fig_format


    def output_file_with_supported_extension(self, output_file_name):
        # Defines the file format of the generated figure according to what your system supports.
        allowed_filetypes = self.allowed_fig_formats()
        file_name = splitext(basename(output_file_name))[0]
        file_extension = splitext(output_file_name)[1][1:]
        if file_extension.lower() not in allowed_filetypes:
            print ('WARNING: The file type you requested for the output figure (%s) is not '
                   "supported by your system and has been changed to '%s'." 
                   % (file_extension, self.default_fig_format))
            print ('         * The choices supported by your system are: %s' 
                   % ", ".join(map(str,allowed_filetypes)))
            output_file_name = file_name + '.' + default_fig_format
        return output_file_name



class BandUpPlot():
    def __init__(self, plot_options):
        self.plot_options = plot_options
        self.indent = plot_options.indent
        self.args = plot_options.args
        args = self.args
        self.KptsCoords, self.energies, self.delta_Ns, self.spin_projections = self.read_BandUP_output()
        self.dE_for_hist2d = self.get_dE_for_grid()
        self.kmin, self.kmax, self.emin, self.emax = self.define_plot_boundaries()
        self.KptsCoords, self.energies, self.delta_Ns, self.spin_projections = self.reduced_read_data()
        self.cmap_name = args.colormap
        self.cmap = plot_options.cmaps[self.cmap_name]

        self.pos_high_symm_points = None
        if(not args.no_symm_lines or not args.no_symm_labels):
            self.pos_high_symm_points, self.labels_high_symm_lines = self.get_pos_and_labels_high_symm_points()
            if(not args.no_symm_lines):
                print ('Vertical lines will be automatically drawn at: k = %s' \
                      % ', '.join(map("{:9.5f}".format, self.pos_high_symm_points)))
        
        self.title = ''
        self.x_axis_label = ''
        self.y_axis_label = '$\epsilon \hspace{0.25} - \hspace{0.25} \epsilon _{_F} (eV)$'

        self.scaling_factor_fig_size = 2.0 # 1 gives a width of 8.6cm for the final figure.
        self.title_size = 10*self.scaling_factor_fig_size
        self.yaxis_labels_size = 14*self.scaling_factor_fig_size
        self.xaxis_labels_size = 10*self.scaling_factor_fig_size
        self.tick_marks_size = args.tick_marks_size * self.scaling_factor_fig_size
        self.colorbar_label_size = 8*self.scaling_factor_fig_size
        self.colorbar_tick_marks_size = 7*self.scaling_factor_fig_size

        if(args.line_style_high_symm_points is None):
            self.line_style_high_symm_points = 'solid' # Change here if you want another default
        else:
            self.line_style_high_symm_points = args.line_style_high_symm_points
        if(args.line_width_high_symm_points is None):
            self.line_width_high_symm_points = 0.6*self.scaling_factor_fig_size
        else:
            self.line_width_high_symm_points = args.line_width_high_symm_points

        if(args.line_style_E_f is None):
            self.line_style_E_f = 'dotted' # Change here if you want another default
        else:
            self.line_style_E_f = args.line_style_E_f
        if(args.line_width_E_f is None):
            self.line_width_E_f = 0.8*self.scaling_factor_fig_size
        else:
            self.line_width_E_f = args.line_width_E_f

        # Position of the line y = E_F in the figure
        self.E_f = 0.000

        if(args.landscape):
            self.cb_orientation = 'vertical'
            self.offset_x_text_colorbar = 3.5
            self.offset_y_text_colorbar = 0.5
        else:
            self.cb_orientation = 'horizontal'
            self.offset_x_text_colorbar = 0.5
            self.offset_y_text_colorbar = -2.5
            self.xaxis_labels_size = self.xaxis_labels_size/args.aspect_ratio
            self.yaxis_labels_size = self.yaxis_labels_size/args.aspect_ratio
            self.tick_marks_size = self.tick_marks_size/args.aspect_ratio
            self.colorbar_label_size = self.colorbar_label_size/args.aspect_ratio
            self.colorbar_tick_marks_size = self.colorbar_tick_marks_size/args.aspect_ratio
            self.line_width_high_symm_points = self.line_width_high_symm_points/args.aspect_ratio
            self.line_width_E_f = self.line_width_E_f/args.aspect_ratio


        self.std_column_widths_cm = [8.6, 17.8]
        self.std_column_widths_inch = [item/2.54 for item in self.std_column_widths_cm]
        self.fig_widths_options_in_inches = [width * self.scaling_factor_fig_size for width in self.std_column_widths_inch] 
        self.fig_width_inches = self.fig_widths_options_in_inches[0]
        if(args.landscape):
            self.fig_height_inches = self.fig_width_inches * args.aspect_ratio
        else:
            self.fig_height_inches = self.fig_width_inches / args.aspect_ratio
            if(not args.portrait and not args.running_from_GUI):
                print ('Assuming portrait orientation for the output figure.')

    def __color_E_f_and_high_symm_lines(self, image, args_linecolor):
        # Defining the colors of the vertical and horizontal lines.
        # Don't take this part too seriously. Change it if you need/want.
        if(args_linecolor is not None):
            color = args_linecolor
        else:
            dist_black = np.linalg.norm(np.array(self.cmap(image.norm.vmin)[:3]))
            dist_white = np.linalg.norm(np.array(self.cmap(image.norm.vmin)[:3]) - np.array([1,1,1]))
            if(dist_black < dist_white):
                color = 'white'
            else:
                color = 'black'
        return color
    def color_E_f_line(self, image):
        return self.__color_E_f_and_high_symm_lines(image, self.args.e_fermi_linecolor)
    def color_high_symm_lines(self, image):
        return self.__color_E_f_and_high_symm_lines(image, self.args.high_symm_linecolor)
        

    def read_BandUP_output(self):
        print ('Reading input file "%s"' % self.args.input_file)
        
        spin_projections = None
        spin_projections_read = False
        failed_reading_spin_projections = False
        try:
            if(self.args.plot_sigma_x):
                KptsCoords, energies, delta_Ns, spin_projections = np.loadtxt(self.args.input_file, usecols=(0,1,2,3),unpack=True)
                spin_projections_read = True
            elif(self.args.plot_sigma_y):
                KptsCoords, energies, delta_Ns, spin_projections = np.loadtxt(self.args.input_file, usecols=(0,1,2,4),unpack=True)
                spin_projections_read = True
            elif(self.args.plot_sigma_z):
                KptsCoords, energies, delta_Ns, spin_projections = np.loadtxt(self.args.input_file, usecols=(0,1,2,5),unpack=True)
                spin_projections_read = True
            elif(self.args.plot_spin_perp):
                KptsCoords, energies, delta_Ns, spin_projections = np.loadtxt(self.args.input_file, usecols=(0,1,2,6),unpack=True)
                spin_projections_read = True
            elif(self.args.plot_spin_para):
                KptsCoords, energies, delta_Ns, spin_projections = np.loadtxt(self.args.input_file, usecols=(0,1,2,7),unpack=True)
                spin_projections_read = True
        except IndexError:
            failed_reading_spin_projections = True
            pass
        if(not spin_projections_read):
            KptsCoords, energies, delta_Ns = np.loadtxt(self.args.input_file, usecols=(0,1,2),unpack=True)
            if(failed_reading_spin_projections):
                print (self.indent + 'WARNING: Could not read spin info.')


        print (self.indent + '* Max. delta_N:          ' + "{:7.3f}".format(np.max(delta_Ns)))
        print (self.indent + '* Min. non-zero delta_N: ' + "{:7.3f}".format(min([value for value in delta_Ns if abs(value) > 0.95E-3])))
        if(spin_projections is not None):
            print (self.indent + 'Plotting also spin info (projections of the expectation values of Pauli matrices).')
            print (2 * self.indent + '* Maxval of the projections: ' + \
                   "{:7.3f}".format(np.max(spin_projections)))
            try:
                min_proj = min([value for value in spin_projections if abs(value) > 0.95E-3])
            except ValueError:
                min_proj = 0.0 
            print (2 * self.indent + '* Minval of the projections:' + \
                   "{:7.3f}".format(min_proj))
            if(abs(min_proj) > 0.95E-3 or abs(np.max(spin_projections)) > 0.95E-3):
                print (2 * self.indent + '* Positive projection values will be shown in blue.')
                print (2 * self.indent + '* Negative projection values will be shown in red.')

        print ('Eliminating duplicated points...')
        sorted_index_data = np.lexsort((energies, KptsCoords))
        KptsCoords = KptsCoords[sorted_index_data] + self.args.shift_kpts_coords
        energies = energies[sorted_index_data] + self.args.shift_energy
        delta_Ns = delta_Ns[sorted_index_data]
        if(spin_projections is not None):
            spin_projections = spin_projections[sorted_index_data]

        new_KptsCoords = []
        new_energies = []
        new_delta_Ns = []
        new_spin_projections = []
        for line_number in range(len(KptsCoords)-1):
            skip_line = False
            if(KptsCoords[line_number]==KptsCoords[line_number+1] and 
               energies[line_number]==energies[line_number+1]):
               delta_Ns[line_number+1] = max(delta_Ns[line_number],delta_Ns[line_number+1])
               try:
                   spin_projections[line_number+1] = max(spin_projections[line_number],spin_projections[line_number+1])
               except TypeError:
                   pass
               skip_line = True
            if(not skip_line):
                new_energies.append(energies[line_number])
                new_KptsCoords.append(KptsCoords[line_number])
                new_delta_Ns.append(delta_Ns[line_number])
                try:
                    new_spin_projections.append(spin_projections[line_number])
                except TypeError:
                   pass
        # Appending the last point
        new_energies.append(energies[-1])
        new_KptsCoords.append(KptsCoords[-1])
        new_delta_Ns.append(delta_Ns[-1])
        try:
            new_spin_projections.append(spin_projections[-1])
        except TypeError:
           pass

        KptsCoords = np.array(new_KptsCoords)
        energies = np.array(new_energies)
        delta_Ns = np.array(new_delta_Ns)
        if(new_spin_projections): 
            spin_projections = np.array(new_spin_projections)
        else:
            spin_projections = None
        print (self.indent + '* Done with eliminating duplicated points.')

        return KptsCoords, energies, delta_Ns, spin_projections

    def reduced_read_data(self):
        KptsCoords, energies, delta_Ns, spin_projections = self.KptsCoords, self.energies, self.delta_Ns, self.spin_projections
        size_of_old_data = float(len(self.KptsCoords))
        if(self.kmin>self.KptsCoords.min() or self.kmax<self.KptsCoords.max() or 
           self.emin>self.energies.min() or self.emax<self.energies.max()):
            print ('Trying to reduce data such that only the needed part of it is parsed.')
            new_KptsCoords = []
            new_energies = []
            new_delta_Ns = []
            new_spin_projections = []
            for iline in range(len(self.KptsCoords)):
                if (self.KptsCoords[iline] >= self.kmin and self.KptsCoords[iline] <= self.kmax and 
                    self.energies[iline] >= self.emin and self.energies[iline] <= self.emax):
                    new_KptsCoords.append(self.KptsCoords[iline])
                    new_energies.append(self.energies[iline])
                    new_delta_Ns.append(self.delta_Ns[iline])
                    if(spin_projections is not None):
                        new_spin_projections.append(self.spin_projections[iline])
            KptsCoords = np.array(new_KptsCoords)
            energies = np.array(new_energies)
            delta_Ns = np.array(new_delta_Ns)
            if(new_spin_projections):
                spin_projections = np.array(new_spin_projections)
            size_of_new_data = float(len(KptsCoords))
            print (self.indent + '* Done. Working with %.2f' % (100.0*size_of_new_data/size_of_old_data), \
                   '% of the data read in.')
        return KptsCoords, energies, delta_Ns, spin_projections


    def define_plot_boundaries(self):

        args = self.args
        kmin = args.kmin
        kmax = args.kmax
        emin = args.emin
        emax = args.emax
        if(any(value is None for value in [kmin, kmax, emin, emax])):
            # If the energy_info_file file exists, then it overrides the 
            # defaults based on the EBS file
            try:
                with open(args.energy_info_file) as energy_info_file:
                    energy_info_file_lines = energy_info_file.readlines()
            except IOError:
                pass

        if(kmin is None):
            try:
                kmin = float(energy_info_file_lines[1].split()[1])
            except:
                kmin = self.KptsCoords.min()
        if(kmax is None):
            try:
                kmax = float(energy_info_file_lines[2].split()[1])
            except:
                kmax = self.KptsCoords.max()

        if(emin is None):
            try:
                emin = float(energy_info_file_lines[1].split()[0])
            except:
                emin = self.energies.min()
        if(emax is None):
            try:
                emax = float(energy_info_file_lines[2].split()[0])
            except:
                emax = self.energies.max()

        if kmin > kmax:
            kmin, kmax = kmax, kmin
        if(kmin < self.KptsCoords.min()): 
            kmin = self.KptsCoords.min()
            print ('WARNING: Resetting k_min to %s.' % "{:9.5f}".format(kmin))
        if(kmax > self.KptsCoords.max()): 
            kmax = self.KptsCoords.max()
            print ('WARNING: Resetting k_max to %s.' % "{:9.5f}".format(kmax))
        if emin > emax:
            emin, emax = emax, emin
        # Setting emin and emax to the nearest values in the energy grid
        energies = sorted(set(self.energies))
        dE = energies[1] - energies[0]
        n_emin = int(round((emin - energies[0]) / dE))
        n_emax = int(round((emax - energies[0]) / dE))
        emin = energies[0] + n_emin * dE
        emax = energies[0] + n_emax * dE

        return kmin, kmax, emin, emax


    def get_dE_for_grid(self):
        if(self.plot_options.args.dE is not None):
            # Priority #1 for the command line arguments
            dE_for_hist2d = self.plot_options.args.dE
        else:
            try:
                # Priority #2 for the energy_info_file 
                with open(self.plot_options.args.energy_info_file) as energy_info_file:
                    energy_info_file_lines = energy_info_file.readlines()
                    dE_for_hist2d = abs(float(energy_info_file_lines[3].split()[0]))
                    if dE_for_hist2d == 0:
                        raise ValueError
            except(IOError, ValueError):
            # Priority #3 for the energy grid found on the EBS file
                energies = sorted(set(self.energies))
                dE_for_hist2d = energies[1] - energies[0]

        return dE_for_hist2d


    def get_pos_and_labels_high_symm_points(self):
        args = self.plot_options.args
        # Determining the positions of high-symmetry BZ points on the plot
        kpts_file_lines = []
        with open(args.kpoints_file,'r') as kpts_file:
            for line in kpts_file:
                if line.strip(): # Skipping blank lines
                    kpts_file_lines.append(line)
        # Checking if the parameter a0 has been passed in the first line of the k-points file (old format)
        try:
            a0_informed_in_old_format = True
            latt_param_kpts_file_old_format = np.float(kpts_file_lines[0].split()[0])
        except:
            a0_informed_in_old_format = False

        try:
            zero_of_kpts_line = np.float(kpts_file_lines[2].split()[1])
        except:
            zero_of_kpts_line = 0.0

        coords_type = kpts_file_lines[3].strip()
        # Checking if the parameter a0 has been passed after the "cartesian" flag (new format)
        try:
            a0_informed_in_new_format = True
            latt_param_kpts_file_new_format = np.float(kpts_file_lines[3].split()[1]) 
        except:
            a0_informed_in_new_format = False

        a0_informed_in_kpts_file = a0_informed_in_new_format or a0_informed_in_old_format
        # If both old and new formats to inform a0 in the k-points file are used, then the new format holds
        if(a0_informed_in_old_format):
            latt_param_kpts_file = latt_param_kpts_file_old_format
        if(a0_informed_in_new_format):
            latt_param_kpts_file = latt_param_kpts_file_new_format

        read_k_start = []
        label_k_start = []
        for line in kpts_file_lines[4::2]:
            read_k_start.append(np.array([np.float(component) for component in line.split()[:3]]))
            if (line.find('!') >= 0):
                line_after_exclm_mark = line[line.find('!')+1:].strip()
                label_k_start.append(line_after_exclm_mark) 
            else:
                label_k_start.append('')
        read_k_start = np.array(read_k_start)
        ndirections = len(read_k_start)

        read_k_end = []
        label_k_end = []
        for line in kpts_file_lines[5::2]:
            read_k_end.append(np.array([np.float(component) for component in line.split()[:3]]))
            if (line.find('!') >= 0):
                line_after_exclm_mark = line[line.find('!')+1:].strip()
                label_k_end.append(line_after_exclm_mark)
            else:
                label_k_end.append('')
        read_k_end = np.array(read_k_end)

        with open(args.prim_cell_file) as prim_cell_lattice_file:
            prim_cell_lattice_file_lines = prim_cell_lattice_file.readlines()
        latt_param = np.float(prim_cell_lattice_file_lines[1])
        a1 = [latt_param*float(comp) for comp in prim_cell_lattice_file_lines[2].split()[:3]]
        a2 = [latt_param*float(comp) for comp in prim_cell_lattice_file_lines[3].split()[:3]]
        a3 = [latt_param*float(comp) for comp in prim_cell_lattice_file_lines[4].split()[:3]]
        cell_vectors = [a1,a2,a3]
        cell_volume = np.fabs(np.dot(a1,np.cross(a2,a3)))
        b1 = (2.0*np.pi)*np.cross(a2,a3)/cell_volume
        b2 = (2.0*np.pi)*np.cross(a3,a1)/cell_volume
        b3 = (2.0*np.pi)*np.cross(a1,a2)/cell_volume
        b_matrix_pc = [b1, b2, b3]

        k_start = [np.zeros(3) for dir in range(ndirections)]
        k_end = [np.zeros(3) for dir in range(ndirections)]
        for idir in range(ndirections):
            if(coords_type[0].upper() == 'C'):
                if(not a0_informed_in_kpts_file):
                    print ('ERROR: You have selected cartesian coordinates in your input k-points file, ' \
                           'but you have not passed a scaling parameter "a0".')
                    print ('       The actuall coordiates of the k-points are given by: '\
                           'ki[actual] = two_pi*ki[passed in file]/a0.')
                    print ('       Please write the value of a0 after your tag "' + \
                           coords_type + '", and run the code again.')
                    print ('Stopping now.')
                    sys.exit(0)
                k_start[idir] = (2.0*np.pi)*read_k_start[idir]/latt_param_kpts_file
                k_end[idir] = (2.0*np.pi)*read_k_end[idir]/latt_param_kpts_file
            else:
                if((coords_type[0].upper() != 'R') and (idir==0)):
                    print ('WARNING: Assuming that the pc-kpts have been informed in '\
                           'fractional (reciprocal) coordinates.')
                for i in range(3):
                    k_start[idir] += read_k_start[idir,i]*b_matrix_pc[i]
                    k_end[idir] += read_k_end[idir,i]*b_matrix_pc[i]


        pos_high_symm_points = [zero_of_kpts_line]
        for idir in range(0,ndirections):
            pos_high_symm_points.append(pos_high_symm_points[-1] + 
                                       np.linalg.norm(k_end[idir] - k_start[idir]))


        labels_high_symm_lines = [label_k_start[0]]
        for idir in range(1,ndirections):
            if(label_k_start[idir] == label_k_end[idir-1]):
                labels_high_symm_lines.append(label_k_start[idir])
            else:
                labels_high_symm_lines.append(label_k_end[idir-1]+','+label_k_start[idir])
        labels_high_symm_lines += [label_k_end[-1]]
        labels_high_symm_lines = self.greek_letters(labels_high_symm_lines)


        return pos_high_symm_points, labels_high_symm_lines


    def greek_letters(self, input_list):
        symb_dict = {'G':'$\Gamma$', 'GAMMA':'$\Gamma$',
                     'G-':'$\overline{\Gamma}$', 'GAMMA-':'$\overline{\Gamma}$',
                     'DELTA':'$\Delta$', 'DELTA-':'$\overline{\Delta}$',
                     'LAMBDA':'$\Lambda$', 'LAMBDA-':'$\overline{\Lambda}$',
                     'SIGMA':'$\Sigma$', 'SIGMA-':'$\overline{\Sigma}$'}
        output_list = input_list[:]
        for i in range(len(input_list)):
            if symb_dict.has_key(input_list[i].upper()):
                output_list[i] = symb_dict[input_list[i].upper()]
            if output_list[i].endswith('-'):
                output_list[i] = '$\overline{\mathrm{' + output_list[i][:-1] + '}}$'
        return output_list


def round(f,n=0):
    return '%.*f' % (n, f)

def print_opening_message():
    print ('                                                                                      \n'
           '===================================================================================== \n'
           '             BandUP: Band Unfolding code for Plane-wave based calculations            \n'
           '                  Copyright (C) 2013-2015 Paulo V. C. Medeiros                        \n'
           '                                                                                      \n'
           '              Post-processing utility "plot_unfolded_EBS_BandUP.py"                   \n'
           '            >>> Visualizing the unfolded EBSs produced by BandUP <<<                  \n'
           '===================================================================================== \n'
           'Copyright (C) 2013-2015 Paulo V. C. Medeiros, Jonas Bjork                             \n'
           '                        paume@ifm.liu.se, jonas.bjork@liu.se                          \n'
           '                        Computational Physics Division                                \n'
           '                        Department of Physics, Chemistry and Biology - IFM            \n'
           '                        Linkoping University                                          \n'
           '                        Sweden                                                        \n'
           'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                           \n'
           '===================================================================================== \n'
           '                                                                                      \n')

def make_plot(plot):
    indent = plot.plot_options.indent
    args = plot.plot_options.args
    # Creating the plot
    print ('Generating the plot...')
    fig = plt.figure(figsize=(plot.fig_width_inches,plot.fig_height_inches))
    ax = fig.add_subplot(111)
    # Defining the color schemes.
    print (indent + '>>> Using the "' + plot.cmap_name + '" colormap.')
    if(plot.plot_options.using_default_cmap and not args.running_from_GUI):
        print (2 * indent + 'Tip: You can try different colormaps by either:')
        print (2 * indent + '     * Running the plot tool with the option -icmap n, ' \
               'with n in the range from 0 to', len(plot.plot_options.cmaps) - 1)
        print (2 * indent + '     * Running the plot tool with the option "-cmap cmap_name".')
        print (2 * indent + '> Take a look at')
        print (4 * indent + '<http://matplotlib.org/examples/color/colormaps_reference.html>')
        print (2 * indent + '  for a list of colormaps, or run')
        print (4 * indent + '"./plot_unfolded_EBS_BandUP.py --help".')

    # Building the countour plot from the read data
    # Defining the (ki,Ej) grid.
    if(args.interpolation is not None):
        ki = np.linspace(plot.kmin, plot.kmax, 2 * len(set(plot.KptsCoords)) + 1, endpoint=True)
        Ei = np.arange(plot.emin, plot.emax + plot.dE_for_hist2d, plot.dE_for_hist2d)
        # Interpolating
        grid_freq = griddata((plot.KptsCoords, plot.energies), plot.delta_Ns, (ki[None,:], Ei[:,None]), 
                             method=args.interpolation, fill_value=0.0)
    else:
        ki = np.unique(np.clip(plot.KptsCoords, plot.kmin, plot.kmax))
        Ei = np.unique(np.clip(plot.energies, plot.emin,  plot.emax))
        grid_freq = griddata((plot.KptsCoords, plot.energies), plot.delta_Ns, (ki[None,:], Ei[:,None]), 
                             method='nearest', fill_value=0.0)

    if(not args.skip_grid_freq_clip):
        grid_freq = grid_freq.clip(0.0) # Values smaller than zero are just noise.
    # Normalizing and building the countour plot
    manually_normalize_colorbar_min_and_maxval = False
    if((args.maxval_for_colorbar is not None) or (args.minval_for_colorbar is not None)):
        manually_normalize_colorbar_min_and_maxval = True
        args.disable_auto_round_vmin_and_vmax = True
        maxval_for_colorbar = args.maxval_for_colorbar
        minval_for_colorbar = args.minval_for_colorbar
    else:
        if not args.disable_auto_round_vmin_and_vmax:
            minval_for_colorbar = float(round(np.min(grid_freq)))
            maxval_for_colorbar = float(round(np.max(grid_freq)))
            args.round_cb = 0
    if(manually_normalize_colorbar_min_and_maxval or not args.disable_auto_round_vmin_and_vmax):
        modified_vmin_or_vmax = False
        if not args.disable_auto_round_vmin_and_vmax and not args.running_from_GUI:
            print (plot.indent + '* Automatically renormalizing color scale '\
                   '(you can disable this with the option --disable_auto_round_vmin_and_vmax):')
        if manually_normalize_colorbar_min_and_maxval:
            print (plot.indent + '* Manually renormalizing color scale')
        if(minval_for_colorbar is not None):
            previous_vmin = np.min(grid_freq)
            if(abs(previous_vmin - minval_for_colorbar) >= 0.1):
                modified_vmin_or_vmax = True
                print (2 * indent + 'Previous vmin = %.1f, new vmin = %.1f' % (previous_vmin, 
                                                                               minval_for_colorbar))
        else:
            minval_for_colorbar = np.min(grid_freq)
        if(maxval_for_colorbar is not None):
            previous_vmax = np.max(grid_freq)
            if(abs(previous_vmax - maxval_for_colorbar) >= 0.1):
                modified_vmin_or_vmax = True
                print (2 * indent + 'Previous vmax = %.1f, new vmax = %.1f' % (previous_vmax, 
                                                                               maxval_for_colorbar))
        else:
            maxval_for_colorbar = np.max(grid_freq)
        if(modified_vmin_or_vmax):
            print (2 * indent + 'The previous vmin and vmax might be slightly different from '
                                'the min and max delta_Ns '
                                'due to the interpolation scheme used for the plot.')
        # values > vmax will be set to vmax, and #<vmin will be set to vmin 
        grid_freq = grid_freq.clip(minval_for_colorbar, maxval_for_colorbar)
        v = np.linspace(minval_for_colorbar, maxval_for_colorbar, args.n_levels, endpoint=True)
    else: 
        v = np.linspace(np.min(grid_freq), np.max(grid_freq), args.n_levels, endpoint=True)
    print (indent + '* Drawing contour plot...')
    print (2 * indent + '> Using %i color levels. Use the option "--n_levels" to choose a different number.' %args.n_levels)
    image = ax.contourf(ki, Ei, grid_freq, levels=v, cmap=plot.cmap)

    plot_spin_proj_requested = args.plot_spin_perp or args.plot_spin_para or args.plot_sigma_x or args.plot_sigma_y or args.plot_sigma_z
    if(plot_spin_proj_requested and plot.spin_projections is not None):
        print (indent + '* Drawing spin projection info')
        cmap_for_spin_plot = [plt.cm.bwr, plt.cm.RdBu, plt.cm.seismic_r][0]

        if(args.clip_spin is None):
            vmin_spin = np.min(plot.spin_projections)
            vmax_spin = np.max(plot.spin_projections)
        else:
            vmax_spin = abs(args.clip_spin)
            vmin_spin = -1.0 * abs(args.clip_spin)
            print (2 * indent + '* New maxval for spin: %.2f' % vmax_spin)
            print (2 * indent + '* New minval for spin: %.2f' % vmin_spin)
 
        spin_projections = np.clip(plot.spin_projections, vmin_spin, vmax_spin)
        grid_freq_spin = griddata((plot.KptsCoords, plot.energies), spin_projections, (ki[None,:], Ei[:,None]), 
                                  method='nearest', fill_value=0.0)

        k_for_scatter = [] 
        E_for_scatter = [] 
        spin_projections_for_scatter = [] 
        for iener in range(len(Ei)):
            for ikpt in range(len(ki)):
                if(abs(grid_freq_spin[iener, ikpt]) > 1E-3):
                    k_for_scatter.append(ki[ikpt])
                    E_for_scatter.append(Ei[iener])
                    spin_projections_for_scatter.append(grid_freq_spin[iener, ikpt])

        if(spin_projections_for_scatter):
            if(args.spin_marker=='o'):
                image2 = ax.scatter(k_for_scatter, E_for_scatter, marker='o', 
                                    s=[10.0 * abs(item) for item in spin_projections_for_scatter], 
                                    c=spin_projections_for_scatter, cmap=cmap_for_spin_plot)
            else:
                image2 = ax.scatter(k_for_scatter, E_for_scatter, marker='_', 
                                    s=[500.0 * (ki[1] - ki[0]) for item in spin_projections_for_scatter], 
                                    linewidth=[100.0 * plot.dE_for_hist2d * (item ** 2) for item in spin_projections_for_scatter], 
                                    c=spin_projections_for_scatter, cmap=cmap_for_spin_plot)
        else:
            print (2 * indent + '* The abs values of the spin projections were all < 1E-3.')

    #Preparing the plot
    ax.set_xlim(plot.kmin, plot.kmax)
    ax.set_ylim(plot.emin, plot.emax)
    ax.set_title(plot.title, fontsize=plot.title_size)
    ax.set_ylabel(plot.y_axis_label, fontsize=plot.yaxis_labels_size)
    plt.yticks(fontsize=plot.tick_marks_size)

    # Fermi energy line
    show_E_f = not args.no_ef
    if(show_E_f and plot.E_f >= plot.emin and plot.E_f <= plot.emax):
        E_f_line = plt.axhline(y=plot.E_f, c=plot.color_E_f_line(image), linestyle=plot.line_style_E_f, lw=plot.line_width_E_f)
    # High symmetry points lines
    if(plot.pos_high_symm_points):
        x_tiks_positions = [kx for kx in plot.pos_high_symm_points if kx - plot.kmax <= 1E-2 and kx >= plot.kmin]
    if(args.no_symm_labels):
        x_tiks_labels = []
    else:
        x_tiks_labels = [plot.labels_high_symm_lines[i] for i in range(len(plot.labels_high_symm_lines)) if 
                         plot.pos_high_symm_points[i] in x_tiks_positions]
        x_tiks_labels = [xlabel for xlabel in x_tiks_labels if xlabel]
    if x_tiks_labels:
        print (indent + '* K-point labels read from the "' + args.kpoints_file + '" file:')
        for ilabel in range(len(x_tiks_labels)):
            print(2 * indent + "k = {:9.5f}".format(x_tiks_positions[ilabel]) + ', label =',\
                   x_tiks_labels[ilabel])
        plt.xticks(x_tiks_positions, x_tiks_labels, fontsize=plot.tick_marks_size)
    else:
        plot.x_axis_label = '$k \hspace{0.25} (\AA^{-1})$'
        plt.locator_params(axis = 'x', nbins = 5)
        ax.set_xlabel(plot.x_axis_label, fontsize=plot.xaxis_labels_size)
        plt.xticks(fontsize=plot.tick_marks_size)
    ax.tick_params(axis='x', pad=10)

    # Drawing vertical lines at the positions of the high-symmetry points
    if(not args.no_symm_lines):
        for line_position in [pos for pos in plot.pos_high_symm_points if float(round(pos, 3)) > float(round(plot.kmin, 3)) and 
                                                                          float(round(pos, 3)) < float(round(plot.kmax, 3))]:
            hs_lines = plt.axvline(x=line_position, c=plot.color_high_symm_lines(image), linestyle=plot.line_style_high_symm_points, 
                                   lw=plot.line_width_high_symm_points)

    # Color bar
    show_colorbar = not args.no_cb
    if show_colorbar:
        if plot.cb_orientation=='vertical':
            cb_pad=0.005
        else:
            cb_pad=0.06
        if(not x_tiks_labels):
            cb_pad += 0.08 # To prevent the cb from overlapping with the numbers.

        cb_yticks = np.arange(int(image.norm.vmin), int(image.norm.vmax) + 1, 1)

        cb_ytick_labels = [round(item,abs(args.round_cb)) for item in cb_yticks]
        cb = plt.colorbar(image, ax=ax, ticks=cb_yticks, orientation=plot.cb_orientation, pad=cb_pad)
        cb.set_ticklabels(cb_ytick_labels)
        cb.ax.tick_params(labelsize=plot.colorbar_tick_marks_size)
       
        color_bar_label = None
        if args.cb_label: 
            color_bar_label = ('$Color scale: \hspace{0.5} \delta N(\\vec{k}; ' + 
                               '\hspace{0.25} \epsilon)$ ') 
        if args.cb_label_full: 
            color_bar_label = ('$Colors cale: \hspace{0.5} \delta N(\\vec{k}; ' +
                               '\hspace{0.25} \epsilon);$ '+ 
                               '$\delta\epsilon=' + round(1000.0*plot.dE_for_hist2d,0) +
                               '\\hspace{0.25} meV.$')

        if plot.cb_orientation=='vertical':
            cb_label_rotation = 90
        else:
            cb_label_rotation = 0
        if color_bar_label:
            cb.ax.text(plot.offset_x_text_colorbar, plot.offset_y_text_colorbar,
                       color_bar_label, rotation=cb_label_rotation, ha='center',
                       va='center', fontsize=plot.colorbar_label_size)

    # Saving/showing the results
    plt.tick_params(which='both', bottom='off', top='off', left='off', right='off', 
                    labelbottom='on')


    default_out_basename = "_".join([splitext(basename(args.input_file))[0], 'E_from', str(plot.emin), 'to', 
                                    str(plot.emax), 'eV_dE', 
                                    str(plot.dE_for_hist2d), 'eV'])
    if(args.save):
        if(args.output_file is None):
            args.output_file = abspath(default_out_basename + '.' + args.file_format)

        print ('Savig figure to file "%s" ...' % args.output_file)
        if(args.fig_resolution[0].upper() == 'H'):
            print (indent + '* High-resolution figure (600 dpi).')
            fig_resolution_in_dpi = 600
        elif (args.fig_resolution[0].upper() == 'M'):
            print (indent + '* Medium-resolution figure (300 dpi).')
            fig_resolution_in_dpi = 300
        elif (args.fig_resolution[0].upper() == 'L'):
            print (indent + '* Low-resolution figure (100 dpi).')
            fig_resolution_in_dpi = 100
        else:
            print (indent + 'Assuming medium-resolution (300 dpi) for the figure.')
            fig_resolution_in_dpi = 300
        plt.savefig(args.output_file, dpi=fig_resolution_in_dpi, bbox_inches='tight')
        print (indent + '* Done saving figure (%s).' % args.output_file)

    if args.saveshow:
        print ('Opening saved figure (%s)...' % default_out_basename)
        # 'xdg-open' might fail to find the defualt program in some systems
        # For such cases, one can try to use other alternatives (just add more to the list below)
        image_viewer_list = ['xdg-open', 'eog']
        for image_viewer in image_viewer_list:
            open_saved_fig = Popen([image_viewer, args.output_file], stdout=PIPE, stderr=PIPE)
            std_out, std_err = open_saved_fig.communicate()
            success_opening_file = std_err.strip() == ''
            if(success_opening_file):
                break
        if(not success_opening_file):
            print (indent + '* Failed (%s): no image viewer detected.' % default_out_basename)

    if args.show:
        print ('Showing figure (%s)...' % default_out_basename)
        plt.show()
        print (indent + '* Done showing figure (%s).' % default_out_basename)


if __name__ == '__main__':
    print_opening_message()
    plot_options = BandUpPlotOptions()
    plot = BandUpPlot(plot_options)
    make_plot(plot)
    sys.exit(0)
