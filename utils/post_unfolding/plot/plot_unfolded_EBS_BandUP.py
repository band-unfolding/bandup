#! /usr/bin/env python
# Copyright (C) 2013, 2014 Paulo V. C. Medeiros, Jonas Bjork
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
import argparse
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import getcwd
from os.path import join, realpath, dirname, splitext
import sys
import MyColMaps 

# Defining the options passed in the command line
parser = argparse.ArgumentParser()
parser.add_argument('input_file', default='unfolded_EBS_symmetry-averaged.dat', nargs='?', 
                    help='Name of the input file.')
parser.add_argument('output_file', nargs='?', default=None, 
                    help='Optional: Name of the output file. If not given, it will be based on ' 
                         'the name of the input file.')
parser.add_argument('-kpts', '--kpoints_file', default='KPOINTS_prim_cell.in', 
                    help='Name of the file containing information'
                    'about the primitive cell k-points. Default: KPOINTS_prim_cell.in')
parser.add_argument('-pc_file', '--prim_cell_file', default='prim_cell_lattice.in', 
                    help='Name of the file containing information'
                    'about the primitive cell lattice vectors. Default: prim_cell_lattice.in')
parser.add_argument('-efile', '--energy_info_file', default='energy_info.in', 
                    help='Name of the file containing information about the energy grid and '
                         'Fermi energy to be used. Default: energy_info.in')
possible_cmap_choices = parser.add_mutually_exclusive_group()
possible_cmap_choices.add_argument('-cmap', '--colormap', default=None,
                                   help='Choice of colormap for the plots (name of the colormap).v'
                                   'You might want to try a few.')
possible_cmap_choices.add_argument('-icmap', '--icolormap', type=int, default=None,
                                   help='Choice of colormap for the plots (integer number).')
parser.add_argument('--save', action='store_true', default=False, 
                    help='Saves the figue to a file. Default: False')
parser.add_argument('--show', action='store_true', default=False, 
                    help='Shows the figue. Default: False if --save is selected, True otherwise.')
parser.add_argument('-res', '--fig_resolution', default='m', choices=['l','m','h'], 
                    help='Resolution of the figure: l = 100 dpi, m = 300 dpi, h = 600 dpi.'
                    'Default: m = 300 dpi')
parser.add_argument('--disable_auto_round_vmin_and_vmax', action='store_true', 
                    help='Disable normalization of the vmin and vmax of the color scale to ' 
                    'integer numbers.')
# Defining the default file format of the generated figure.
# If the default file format is not supported, it will be changed.
default_fig_format = 'png' # Change here if you want another format to be the default one. 
temp_fig = plt.figure()
allowed_filetypes = [format.lower() for format in temp_fig.canvas.get_supported_filetypes().keys()]
plt.close(temp_fig)
del temp_fig
if default_fig_format not in allowed_filetypes:
    default_fig_format = allowed_filetypes[0]
# Done defining the default file format of the generated figure.
parser.add_argument('-fmt', '--file_format', default=default_fig_format,
                    help='File format of the figure. Default: ' + default_fig_format)
possible_fig_orientations = parser.add_mutually_exclusive_group()
possible_fig_orientations.add_argument('--landscape', action='store_true')
possible_fig_orientations.add_argument('--portrait', action='store_true')
parser.add_argument('-ar', '--aspect_ratio', default='3.0/4.0', 
                    help='Aspect ratio of the generated plot. Default: 3/4')
parser.add_argument('--no_ef', action='store_true',help='Hides the E-Fermi line.')
parser.add_argument('--no_cb', action='store_true',help='Hides the colorbar.')
parser.add_argument('-nlev', '--n_levels', type=int, default=1001, 
                    help='Number of different levels used in the contour plot. Default: 1001')
parser.add_argument('-shift_k', '--shift_kpts_coords', type=float, default=0.0, 
                    help='Shift in the k-points. Default: 0.0')
parser.add_argument('-shift_e', '--shift_energy', type=float, default=0.0,
                    help='Shift in the energy grid. Default: 0.0')
parser.add_argument('-vmin', '--minval_for_colorbar', type=float, 
                    help='Value to which the first color of the colormap will be normalized.')
parser.add_argument('-vmax', '--maxval_for_colorbar', type=float, 
                    help='Value to which the last color of the colormap will be normalized.')
parser.add_argument('--round_cb', type=int, default=1, 
                    help='Number of decimal digits displayed in the colobar ticks.')
show_cb_label_options = parser.add_mutually_exclusive_group()
show_cb_label_options.add_argument('--cb_label', action='store_true', 
                                   help='Show the colorbar label.')
show_cb_label_options.add_argument('--cb_label_full', action='store_true', 
                                   help='Show the colorbar label (full).')
args = parser.parse_args()
# Opening message
print ('                                                                                      \n'
       '===================================================================================== \n'
       '             BandUP: Band Unfolding code for Plane-wave based calculations            \n'
       '                                                                                      \n'
       '              Post-processing utility "plot_unfolded_EBS_BandUP.py"                   \n'
       '            >>> Visualizing the unfolded EBSs produced by BandUP <<<                  \n'
       '===================================================================================== \n'
       'Copyright (C) 2013, 2014 Paulo V. C. Medeiros, Jonas Bjork                            \n'
       '                         paume@ifm.liu.se, jonas.bjork@liu.se                         \n'
       '                         Computational Physics Division                               \n'
       '                         Department of Physics, Chemistry and Biology - IFM           \n'
       '                         Linkoping University                                         \n'
       '                         Sweden                                                       \n'
       'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                           \n'
       '===================================================================================== \n'
       '                                                                                      \n')

indent='    '
###################################################################################################

title_of_the_plot = ''
x_axis_label = ''
y_axis_label = '$\epsilon \hspace{0.25} - \hspace{0.25} \epsilon _{_F} (eV)$'

# Defining the input file
input_file = args.input_file


scaling_factor_fig_size = 2.0 # 1 gives a width of 8.6cm for the final figure.
title_size = 10*scaling_factor_fig_size
yaxis_labels_size = 14*scaling_factor_fig_size
xaxis_labels_size = 12*scaling_factor_fig_size
tick_marks_size = 8*scaling_factor_fig_size
colorbar_label_size = 8*scaling_factor_fig_size
colorbar_tick_marks_size = 7*scaling_factor_fig_size
line_width_high_symm_points = 0.6*scaling_factor_fig_size
line_width_E_f = 0.8*scaling_factor_fig_size
# Position of the line y = E_F in the figure
E_f = 0.000
# Figure parameters
# Change the value of 'default_colormap' if you would like to have another colormap as default. 
# I used 'gist_ncar' in my paper. 
default_colormap = 'gist_ncar'
aspect_ratio = eval(args.aspect_ratio)
if args.portrait:
    orientation = 'portrait'
elif args.landscape:
    orientation = 'landscape'
else:
    orientation = 'portrait'

if orientation=='landscape':
    cb_orientation = 'vertical'
    offset_x_text_colorbar = 3.5
    offset_y_text_colorbar = 0.5
else:
    cb_orientation = 'horizontal'
    offset_x_text_colorbar = 0.5
    offset_y_text_colorbar = -2.5
    xaxis_labels_size = xaxis_labels_size/aspect_ratio
    yaxis_labels_size = yaxis_labels_size/aspect_ratio
    tick_marks_size = tick_marks_size/aspect_ratio
    colorbar_label_size = colorbar_label_size/aspect_ratio
    colorbar_tick_marks_size = colorbar_tick_marks_size/aspect_ratio
    line_width_high_symm_points = line_width_high_symm_points/aspect_ratio
    line_width_E_f = line_width_E_f/aspect_ratio


# Reading data from the input file
mypath = realpath(join(getcwd(), dirname(__file__)))
print 'Reading input file "%s"' % input_file
KptsCoords, energies, delta_Ns = np.loadtxt(join(mypath,input_file),usecols=(0,1,2),unpack=True)
print indent + '* Max. delta_N:          ',np.max(delta_Ns)
print indent + '* Min. non-zero delta_N: ',np.min(delta_Ns[np.nonzero(delta_Ns)])
print 'Eliminating duplicated points...'
sorted_index_data = np.lexsort((energies, KptsCoords))
KptsCoords = KptsCoords[sorted_index_data] + args.shift_kpts_coords
energies = energies[sorted_index_data] + args.shift_energy
delta_Ns = delta_Ns[sorted_index_data]

new_KptsCoords = []
new_energies = []
new_delta_Ns = []
for line_number in range(len(KptsCoords)-1):
    skip_line = False
    if(KptsCoords[line_number]==KptsCoords[line_number+1] and 
       energies[line_number]==energies[line_number+1]):
       delta_Ns[line_number+1] = max(delta_Ns[line_number],delta_Ns[line_number+1])
       skip_line = True
    if(not skip_line):
        new_energies.append(energies[line_number])
        new_KptsCoords.append(KptsCoords[line_number])
        new_delta_Ns.append(delta_Ns[line_number])
# Appending the last point
new_energies.append(energies[-1])
new_KptsCoords.append(KptsCoords[-1])
new_delta_Ns.append(delta_Ns[-1])

KptsCoords = np.array(new_KptsCoords)
energies = np.array(new_energies)
delta_Ns = np.array(new_delta_Ns)
print indent + '* Done with eliminating duplicated points.'

# Defining the boundaries of the plot
with open(args.energy_info_file) as energy_info_file:
    energy_info_file_lines = energy_info_file.readlines()
    ymin = float(energy_info_file_lines[1].split()[0])
    ymax = float(energy_info_file_lines[2].split()[0])
    if ymin > ymax:
        ymin, ymax = ymax, ymin
    try:
        xmin = float(energy_info_file_lines[1].split()[1])
    except:
        xmin = KptsCoords.min()
    try:
        xmax = float(energy_info_file_lines[2].split()[1])
    except:
        xmax = KptsCoords.max()
    if xmin > xmax:
        xmin, xmax = xmax, xmin

try:
    energy_tolerance_for_hist2d = abs(float(energy_info_file_lines[3].split()[0]))
    if energy_tolerance_for_hist2d == 0:
        raise ValueError
except:
    energy_tolerance_for_hist2d = (0.4E-2)*(ymax - ymin)
    print 'WARNING: Could not read the size of the energy intervals. Setting it to ' \
          '0.4E-2*(Emax - Emin) = ', energy_tolerance_for_hist2d, ' eV.'

size_of_old_data = float(len(KptsCoords))
if xmin>KptsCoords.min() or xmax<KptsCoords.max() or ymin>energies.min() or ymax<energies.max():
    print 'Trying to reduce data such that only the needed part of it is parsed.'
    new_KptsCoords = []
    new_energies = []
    new_delta_Ns = []
    for iline in range(len(KptsCoords)):
        if (KptsCoords[iline] >= xmin and KptsCoords[iline] <= xmax and 
            energies[iline] >= ymin and energies[iline] <= ymax):
            new_KptsCoords.append(KptsCoords[iline])
            new_energies.append(energies[iline])
            new_delta_Ns.append(delta_Ns[iline])
    KptsCoords = np.array(new_KptsCoords)
    energies = np.array(new_energies)
    delta_Ns = np.array(new_delta_Ns)
    size_of_new_data = float(len(KptsCoords))
    print indent + '* Done. Working with %.2f' % (100.0*size_of_new_data/size_of_old_data), \
          '% of the data read in.'

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
            print 'ERROR: You have selected cartesian coordinates in your input k-points file, ' \
                  'but you have not passed a scaling parameter "a0".'
            print '       The actuall coordiates of the k-points are given by: '\
                  'ki[actual] = two_pi*ki[passed in file]/a0.'
            print '       Please write the value of a0 after your tag "' + \
                  coords_type + '", and run the code again.'
            print 'Stopping now.'
            sys.exit(0)
        k_start[idir] = (2.0*np.pi)*read_k_start[idir]/latt_param_kpts_file
        k_end[idir] = (2.0*np.pi)*read_k_end[idir]/latt_param_kpts_file
    else:
        if((coords_type[0].upper() != 'R') and (idir==0)):
            print 'WARNING: Assuming that the pc-kpts have been informed in '\
                  'fractional (reciprocal) coordinates.'
        for i in range(3):
            k_start[idir] += read_k_start[idir,i]*b_matrix_pc[i]
            k_end[idir] += read_k_end[idir,i]*b_matrix_pc[i]


pos_high_symm_lines = [zero_of_kpts_line]
for idir in range(0,ndirections):
    pos_high_symm_lines.append(pos_high_symm_lines[-1] + 
                               np.linalg.norm(k_end[idir] - k_start[idir]))
print 'Vertical lines will be automatically drawn at: k = %s' \
      % ', '.join(map("{:9.5f}".format,pos_high_symm_lines))
# End of determining the positions of high-symmetry BZ points on the plot
# Labeling the high-symmetry BZ points
labels_high_symm_lines = [label_k_start[0]]
for idir in range(1,ndirections):
    if (label_k_start[idir] == label_k_end[idir-1]):
        labels_high_symm_lines.append(label_k_start[idir])
    else:
        labels_high_symm_lines.append(label_k_end[idir-1]+','+label_k_start[idir])
labels_high_symm_lines += [label_k_end[-1]]
def convert_to_symbols(list):
    dict = {'G':'$\Gamma$', 'GAMMA':'$\Gamma$','DELTA':'$\Delta$',
            'LAMBDA':'$\Lambda$','SIGMA':'$\Sigma$'}
    for i in range(len(list)):
        if dict.has_key(list[i].upper()):
            list[i] = dict[list[i].upper()]
convert_to_symbols(labels_high_symm_lines)

std_column_widths_cm = [8.6, 17.8]
std_column_widths_inch = [item/2.54 for item in std_column_widths_cm]
fig_widths_options_in_inches = [std_column_widths_inch[0]*scaling_factor_fig_size, 
                                5.0*std_column_widths_inch[0], std_column_widths_inch[1], 18.5]
fig_width_inches = fig_widths_options_in_inches[0]
if aspect_ratio>1:
    aspect_ratio = 1.0/aspect_ratio
if orientation == 'landscape':
    fig_height_inches = fig_width_inches*aspect_ratio
else:
    if orientation != 'portrait':
        print 'Assuming portrait orientation for the output.'
    fig_height_inches = fig_width_inches/aspect_ratio

# Creating the plot
print 'Generating the plot...'
fig = plt.figure(figsize=(fig_width_inches,fig_height_inches))
ax = fig.add_subplot(111)

# Defining the color schemes.
# Custom colormaps
try: 
    custom_cm_flame = MyColMaps.LoadColMap('Flame.cmap','') # It's still a test
except:
    custom_cm_flame = None
all_custom_cmaps = [custom_cm_flame]
# Matplotlib colormaps
mpl_colormaps_I_tested = [plt.cm.gist_ncar,plt.cm.jet, plt.cm.jet_r, 
                      plt.cm.seismic, plt.cm.seismic_r, plt.cm.hot, 
                      plt.cm.afmhot, plt.cm.gist_heat, plt.cm.Greys_r, 
                      plt.cm.cubehelix_r, plt.cm.gnuplot2_r, plt.cm.terrain_r, 
                      plt.cm.gist_stern_r,plt.cm.hot_r]
all_other_mpl_cmaps = [plt.get_cmap(cmap) for cmap in plt.cm.datad if not 
                       plt.get_cmap(cmap) in mpl_colormaps_I_tested] 
valid_custom_cmaps = [cmap for cmap in all_custom_cmaps if cmap is not None]
# Gathering all colormaps
colormaps = mpl_colormaps_I_tested + all_other_mpl_cmaps + valid_custom_cmaps

valid_cmap_choice = True
user_did_not_specify_cmap = (args.icolormap == None and args.colormap == None)
if(user_did_not_specify_cmap):
    chosen_colormap = plt.get_cmap(default_colormap)
else:
    if(args.icolormap != None):
        if(args.icolormap<=len(colormaps)-1):
            chosen_colormap = colormaps[args.icolormap]
        else:
            valid_cmap_choice = False
            print indent + 'WARNING: icolormap options range from 0 to ', len(colormaps)-1,'.'
    if(args.colormap != None):
        for cmap in plt.cm.datad:
            if cmap.upper() == args.colormap.upper():
                args.colormap = cmap
        if(plt.get_cmap(args.colormap) in colormaps):
            chosen_colormap = plt.get_cmap(args.colormap)
        else:
            valid_cmap_choice = False
            print indent + 'WARNING: "' + args.colormap + '" is not a valid choice of colormap.'
    if(not valid_cmap_choice):
        chosen_colormap = plt.get_cmap(default_colormap)
        print 3 * indent + ' * The default color scheme will be used.'
for cmap in plt.cm.datad:
    if plt.get_cmap(cmap) == chosen_colormap:
        print indent + '>>> Using the "' + cmap + '" colormap.'
        break
if(user_did_not_specify_cmap):
    print 2 * indent + 'Tip: You can try different colormaps by either:'
    print 2 * indent + '     * Running the plot tool with the option -icmap n, ' \
          'with n in the range from 0 to', len(colormaps)-1
    print 2 * indent + '     * Running the plot tool with the option -cmap cmap_name.'
    print 2 * indent + '> Take a look at ' \
          '<http://matplotlib.org/examples/color/colormaps_reference.html> for '\
          'a list of colormaps.'

# Building the countour plot from the read data
# Defining the (ki,Ej) grid.
ki = np.linspace(xmin,xmax,2*len(set(KptsCoords))+1,endpoint=True)
Ei = np.arange(ymin,ymax + energy_tolerance_for_hist2d,energy_tolerance_for_hist2d)
# Interpolating
grid_freq = griddata((KptsCoords, energies), delta_Ns, (ki[None,:], Ei[:,None]), 
                     method='linear',fill_value=0.0)
grid_freq = grid_freq.clip(0.0) # Values smaller than zero are just noise.
# Normalizing and building the countour plot
n_levels = args.n_levels
manually_normalize_colorbar_min_and_maxval = False
if(args.maxval_for_colorbar != None or args.minval_for_colorbar != None):
    manually_normalize_colorbar_min_and_maxval = True
    args.disable_auto_round_vmin_and_vmax = True
    maxval_for_colorbar = args.maxval_for_colorbar
    minval_for_colorbar = args.minval_for_colorbar
else:
    if not args.disable_auto_round_vmin_and_vmax:
        minval_for_colorbar = round(np.min(grid_freq))
        maxval_for_colorbar = round(np.max(grid_freq))
        args.round_cb = 0
if manually_normalize_colorbar_min_and_maxval or not args.disable_auto_round_vmin_and_vmax:
    if not args.disable_auto_round_vmin_and_vmax:
        print indent + '* Automatically renormalizing color scale '\
              '(you can disable this with the option --disable_auto_round_vmin_and_vmax):'
    if manually_normalize_colorbar_min_and_maxval:
        print indent + '* Manually renormalizing color scale:'
    if(minval_for_colorbar != None):
        print 2 * indent + 'Previous vmin = %.1f, new vmin = %.1f' % (np.min(grid_freq), 
                                                                      minval_for_colorbar)
    else:
        minval_for_colorbar = np.min(grid_freq)
    if(maxval_for_colorbar != None):
        print 2 * indent + 'Previous vmax = %.1f, new vmax = %.1f' % (np.max(grid_freq), 
                                                                      maxval_for_colorbar)
    else:
        maxval_for_colorbar = np.max(grid_freq)
    print (2 * indent + 'The previous vmin and vmax might be slightly different from '
                        'the min and max delta_Ns '
                        'due to the interpolation scheme used for the plot.')
    # values > vmax will be set to vmax, and #<vmin will be set to vmin 
    grid_freq = grid_freq.clip(minval_for_colorbar, maxval_for_colorbar)
 
    v = np.linspace(minval_for_colorbar, maxval_for_colorbar, n_levels, endpoint=True)
    image = ax.contourf(ki,Ei,grid_freq,levels=v,cmap=chosen_colormap)
else: 
    image = ax.contourf(ki,Ei,grid_freq,n_levels,cmap=chosen_colormap)

#Preparing the plot

# Defining the colors of the vertical and horizontal lines.
# Don't take this part too seriously. Change it if you need/want.
dist_black = np.linalg.norm(np.array(chosen_colormap(image.norm.vmin)[:3]))
dist_white = np.linalg.norm(np.array(chosen_colormap(image.norm.vmin)[:3]) - np.array([1,1,1]))
if(dist_black < dist_white):
    color_v_and_h_lines = 'white'
else:
    color_v_and_h_lines = 'black'

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_title(title_of_the_plot, fontsize=title_size)
ax.set_ylabel(y_axis_label, fontsize=yaxis_labels_size)
plt.yticks(fontsize=tick_marks_size)

# Fermi energy line
show_E_f = not args.no_ef
if show_E_f and E_f>=ymin and E_f<=ymax:
    plt.axhline(y=E_f, c=color_v_and_h_lines, linestyle=':', lw=line_width_E_f)
# High symmetry points lines
x_tiks_positions = [kx for kx in pos_high_symm_lines if kx-xmax<=1E-2 and kx >= xmin]
x_tiks_labels = [labels_high_symm_lines[i] for i in range(len(labels_high_symm_lines)) if 
                 pos_high_symm_lines[i] in x_tiks_positions]

x_tiks_labels = [xlabel for xlabel in x_tiks_labels if xlabel]
if x_tiks_labels:
    print indent + '* K-point labels read from the "' + args.kpoints_file + '" file:'
    for ilabel in range(len(x_tiks_labels)):
        print 2 * indent + "k = {:9.5f}".format(x_tiks_positions[ilabel]) + ', label =',\
              x_tiks_labels[ilabel]
    plt.xticks(x_tiks_positions, x_tiks_labels, fontsize=tick_marks_size)
else:
    x_axis_label = '$k \hspace{0.25} (\AA^{-1})$'
    plt.locator_params(axis = 'x', nbins = 5)
    ax.set_xlabel(x_axis_label, fontsize=xaxis_labels_size)
    plt.xticks(fontsize=tick_marks_size)
ax.tick_params(axis='x', pad=10)
for line_position in [pos for pos in pos_high_symm_lines if pos > xmin and pos < xmax]:
    plt.axvline(x=line_position, c=color_v_and_h_lines, linestyle='-', 
                lw=line_width_high_symm_points)
# Color bar
show_colorbar = not args.no_cb
if show_colorbar:
    if cb_orientation=='vertical':
        cb_pad=0.005
    else:
        cb_pad=0.06

    cb_yticks = np.arange(int(image.norm.vmin), int(image.norm.vmax) + 1, 1)
    def round(f,n):
                return '%.*f' % (n, f)
    cb_ytick_labels = [round(item,abs(args.round_cb)) for item in cb_yticks]
    cb = plt.colorbar(image, ax=ax, ticks=cb_yticks, orientation=cb_orientation, pad=cb_pad)
    cb.set_ticklabels(cb_ytick_labels)
    cb.ax.tick_params(labelsize=colorbar_tick_marks_size)
   
    color_bar_label = None
    if args.cb_label: 
        color_bar_label = ('$Color scale: \hspace{0.5} \delta N(\\vec{k}; ' + 
                           '\hspace{0.25} \epsilon)$ ') 
    if args.cb_label_full: 
        color_bar_label = ('$Colors cale: \hspace{0.5} \delta N(\\vec{k}; ' +
                           '\hspace{0.25} \epsilon);$ '+ 
                           '$\delta\epsilon='+round(1000.0*energy_tolerance_for_hist2d,0) +
                           '\\hspace{0.25} meV.$')

    if cb_orientation=='vertical':
        cb_label_rotation = 90
    else:
        cb_label_rotation = 0
    if color_bar_label:
        cb.ax.text(offset_x_text_colorbar,offset_y_text_colorbar,
                   color_bar_label,rotation=cb_label_rotation,ha='center',
                   va='center',fontsize=colorbar_label_size)
# Showing the results
plt.tick_params(\
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    left='off',     
    right='off',         
    labelbottom='on')

def validated_output_file_format(output_file_name,allowed_filetypes,default_fig_format):
    # Defines the file format of the generated figure according to what your system supports.
    file_name = splitext(output_file_name)[0]
    file_extension = splitext(output_file_name)[1][1:]
    if file_extension.lower() not in allowed_filetypes:
        print ('WARNING: The file type you requested for the output figure (%s) is not '
               "supported by your system and has been changed to '%s'." 
               % (file_extension, default_fig_format))
        print ('         * The choices supported by your system are: %s' 
               % ", ".join(map(str,allowed_filetypes)))
        output_file_name = file_name + '.' + default_fig_format
    return output_file_name

output_file_base = "_".join([splitext(input_file)[0], 'E_from', str(ymin), 'to', 
                            str(ymax), 'eV_dE', 
                            str(energy_tolerance_for_hist2d), 'eV'])
if(not args.save and not args.show):
    args.save = False
    args.show = True
if args.output_file:
    args.save = True
    output_file_name = args.output_file
else:
    output_file_name = output_file_base + '.' + args.file_format
if (args.save):
    output_file_name = validated_output_file_format(output_file_name,
                                                    allowed_filetypes,
                                                    default_fig_format)
    print 'Savig figure to file "%s" ...' % output_file_name
    fig_resolution = args.fig_resolution
    if(fig_resolution[0].upper() == 'H'):
        print indent + '* High-resolution figure (600 dpi).'
        fig_resolution_in_dpi = 600
    elif (fig_resolution[0].upper() == 'M'):
        print indent + '* Medium-resolution figure (300 dpi).'
        fig_resolution_in_dpi = 300
    elif (fig_resolution[0].upper() == 'L'):
        print indent + '* Low-resolution figure (100 dpi).'
        fig_resolution_in_dpi = 100
    else:
        print indent + 'Assuming medium-resolution (300 dpi) for the figure.'
        fig_resolution_in_dpi = 300
    plt.savefig(output_file_name, dpi=fig_resolution_in_dpi, bbox_inches='tight')
    print indent + '* Done saving figure (%s).' % output_file_name

if args.show:
    print 'Showing figure (%s)...' % output_file_base
    plt.show()
    print indent + '* Done showing figure (%s).' % output_file_base
