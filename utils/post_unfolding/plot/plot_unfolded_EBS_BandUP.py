#! /usr/bin/env python
# Copyright (C) 2013 Paulo V. C. Medeiros, Jonas Bjork
# This file is part of the BandUP code: Band Unfolding code for Plane-wave based calculations.
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
from matplotlib.colors import LinearSegmentedColormap
from os import getcwd
from os.path import join, realpath, dirname
import sys

parser = argparse.ArgumentParser()
parser.add_argument('input_file', default='unfolded_band_structure.dat', nargs='?', help='')
parser.add_argument('output_file', nargs='?', default=None, help='')
parser.add_argument('-kpts', '--kpoints_file', default='KPOINTS_prim_cell.in')
parser.add_argument('-efile', '--energy_info_file', default='energy_info.in')
parser.add_argument('-shift_k', '--shift_kpts_coords', type=float, help='', default=0.0)
parser.add_argument('-shift_e', '--shift_energy', type=float, help='', default=0.0)
parser.add_argument('--no_ef', help='Hides the E-Fermi line.', action='store_true')
parser.add_argument('--no_cb', help='Hides the colorbar.', action='store_true')
parser.add_argument('-vmax', '--maxval_for_colorbar', type=float, help='Value the last color of the colormap will be normalized to.')
parser.add_argument('-vmin', '--minval_for_colorbar', type=float, help='Value the first color of the colormap will be normalized to.')
parser.add_argument('-nlev', '--n_levels', type=int, default=1001, help='Number of different levels used.')
parser.add_argument('-ar', '--aspect_ratio', default='3.0/4.0', help='Aspect ratio for the generated plot.')
parser.add_argument('-res', '--fig_resolution', default='m', choices=['l','m','h'], help='Resolution for the generated plot: l = 100 dpi, m = 300 dpi, h = 600 dpi.')
parser.add_argument('-icmap', '--icolormap', type=int, default=None)
possible_fig_orientations = parser.add_mutually_exclusive_group()
possible_fig_orientations.add_argument('--landscape', action='store_true')
possible_fig_orientations.add_argument('--portrait', action='store_true')
parser.add_argument('-fmt', '--file_format', default='png', choices=['png', 'eps', 'pdf', 'jpg', 'tif'])
parser.add_argument('--round_cb', type=int, default=1, help='Number of digits displayed in the colobar ticks.')
show_cb_label_options = parser.add_mutually_exclusive_group()
show_cb_label_options.add_argument('--cb_label', action='store_true', help='Show the colorbar label.')
show_cb_label_options.add_argument('--cb_label_full', action='store_true', help='Show the colorbar label (full).')
parser.add_argument('--save', help='Saves the figue to a file.', action='store_true', default=False)
parser.add_argument('--show', help='Shows the figue.', action='store_true', default=False)
args = parser.parse_args()

print ('===================================================================================== \n'
       '             BandUP: Band Unfolding code for Plane-wave based calculations            \n'
       '===================================================================================== \n'
       'Copyright (C) 2013, 2014 Paulo V. C. Medeiros, Jonas Bjork                            \n'
       '                         paume@ifm.liu.se, jonas.bjork@liu.se                         \n'
       '                         Computational Physics Division                               \n'
       '                         Department of Physics, Chemistry and Biology - IFM           \n'
       '                         Linkoping University                                         \n'
       '                         Sweden                                                       \n'
       'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                           \n'
       '===================================================================================== \n'
       '                                                                                      \n'
       '              Post-processing utility "plot_unfolded_EBS_BandUP.py"                   \n'
       '            >>> Visualizing the unfolded EBSs produced by BandUP <<<                  \n'
       '                                                                                      \n')

indent='    '
############################################################################################################

title_of_the_plot = ''
x_axis_label = ''
y_axis_label = '$\epsilon \hspace{0.25} - \hspace{0.25} \epsilon _{_F} (eV)$'

# Defining input and output files
myfile = args.input_file


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


# Reading data from input file
mypath = realpath(join(getcwd(), dirname(__file__)))
print 'Reading input file "%s"' % myfile
KptsCoords, energies, delta_Ns = np.loadtxt(join(mypath,myfile),usecols=(0,1,2),unpack=True)
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
    energy_tolerance_for_hist2d = float(energy_info_file_lines[3].split()[0])  # 0.175 seemed fine.
except:
    energy_tolerance_for_hist2d = (0.4E-2)*(ymax - ymin)
    print 'Automatically setting the size of the energy intervals to 0.4E-2*(Emax - Emin) = ', energy_tolerance_for_hist2d,' eV.'

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
print indent + '* Done. Working with %.2f' % (100.0*size_of_new_data/size_of_old_data), '% of the data read in.' 

# Determining the positions of high-symmetry BZ points on the plot
kpts_file_lines = []
with open(args.kpoints_file,'r') as kpts_file:
    for line in kpts_file:
        if line.strip(): 
            # Skipping blank lines
            kpts_file_lines.append(line)
latt_param_kpts_file = np.float(kpts_file_lines[0].split()[0])
try:
    zero_of_kpts_line = np.float(kpts_file_lines[2].split()[1])
except:
    zero_of_kpts_line = 0.0

coords_type = kpts_file_lines[3].strip()

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

with open('prim_cell_lattice.in') as prim_cell_lattice_file:
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
    if(coords_type[0].upper() == 'R'):
        for i in range(3):
            k_start[idir] += read_k_start[idir,i]*b_matrix_pc[i]
            k_end[idir] += read_k_end[idir,i]*b_matrix_pc[i]
    else:
        if((coords_type[0].upper() != 'C') and (idir==0)):
            print 'WARNING: Assuming that the pc-kpts have been informed in cartesian coordinates.'
        k_start[idir] = (2.0*np.pi)*read_k_start[idir]/latt_param_kpts_file
        k_end[idir] = (2.0*np.pi)*read_k_end[idir]/latt_param_kpts_file

pos_high_symm_lines = [zero_of_kpts_line]
for idir in range(0,ndirections):
    pos_high_symm_lines.append(pos_high_symm_lines[-1] + np.linalg.norm(k_end[idir] - k_start[idir]))
print 'Positions of the high-symmetry lines in the  k-axis: ', pos_high_symm_lines
# End of determining the positions of high-symmetry BZ points on the plot
# Labeling the high-symmetry BZ points
labels_high_symm_lines = [label_k_start[0]]
for idir in range(1,ndirections):
    if (label_k_start[idir] == label_k_end[idir-1]):
        labels_high_symm_lines.append(label_k_start[idir])
    else:
        labels_high_symm_lines.append(label_k_end[idir-1]+', '+label_k_start[idir])
labels_high_symm_lines += [label_k_end[-1]]
def convert_to_symbols(list):
    dict = {'G':'$\Gamma$', 'GAMMA':'$\Gamma$'}
    for i in range(len(list)):
        if dict.has_key(list[i].upper()):
            list[i] = dict[list[i].upper()]
convert_to_symbols(labels_high_symm_lines)

std_column_widths_cm = [8.6, 17.8]
std_column_widths_inch = [item/2.54 for item in std_column_widths_cm]
fig_widths_options_in_inches = [std_column_widths_inch[0]*scaling_factor_fig_size, 5.0*std_column_widths_inch[0], std_column_widths_inch[1], 18.5]
fig_width_inches = fig_widths_options_in_inches[0]
if aspect_ratio>1:
    aspect_ratio = 1.0/aspect_ratio
if orientation == 'landscape':
    fig_height_inches = fig_width_inches*aspect_ratio
else:
    if orientation != 'portrait':
        print 'Assuming portrait orientation for the output.'
    fig_height_inches = fig_width_inches/aspect_ratio
fig = plt.figure(figsize=(fig_width_inches,fig_height_inches))
ax = fig.add_subplot(111)

# Color schemes
colormaps = [plt.cm.jet, plt.cm.jet_r, plt.cm.gist_ncar, 
             plt.cm.seismic, plt.cm.seismic_r,
             plt.cm.hot, plt.cm.afmhot, plt.cm.gist_heat, 
             plt.cm.Greys_r]

chosen_colormap = plt.cm.gist_ncar
if(args.icolormap != None):
    colormap_option = args.icolormap
    if(colormap_option<=len(colormaps)-1):
        chosen_colormap = colormaps[colormap_option]
    else:
        print 'WARNING: Colormap option ranges from 0 to ', len(colormaps)-1,'.'
        print indent + 'Using the default colormap (gist_ncar).'

# Building the countour plot from the read data
# define grid.
print 'Generating the plot...'
ki = np.linspace(xmin,xmax,2*len(set(KptsCoords))+1,endpoint=True)
Ei = np.arange(ymin,ymax+energy_tolerance_for_hist2d,energy_tolerance_for_hist2d)
grid_freq = griddata((KptsCoords, energies), delta_Ns, (ki[None,:], Ei[:,None]), method='cubic',fill_value=0.0)
grid_freq = grid_freq.clip(0.0) # Values smaller than zero are just noise.

n_levels = args.n_levels
manually_normalize_colorbar_min_and_maxval = False
if(args.maxval_for_colorbar != None or args.minval_for_colorbar != None):
    manually_normalize_colorbar_min_and_maxval = True
    maxval_for_colorbar = args.maxval_for_colorbar
    minval_for_colorbar = args.minval_for_colorbar
if manually_normalize_colorbar_min_and_maxval:
    print indent + '* Renormalizing color scale:'
    if(minval_for_colorbar != None):
        print 2 * indent + 'Previous vmin = %.1f, new vmin = %.1f' % (np.min(grid_freq), minval_for_colorbar)
    else:
        minval_for_colorbar = np.min(grid_freq)
    if(maxval_for_colorbar != None):
        print 2 * indent + 'Previous vmax = %.1f, new vmax = %.1f' % (np.max(grid_freq), maxval_for_colorbar)
    else:
        maxval_for_colorbar = np.max(grid_freq)

    grid_freq = grid_freq.clip(minval_for_colorbar, maxval_for_colorbar) # #>vmax will be set to vmax, and #<vmin will be set to vmin 
    v = np.linspace(minval_for_colorbar, maxval_for_colorbar, n_levels, endpoint=True)
    image = ax.contourf(ki,Ei,grid_freq,levels=v,cmap=chosen_colormap)
else: 
    image = ax.contourf(ki,Ei,grid_freq,n_levels,cmap=chosen_colormap)

#Preparing the plot
if chosen_colormap(image.norm.vmin) == (1.0,1.0,1.0,1.0):
    color_v_and_h_lines = chosen_colormap(image.norm.vmax)
else:
    color_v_and_h_lines = (1.0,1.0,1.0,1.0)  # White color for the lines drawn
    
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
x_tiks_labels = [labels_high_symm_lines[i] for i in range(len(labels_high_symm_lines)) if pos_high_symm_lines[i] in x_tiks_positions]

x_tiks_labels = [xlabel for xlabel in x_tiks_labels if xlabel]
if x_tiks_labels:
    print indent + '* x-ticks placed at ', x_tiks_positions
    print indent + 'Labels:', x_tiks_labels
    plt.xticks(x_tiks_positions, x_tiks_labels, fontsize=tick_marks_size)
else:
    x_axis_label = '$k \hspace{0.25} (\AA^{-1})$'
    plt.locator_params(axis = 'x', nbins = 5)
    ax.set_xlabel(x_axis_label, fontsize=xaxis_labels_size)
    plt.xticks(fontsize=tick_marks_size)
ax.tick_params(axis='x', pad=10)
for line_position in pos_high_symm_lines:
    plt.axvline(x=line_position, c=color_v_and_h_lines, linestyle='-', lw=line_width_high_symm_points)
# Color bar
show_colorbar = not args.no_cb
if show_colorbar:
    if cb_orientation=='vertical':
        cb_pad=0.005
    else:
        cb_pad=0.05
 
    cb_yticks = [image.norm.vmin, 0.5*(image.norm.vmin+image.norm.vmax), image.norm.vmax]
    def round(f,n):
                return '%.*f' % (n, f)
    cb_ytick_labels = [round(item,abs(args.round_cb)) for item in cb_yticks]
    cb = plt.colorbar(image, ax=ax, ticks=cb_yticks, orientation=cb_orientation, pad=cb_pad)
    cb.set_ticklabels(cb_ytick_labels)
    cb.ax.tick_params(labelsize=colorbar_tick_marks_size)
   
    color_bar_label = None
    if args.cb_label: 
        color_bar_label = ('$Color scale: \hspace{0.5} \delta N(\\vec{k}; \hspace{0.25} \epsilon)$ ') 
    if args.cb_label_full: 
        color_bar_label = ('$Colors cale: \hspace{0.5} \delta N(\\vec{k}; \hspace{0.25} \epsilon);$ '+ 
                          '$\delta\epsilon='+round(0,1000.0*energy_tolerance_for_hist2d)+'\\hspace{0.25} meV.$')

    if cb_orientation=='vertical':
        cb_label_rotation = 90
    else:
        cb_label_rotation = 0
    if color_bar_label:
        cb.ax.text(offset_x_text_colorbar,offset_y_text_colorbar,color_bar_label,rotation=cb_label_rotation,ha='center',va='center',fontsize=colorbar_label_size)
# Showing the results
plt.tick_params(\
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    left='off',     
    right='off',         
    labelbottom='on')
#plt.tight_layout()

#print fig.canvas.get_supported_filetypes()
if(not args.save and not args.show):
    args.save = False
    args.show = True

if args.output_file:
    args.save = True
    output_file_name = args.output_file
else:
    output_file_name = myfile.strip('.dat') + '_E_from_' + str(ymin) + '_to_' + str(ymax) + '_eV_dE_' + str(energy_tolerance_for_hist2d) + '_eV' + '.'+args.file_format
if (args.save):
    fig_resolution = args.fig_resolution
    if(fig_resolution[0].upper() == 'H'):
        print indent + 'High-resolution figure (600 dpi) requested.'
        fig_resolution_in_dpi = 600
    elif (fig_resolution[0].upper() == 'M'):
        print indent + 'Medium-resolution figure (300 dpi) requested.'
        fig_resolution_in_dpi = 300
    elif (fig_resolution[0].upper() == 'L'):
        print indent + 'Low-resolution figure (100 dpi) requested.'
        fig_resolution_in_dpi = 100
    else:
        print indent + 'Assuming medium-resolution (300 dpi) for the figure.'
        fig_resolution_in_dpi = 300
    print indent + 'Savig figure to file "%s" ...' % output_file_name
    plt.savefig(output_file_name, dpi=fig_resolution_in_dpi, bbox_inches='tight')
    print indent + ' * Done.'

if args.show:
    print indent + 'Showing figure...'
    plt.show()
    print 'Done.'
