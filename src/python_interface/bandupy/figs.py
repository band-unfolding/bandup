# Copyright (C) 2017 Paulo V. C. Medeiros
# Module used by the Python interface to BandUP to deal with figures
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
import os
import json
import matplotlib as mpl
# Below: To be able to create mpl figures without running an X server
#        mpl.use must be called before importing pyplot
mpl.use('Agg')
from matplotlib import pyplot as plt
# Imports from within the package
from .constants import PACKAGE_DIR


def get_matplotlib_color_names():
    return sorted(mpl.colors.cnames.keys())

def get_available_cmaps(only_names=False):
    colormap_names = sorted(plt.cm.datad.keys(), key=lambda s: s.lower())
    # Removing deprecated colormaps for future portability with Python 3
    # Matplotlib and warnings don't seem to like each other in python3 if these 
    # cmaps are not removed
    colormap_names = [name for name in colormap_names if 
                      not name.lower().startswith('spectral') and 
                      not 'Vega' in name]

    colormaps = dict([[cmap_name, plt.get_cmap(cmap_name)] for cmap_name in 
                     colormap_names])
    # Custom colormaps
    try:
        custom_cm_folder = os.path.join(PACKAGE_DIR, 'custom_colormaps')
        os_listdir_full_custom_cm_folder = (
            [os.path.abspath(os.path.join(custom_cm_folder, cmap_file)) for 
             cmap_file in os.listdir(custom_cm_folder)]
        )
        custom_cmap_files = [cmap_file for cmap_file in 
                             os_listdir_full_custom_cm_folder if 
                             cmap_file.endswith('.cmap')]
        custom_cmaps = []
        for cmap_file in custom_cmap_files:
            cmap_name = os.path.splitext(os.path.basename(cmap_file))[0]
            color_dict = json.load(open(cmap_file))
            colormap_names.append(cmap_name)
            colormaps[cmap_name] = (
                plt.cm.colors.LinearSegmentedColormap(cmap_name, color_dict, 2048)
            )
    except(TypeError):
        pass
    if(only_names):
        return sorted(colormap_names)
    else:
        return sorted(colormap_names), colormaps


def allowed_fig_formats():
    temp_fig = plt.figure()
    all_allowed_filetypes = [f.lower() for f in 
                             temp_fig.canvas.get_supported_filetypes().keys()]
    plt.close(temp_fig)
    preferred_filetypes = ['tiff', 'png', 'bmp', 'jpg', 'pdf', 'eps']
    allowed_filetypes = [fmt for fmt in preferred_filetypes if fmt in 
                         all_allowed_filetypes]
    allowed_filetypes += [fmt for fmt in all_allowed_filetypes if fmt not in 
                          preferred_filetypes]
    return allowed_filetypes


def set_default_fig_format(default_fig_format='tiff'):
    allowed_filetypes = allowed_fig_formats()
    if default_fig_format not in allowed_filetypes:
        default_fig_format = allowed_filetypes[0]
    return default_fig_format

def output_file_with_supported_extension(output_file_name, default_fig_format):
    """Defines the file format of the generated figure. 

        This routine checks whether the chosen file type is supported by the system, and,
        if it's not, returns the input file name with the original extension replaced by
         a valid one.
    """
    allowed_filetypes = allowed_fig_formats()
    file_name = splitext(basename(output_file_name))[0]
    file_extension = splitext(output_file_name)[1][1:]
    if file_extension.lower() not in allowed_filetypes:
        print ('WARNING: The file type you requested for the output figure (%s) is not '
               "supported by your system and has been changed to '%s'." 
               % (file_extension, default_fig_format))
        print ('         * The choices supported by your system are: %s' 
               % ", ".join(map(str,allowed_filetypes)))
        output_file_name = file_name + '.' + default_fig_format
    return output_file_name

