# Copyright (C) 2017 Paulo V. C. Medeiros
# A python wrapper to BandUP and its plotting tool
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
import os
from fractions import Fraction
from collections import OrderedDict
import sys
# Imports from within the package
from .environ import (
    bandup_dir,
    working_dir,
)
from .defaults import defaults
from .files import continuation_lines, get_efermi
from .figs import (
    get_matplotlib_color_names,
    get_available_cmaps,
    set_default_fig_format,
)

def get_bandup_registered_clas(bandup_path=bandup_dir):
    cla_file = os.path.join(bandup_path, 'src', 'cla_wrappers_mod.f90')
    cla_register_args = ['key', 'description', 'kkind', 'default']
    bandup_registered_clas = []
    for line in continuation_lines(cla_file):
        line = ' '.join(line.split())
        line = line.replace('// ','//').replace(' //','//')
        line = line.replace('//"','').replace('"//','')
        line = line.replace("//'",'').replace("'//",'')
        line = line.replace("new_line('A')",'\n').replace('new_line("A")','\n')
        arg_definitions_in_line = {}
        if("cla_register" in line):
            for cla_register_arg in cla_register_args:
                try:
                    cla_arg_start = line.index('%s='%(cla_register_arg))
                    cla_arg_end = min([line.index('%s='%(arg)) for
                                       arg in cla_register_args if 
                                       arg!=cla_register_arg and
                                       line.index('%s='%(arg))>cla_arg_start])
                except(ValueError):
                    cla_arg_end = len(line) - 1
                arg_k_arg_v = line[cla_arg_start:cla_arg_end].strip().strip(',')
                arg_k, arg_v = arg_k_arg_v.split('=')
                arg_v = arg_v.strip().strip('"').strip("'")
                arg_definitions_in_line[arg_k] = arg_v
            bandup_registered_clas.append(arg_definitions_in_line)
    return bandup_registered_clas


class BandUpPreUnfoldingArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if('remove_conflicts_with_bandup' in kwargs): 
            del kwargs['remove_conflicts_with_bandup']
        if('formatter_class') not in kwargs:
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        if('add_help' not in kwargs): kwargs['add_help'] = False
        super(BandUpPreUnfoldingArgumentParser, self).__init__(*args, **kwargs)

        # Adding BandUP's supported args to the parser
        symm_args = self.add_argument_group('Symmetry-related options')
        input_files_args = self.add_argument_group('Input file options')
        output_files_args = self.add_argument_group('Output file options')
        self.bandup_registered_clas = get_bandup_registered_clas()
        for cla_dict in self.bandup_registered_clas:
            arg_name = '%s'%(cla_dict['key'])
            arg_help = cla_dict['description']
            if(not arg_help): arg_help = 'No description available.'
            parser_remaining_kwargs = {}
            if(cla_dict['kkind'] in 
               ['cla_int', 'cla_float', 'cla_char', 'cla_xchar', 'cla_logical']):
                parser_remaining_kwargs['default'] = cla_dict['default']
            else:
                parser_remaining_kwargs['action'] = 'store_true'
            if(arg_name in ['-out_sckpts_file']):
                output_files_args.add_argument(arg_name, help=arg_help, 
                                               **parser_remaining_kwargs)
            elif(('no_symm' in arg_name) or ('skip_propose_pc' in arg_name)):
                symm_args.add_argument(arg_name, help=arg_help,
                                          **parser_remaining_kwargs)
            elif(arg_name in ['-pc_file', '-sc_file', '-pckpts_file']):
                input_files_args.add_argument(arg_name, help=arg_help,
                                         **parser_remaining_kwargs)

class BandUpArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if('remove_conflicts_with_bandup' in kwargs): 
            del kwargs['remove_conflicts_with_bandup']
        if('formatter_class') not in kwargs:
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        if('add_help' not in kwargs): kwargs['add_help'] = False
        super(BandUpArgumentParser, self).__init__(*args, **kwargs)

        # Adding BandUP's supported args to the parser
        debug_args = self.add_argument_group('Ideally only for tests')
        spinor_args = self.add_argument_group('For noncollinear mag. calcs.')
        symm_args = self.add_argument_group('Symmetry-related options')
        egrid_args = self.add_argument_group('Energy grid options '+
            '(specify either ONLY "-energy_file" OR *all* the others)')
        input_files_args = self.add_argument_group('Input file options')
        output_files_args = self.add_argument_group('Output file options')
        spin_args = self.add_argument_group('Spin-polarized calculations')
        castep_args = self.add_argument_group('For CASTEP')
        abinit_args = self.add_argument_group('For ABINIT')
        qe_args = self.add_argument_group('For Quantum ESPRESSO')
        self.bandup_registered_clas = get_bandup_registered_clas()
        for cla_dict in self.bandup_registered_clas:
            arg_name = '%s'%(cla_dict['key'])
            arg_help = cla_dict['description']
            if(not arg_help): arg_help = 'No description available.'
            parser_remaining_kwargs = {}
            if(cla_dict['kkind'] in 
               ['cla_int', 'cla_float', 'cla_char', 'cla_xchar', 'cla_logical']):
                parser_remaining_kwargs['default'] = cla_dict['default']
            else:
                parser_remaining_kwargs['action'] = 'store_true'
            if(arg_name in ['-saxis', '-normal_to_proj_plane', 
                              '-origin_for_spin_proj_rec', 
                              '-origin_for_spin_proj_cart']):
                spinor_args.add_argument(arg_name, help=arg_help, 
                                         **parser_remaining_kwargs)
            elif(arg_name in ['-castep', '-seed']):
                castep_args.add_argument(arg_name, help=arg_help, 
                                         **parser_remaining_kwargs)
            elif(arg_name in ['-abinit', '-files_file']):
                abinit_args.add_argument(arg_name, help=arg_help,
                                         **parser_remaining_kwargs)
            elif(arg_name in ['-qe', '-outdir', '-prefix']):
                qe_args.add_argument(arg_name, help=arg_help, 
                                         **parser_remaining_kwargs)
            elif(arg_name in ['-energy_file']):
                arg_help += ' Otherwise, see options below in this group.'
                egrid_args.add_argument(arg_name, help=arg_help, metavar='FILE',
                                          **parser_remaining_kwargs)
            elif('-out' in arg_name):
                if('out_sckpts_file' in arg_name): continue # Pre-unf utility only
                output_files_args.add_argument(arg_name, help=arg_help, metavar='FILE',
                                               **parser_remaining_kwargs)
            elif('file' in arg_name):
                input_files_args.add_argument(arg_name, help=arg_help, metavar='FILE',
                                              **parser_remaining_kwargs)
            elif(('symm' in arg_name) or ('skip_propose_pc' in arg_name)):
                symm_args.add_argument(arg_name, help=arg_help,
                                          **parser_remaining_kwargs)
            elif(arg_name in ['-spin_channel']):
                spin_args.add_argument(arg_name, help=arg_help, metavar='ISPIN',
                                          **parser_remaining_kwargs)
            elif(arg_name in ['-continue_if_not_commensurate', 
                              '--continue_if_npw_smaller_than_expected',
                              '-dont_unfold']):
                debug_args.add_argument(arg_name, help=arg_help,
                                        **parser_remaining_kwargs)
            else:
                self.add_argument(arg_name, help=arg_help, **parser_remaining_kwargs)
        # Energy-related stuff
        egrid_args.add_argument('-efermi', type=float, default=None,
                                help='Fermi energy (eV). Take from self-consist. calc.!')
        egrid_args.add_argument('-emin', type=float, default=None, 
                                help='Minimum energy in the grid')
        egrid_args.add_argument('-emax', type=float, default=None, 
                                help='Maximum energy in the grid')
        egrid_args.add_argument('-dE', type=float, default=None, metavar='dE',
                                help='Energy grid spacing.')

    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one
        args, unknown_args = (
            super(BandUpArgumentParser, self).parse_known_args())
        args.argv = self.get_argv(args)
        return args, unknown_args
    def get_argv(self, args, run_dir=None):
        # Working out which of the supported BandUP options have actually been passed, 
        # and generating an arg list to pass to Popen
        argv = []
        for cla_dict in self.bandup_registered_clas:
            arg_name = '%s'%(cla_dict['key'])
            if(arg_name in sys.argv):
                argv.append(arg_name)
                if(cla_dict['kkind'] in
                   ['cla_int', 'cla_float', 'cla_char', 'cla_xchar', 'cla_logical']):
                    argparse_var_name = arg_name.strip('-').strip('-')
                    val = getattr(args, argparse_var_name) 
                    # Relative paths instead of absolute ones, to get shorter paths
                    if(('file' in argparse_var_name) and (run_dir is not None)):
                        abspath = os.path.abspath(val)
                        path_rel_to_run_dir = os.path.relpath(abspath, run_dir)
                        val = path_rel_to_run_dir
                    argv.append(val)
        return argv

class BandUpPlotArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        remove_conflicts_with_bandup = False
        if('remove_conflicts_with_bandup' in kwargs):
            remove_conflicts_with_bandup = kwargs['remove_conflicts_with_bandup']
            del kwargs['remove_conflicts_with_bandup']
        if(remove_conflicts_with_bandup):
            kwargs['add_help'] = False
        super(BandUpPlotArgumentParser, self).__init__(*args, **kwargs)

        self.aux_settings = {}
        self.add_argument('input_file', default='unfolded_EBS_symmetry-averaged.dat', 
                          nargs='?', help='Name of the input file.')
        self.add_argument('output_file', nargs='?', default=None, 
                           help='%s %s %s'%('Optional: Name of the output file.', 
                                            'If not given, it will be based on',
                                            'the name of the input file.'))
        self.add_argument('-kpts', '--kpoints_file', default='KPOINTS_prim_cell.in', 
                          help=('Name of the file containing information' + 
                                'about the primitive cell k-points. '+
                                 'Default: KPOINTS_prim_cell.in'))
        if(not remove_conflicts_with_bandup):
            self.add_argument('-pc_file', '--prim_cell_file', 
                              default='prim_cell_lattice.in', 
                              help=('Name of the file containing information'
                                    'about the primitive cell lattice vectors. ' +
                                    'Default: prim_cell_lattice.in'))
        self.add_argument('-efile', '--energy_info_file', default='energy_info.in', 
                          help=('Name of the file containing information about ' +
                                'the energy grid and Fermi energy to be used. ' +
                                'Default: energy_info.in. '
                                'This file is optional for the plotting tool.'))
        # Colormaps
        self.aux_settings['cmap_names'],self.aux_settings['cmaps']=get_available_cmaps() 
        # Change self.aux_settings['default_cmap'] if you want 
        # another colormap to be the default one. 
        # I used 'gist_ncar' in my 2014 PRB(R) paper. 
        self.aux_settings['default_cmap'] = 'gist_ncar' 
        if('gist_ncar' not in self.aux_settings['cmap_names']):
            self.aux_settings['default_cmap'] = self.aux_settings['cmap_names'][0]
        self.possible_cmap_choices = self.add_mutually_exclusive_group()
        self.possible_cmap_choices.add_argument('-cmap', '--colormap', default=None, 
                                                choices=self.aux_settings['cmap_names'], 
                        help='Choice of colormap for the plots (name of the colormap).'
                             'You might want to try a few.', metavar='colormap_name')
        self.possible_cmap_choices.add_argument('-icmap', '--icolormap', type=int, 
                             default=None, 
                             choices=range(len(self.aux_settings['cmap_names'])), 
                             help='Choice of colormap for the plots (integer number).', 
                             metavar='0~'+str(len(self.aux_settings['cmap_names']) - 1))
        # E-Fermi and high-symm lines
        mpl_colors = get_matplotlib_color_names()
        self.add_argument('--high_symm_linecolor', default=None, choices=mpl_colors, 
           help='Color of the lines marking the positions of the high-symmetry points.', 
                 metavar='some_matplotlib_color_name')
        self.add_argument('--e_fermi_linecolor', default=None, choices=mpl_colors, 
                          help='Color of the Fermi energy line.', 
                          metavar='some_matplotlib_color_name')
        self.add_argument('--line_width_E_f', type=float, default=None)
        self.add_argument('--line_width_high_symm_points', type=float, default=None)
        linestyles = ['solid', 'dashed', 'dashdot', 'dotted', '-', '--', '-.', ':']
        self.add_argument('--line_style_high_symm_points', default=None, 
                          choices=linestyles)
        self.add_argument('--line_style_E_f', default=None, choices=linestyles)
        self.add_argument('--tick_marks_size', type=float, default=8.0)

        self.add_argument('--save', action='store_true', default=False, 
                          help='Saves the figue to a file. Default: False')
        self.add_argument('--show', action='store_true', default=False, 
           help='Shows the figue. Default: False if --save is selected, True otherwise.')
        self.add_argument('--saveshow', action='store_true', default=False, 
                          help='Saves the figue to a file and opens the saved file'
                               'instead of creating a pyplot window. Default: False')
        self.add_argument('-res', '--fig_resolution', default='m', choices=['l','m','h'],
                          help='Resolution of the figure:' 
                               'l = 100 dpi, m = 300 dpi, h = 600 dpi.'
                               'Default: m = 300 dpi')
        self.add_argument('--disable_auto_round_vmin_and_vmax', action='store_true', 
                          help='Disable normalization of the vmin and vmax of '
                               'the color scale to integer numbers.')

        # File format of the output figure
        # Change the following line if you want another format to be the default.
        self.aux_settings['default_fig_format'] = set_default_fig_format('tiff') 
        self.add_argument('-fmt', '--file_format', 
                          default=self.aux_settings['default_fig_format'],
                          help='File format of the figure. Default: ' +
                               self.aux_settings['default_fig_format'])

        self.possible_fig_orientations = self.add_mutually_exclusive_group()
        self.possible_fig_orientations.add_argument('--landscape', action='store_true', 
                                                    default=False)
        self.possible_fig_orientations.add_argument('--portrait', action='store_true', 
                                                    default=False)

        self.add_argument('-ar', '--aspect_ratio', type=type(Fraction('3/4')), 
                          default=Fraction('3/4'), 
                          help='Aspect ratio of the generated plot. Default: 3/4')
        self.add_argument('--no_ef', action='store_true',help='Hides the E-Fermi line.')

        self.add_argument('--no_symm_lines', action='store_true', 
                          help='Hides the high-symmetry lines.')
        self.add_argument('--no_symm_labels', action='store_true',
                          help='Hides the high-symmetry labels.')
        self.add_argument('--no_cb', action='store_true',help='Hides the colorbar.')


        self.plot_spin_projections = self.add_mutually_exclusive_group()
        self.plot_spin_projections.add_argument('--plot_sigma_x', action='store_true', 
                                                default=False)
        self.plot_spin_projections.add_argument('--plot_sigma_y', action='store_true', 
                                                default=False)
        self.plot_spin_projections.add_argument('--plot_sigma_z', action='store_true', 
                                                default=False)
        self.plot_spin_projections.add_argument('--plot_spin_perp', action='store_true',
                                                default=False)
        self.plot_spin_projections.add_argument('--plot_spin_para', action='store_true',
                                                 default=False)
        self.add_argument('--clip_spin', type=float, default=None)
        self.add_argument('--spin_marker', default='_', choices=['_', 'o'])


        self.add_argument('-nlev', '--n_levels', type=int, default=101, 
                          help='Number of different levels used in the contour plot.'
                               'Default: 101')
        self.add_argument('-kmin', type=float, default=None)
        self.add_argument('-kmax', type=float, default=None)
        self.add_argument('-shift_k', '--shift_kpts_coords', type=float, default=0.0, 
                          help='Shift in the k-points. Default: 0.0')

        self.add_argument('-emin', type=float, default=None)
        self.add_argument('-emax', type=float, default=None)
        self.add_argument('-dE', type=float, default=None)
        self.add_argument('-shift_e', '--shift_energy', type=float, default=0.0,
                          help='Shift in the energy grid. Default: 0.0')

        self.add_argument('-interp', '--interpolation', default=None, 
                          choices=['nearest', 'linear', 'cubic'], 
                          help='Interpolation scheme used. Default: nearest')

        self.add_argument('-vmin', '--minval_for_colorbar', type=float, 
                          help='Value to which the first color of the colormap'
                               'will be normalized.')
        self.add_argument('-vmax', '--maxval_for_colorbar', type=float, 
                          help='Value to which the last color of the colormap '
                               'will be normalized.')
        self.add_argument('--round_cb', type=int, default=1, 
                          help='Number of decimal digits displayed'
                               'in the colobar ticks.')

        self.show_cb_label_options = self.add_mutually_exclusive_group()
        self.show_cb_label_options.add_argument('--cb_label', action='store_true', 
                                                help='Show the colorbar label.')
        self.show_cb_label_options.add_argument('--cb_label_full', action='store_true', 
                                                help='Show the colorbar label (full).')

        
        self.add_argument('--skip_grid_freq_clip', action='store_true', 
                          help=argparse.SUPPRESS)
        # Flag to know if this script has been called by the GUI
        self.add_argument('--running_from_GUI', action='store_true', 
                          help=argparse.SUPPRESS)


    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one
        args, unknown_args = (
            super(BandUpPlotArgumentParser, self).parse_known_args(*fargs, **fkwargs))
        args = self.filter_args_plot(args)
        return args, unknown_args
    def filter_args_plot(self, args):
        args.aux_settings = self.aux_settings
        if(args.icolormap is not None):
            args.colormap = self.aux_settings['cmap_names'][args.icolormap]
        self.using_default_cmap = args.colormap is None
        if(self.using_default_cmap):
            args.colormap = self.aux_settings['default_cmap']

        if(args.saveshow):
            args.save = True

        if(args.output_file is not None):
            args.output_file =(
                abspath(
                    output_file_with_supported_extension(
                        args.output_file, self.aux_settings['default_fig_format']
                                                            )
                       )
            )
            args.save = True

        if(not args.save and not args.show):
            args.save = False
            args.show = True

        args.aspect_ratio = float(args.aspect_ratio)
        if(args.aspect_ratio > 1.0):
            args.aspect_ratio = 1.0 / args.aspect_ratio

        return args

class BandUpPythonArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if('formatter_class' not in kwargs):
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        super(BandUpPythonArgumentParser, self).__init__(*args, **kwargs)

        self.bandup_parser = BandUpArgumentParser(add_help=False)
        self.bandup_plot_parser = BandUpPlotArgumentParser(add_help=False)
        self.bandup_pre_unf_parser = BandUpPreUnfoldingArgumentParser(add_help=False)

        # Defining available tasks and setting the default one to 'unfolding'
        self.allowed_tasks = OrderedDict()
        self.allowed_tasks['pre-unfold'] = {'subparser_name':'bandup',
                                            'help':"Runs BandUP's pre-unfolding tool "+
                                            "to get the SC-KPTs needed for unfolding.",
                                            'parents':[self.bandup_pre_unf_parser],
                                         'parent_class':BandUpPreUnfoldingArgumentParser}
        self.allowed_tasks['unfold'] = {'subparser_name':'bandup',
                                        'help':"Runs BandUP's main code",
                                        'parents':[self.bandup_parser],
                                        'parent_class':BandUpArgumentParser}
        self.allowed_tasks['plot'] = {'subparser_name':'bandup_plot',
                                      'help':"Plots BandUP's output files.",
                                      'parents':[self.bandup_plot_parser],
                                      'parent_class':BandUpPlotArgumentParser}
        self.default_main_task = 'unfold'

        # Implementing each task with a different subparser
        self.subparsers = self.add_subparsers(dest='main_task',
                              help='Task to be performed. (default: unfold)')
        # This is to prevent the .subparsers.add_parser method to fail,
        # attempting to init this class recursively
        self.subparsers._parser_class = argparse.ArgumentParser
        for task, task_info in self.allowed_tasks.iteritems():
            subparser = self.subparsers.add_parser(task, 
                help=task_info['help'], parents=task_info['parents'],
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            setattr(self, '%s_subparser'%(task_info['subparser_name']), subparser)

        # Adding task-dependent options
        # Plot subparser:
        bandup_plot_extra_opts = self.bandup_plot_subparser.add_argument_group(
                                     'Options available'+
                                     ' only through this Python interface')
        bandup_plot_extra_opts.add_argument('-plotdir', 
            default=os.path.join(working_dir,"plots"),
            help='Directory where plot will be saved')
        bandup_plot_extra_opts.add_argument('-results_dir', 
                                            default=defaults['results_dir'],
            help='Dir where BandUP was run.')
        bandup_plot_extra_opts.add_argument('--overwrite', action='store_true',
            help='Overwrite files/directories if copy paths coincide.')
        # BandUP subparser
        bandup_extra_opts = self.bandup_subparser.add_argument_group('Options available'+
                                                 ' only through this Python interface')
        bandup_extra_opts.add_argument('-wavefunc_calc_dir', 
                                       default=defaults['wavefunc_calc_dir'],
            help='Directory where the WFs to be unfolded are stored.')
        bandup_extra_opts.add_argument('-self_consist_calc_dir', 
                                       default=defaults['self_consist_calc_dir'],
            help='Dir containing the self-consistent calc files.')
        bandup_extra_opts.add_argument('-results_dir', default=defaults['results_dir'],
            help='Dir where the BandUP will be run.')
        bandup_extra_opts.add_argument('--overwrite', action='store_true',
            help='Overwrite files/directories if copy paths coincide.')


    def print_help(self, *args, **kwargs):
        calling_script = sys.argv[0]
        print self.format_help()
        extra_help_msg = "Each task has its own help as well, which can be requested\n"
        extra_help_msg+= 'by passing "-h" or "--help" after the name of the task.\n'
        extra_help_msg+= 'Eg.: %s unfold -h'%(calling_script)
        print extra_help_msg
        print ''
        sys.exit(0)

    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one

        # Finding 1st positional arg passed to program. This will be used to handle the 
        # "no subparser selected" and set task to self.default_main_task in this case
        first_positional_arg = None
        i_first_pos_arg = None
        for iarg, arg in enumerate(sys.argv[1:]):
            if(not arg.startswith('-')): 
                first_positional_arg = arg
                i_first_pos_arg = iarg + 1
                break

        # To allow abbreviations in the task names:
        if(first_positional_arg is not None):
            n_compat_choices = 0
            compatible_task = None
            for allowed_task in self.allowed_tasks:
                if(first_positional_arg==allowed_task[:len(first_positional_arg)]):
                    n_compat_choices += 1
                    compatible_task = allowed_task
            if(n_compat_choices!=1): compatible_task = None
            if(compatible_task is not None): 
                first_positional_arg = compatible_task
                sys.argv[i_first_pos_arg] = first_positional_arg

        # Defining 'unfold' as default task
        task_defined = first_positional_arg in self.allowed_tasks
        help_requested = len(set(['-h', '--help']).intersection(sys.argv[1:])) > 0
        if(not task_defined and not help_requested):
            sys.argv.insert(1, self.default_main_task) 

        # Making sure arguments are compatible
        needed_energy_opts = ['-emin', '-emax', '-dE']
        # EFermi will be read from file if not passed
        individual_energy_opts = needed_energy_opts + ['-efermi']

        efile_passed = (
            len(set(['-efile', '--energy_info_file']).intersection(sys.argv[1:]))>0
        )
        individual_energy_opts_passed = (
            len(set(individual_energy_opts).intersection(sys.argv[1:]))>0
        )
        if(efile_passed and individual_energy_opts_passed):
            msg = ('%s: error: Args "-efile" or "--energy_info_file"'%(sys.argv[0]) +
                   ' not allowed with any of:\n'+
                   '"-efermi", "-emin", "-emax", "-dE".')
            self.error(msg)
        elif(individual_energy_opts_passed):
            all_needed_ener_opts_passed = (
                len(set(needed_energy_opts).intersection(sys.argv[1:])) == 
                len(set(needed_energy_opts))
            )
            if(not all_needed_ener_opts_passed):
                msg = 'All options from %s'%(', '.join(needed_energy_opts))
                msg += ' need to be passed if one of them is passed. Alternatively, '
                msg += 'you can also specify these values in a file using the option '
                msg += '"--energy_info_file FILENAME".'
                self.error(msg)


        # Now finally parsing args
        args, unknown_args = (
            super(BandUpPythonArgumentParser, self).parse_known_args(*fargs, **fkwargs))
        args = self.filter_args(args) 

        return args, unknown_args

    def filter_args(self, args, *fargs, **fkwargs):
        if(args.main_task == 'unfold'):
            subparser = self.bandup_parser
            # args.argv will be passed to Popen to run BandUP's Fortran core code
            args.argv = subparser.get_argv(args, run_dir=args.results_dir)
            if(args.castep):
                if(not args.seed):
                    # Working out the seed
                    bands_file = [os.path.join(args.self_consist_calc_dir,fname) for 
                                  fname in
                                  os.listdir(args.self_consist_calc_dir) if
                                  fname.endswith('.bands')][0]
                    args.seed = os.path.splitext(os.path.basename(bands_file))[0]
                    print 'WARNING: Seed not passed as argument. Using "%s".'%(
                           args.seed)
                    args.argv.append('-seed')
                    args.argv.append(args.seed)
        elif(args.main_task == 'plot'):
            subparser = self.bandup_plot_parser
            args = subparser.filter_args_plot(args)
        elif(args.main_task == 'pre-unfold'):
            subparser = self.bandup_pre_unf_parser

        try:
            if(args.efermi is None): 
                args.efermi = get_efermi(args)
                if(args.efermi is None):
                    warnings.warn('Could not get E-Fermi!') 
        except(AttributeError):
            pass

        args.default_values = {}
        for arg in dir(args):
            if(arg.startswith('_')): continue
            try:
                args.default_values[arg] = subparser.get_default(arg)
            except(TypeError):
                pass
        return args
