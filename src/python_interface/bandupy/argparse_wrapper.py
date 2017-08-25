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
from six import iteritems
# Imports from within the package
from .version import __version__
from .constants import (
    BANDUP_DIR,
    WORKING_DIR,
    BANDUP_BIN,
    BANDUP_PRE_UNFOLDING_BIN,
)
from .defaults import defaults
from .files import (
    continuation_lines, 
    get_efermi, 
    valid_path, 
    assert_exec_exists,
    guess_castep_seed,
)
from .figs import (
    get_matplotlib_color_names,
    get_available_cmaps,
    set_default_fig_format,
)
from .warnings_wrapper import warnings
from .sysargv import arg_passed
from .lists import str_list_to_int_range
from .build import castep_interface_available

def store_abs_path_action_gen(assert_existence=False, rel_path_start=WORKING_DIR):
    class StoreAbsPath(argparse.Action):
        """ Action to store an absolute path derived from the passed relative path. 

            It is useful to remember that
            argparse does not use the action when applying the default

        """
        def __init__(self, option_strings, dest, **kwargs):
            super(StoreAbsPath, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # If the path argument is explicitly passed, then it is assumed
            # to be relative to WORKING_DIR
            new_path = values
            if(values is not None):
                new_path = os.path.abspath(os.path.relpath(values, rel_path_start))
                if(assert_existence):
                    if(not os.path.exists(new_path)):
                        msg = '"%s" argument:\n'%(self.dest)
                        msg += 'Path\n\n'
                        msg += 4*' ' + "%s\n\n"%(new_path)
                        msg += 'does not exist or is not accessible.'
                        parser.error(msg)
            setattr(namespace, self.dest, new_path)
    return StoreAbsPath

class AssertCastepInterfaceAvailable(argparse.Action):
    """ Produce an error if using "-castep" without compiling BandUP with CASTEP support

    """
    def __init__(self, option_strings, dest, default=False, required=False, help=None):
        super(AssertCastepInterfaceAvailable, self).__init__(
            option_strings, dest, nargs=0, default=default, required=required, help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        if(castep_interface_available()):
            setattr(namespace, self.dest, True)
        else:
            msg = 'BandUP has been compiled without CASTEP support!'
            parser.error(msg)

def obsolete_arg_action_gen(alternative_option=None):
    class RefuseObsoleteArgs(argparse.Action):
        """ This action will produce an error is an obsolete option has been passed 

            This assumes that the default value of such options is argparse.SUPPRESS

        """
        def __init__(self, option_strings, dest, **kwargs):
            super(RefuseObsoleteArgs, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            if(values != argparse.SUPPRESS):
                msg = 'The option "%s" is no longer accepted.'%(self.dest)
                if(alternative_option is not None):
                    msg += '\n' + 22*' '
                    msg += 'Please use "%s %s" instead.'%(alternative_option, values)
                parser.error(msg)
    return RefuseObsoleteArgs

def str_list_to_int_range_arg_action_gen(convert_fortran2python_indexing,
                                 include_last_from_subranges
                                 ):
    class StrList2IntRange(argparse.Action):
        """ This action converts a range pased as string to a list of integers 

            This assumes that the string items are given either as "I0" or "I0-I1"

        """
        def __init__(self, option_strings, dest, **kwargs):
            super(StrList2IntRange, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            try:
                int_range = str_list_to_int_range(
                         values, 
                         convert_fortran2python_indexing=convert_fortran2python_indexing,
                         include_last_from_subranges=include_last_from_subranges
                            )
                setattr(namespace, self.dest, int_range)
            except(ValueError):
                msg = 'Invalid range.'
                parser.error(msg)
    return StrList2IntRange 

def get_bandup_registered_clas(bandup_path=BANDUP_DIR):
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
                output_files_args.add_argument(arg_name, 
                    help=arg_help.replace('(pre-unfolding utility)',''), 
                    action=store_abs_path_action_gen(assert_existence=False), 
                    **parser_remaining_kwargs)
            elif(('no_symm' in arg_name) or ('skip_propose_pc' in arg_name)):
                symm_args.add_argument(arg_name, help=arg_help,
                                          **parser_remaining_kwargs)
            elif(arg_name in ['-pc_file', '-sc_file', '-pckpts_file']):
                input_files_args.add_argument(arg_name, help=arg_help,
                    action=store_abs_path_action_gen(assert_existence=True), 
                    **parser_remaining_kwargs)
    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one
        args, unknown_args = (
            super(BandUpPreUnfoldingArgumentParser, self).parse_known_args())
        args.argv = self.get_argv(args)
        return args, unknown_args
    def get_argv(self, args, run_dir=None):
        # Working out which of the supported BandUP options have actually been passed, 
        # and generating an arg list to pass to Popen
        argv = []
        for cla_dict in self.bandup_registered_clas:
            arg_name = '%s'%(cla_dict['key'])
            if(arg_passed(arg_name)):
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

class BandUpAtOrbProjUnfoldingArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if('formatter_class') not in kwargs:
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        if('add_help' not in kwargs): kwargs['add_help'] = False
        super(BandUpAtOrbProjUnfoldingArgumentParser, self).__init__(*args, **kwargs)
        self.add_argument('-spin_channel', default=1, type=int, choices=[1,2])
        self.add_argument('-orbs', default=None, nargs='+',
                          help=(
                                'The orbitals to be included in the sum. If not '+
                                'passed, then all available orbitals will be used. '+
                                'Ex.: -orbs px'
                               )
                         )
        self.add_argument('-atom_indices', default=None, nargs='+',
                          action=str_list_to_int_range_arg_action_gen(
                                     convert_fortran2python_indexing=True,
                                     include_last_from_subranges=True,
                                 ),
                          help=(
                                'The indices of the atoms to be included, *starting '+
                                'from 1*. If not passed, then all atoms will be used. '+
                                'Ranges are accepted. Ex: -at 1 3 5-18 22-30 41 52'
                               )
                         )
        self.add_argument('-results_dir', 
            action=store_abs_path_action_gen(assert_existence=False), 
            default=defaults['results_dir'],
            help=('Dir where the "orbital_projections" dir and the unfolding-density '+
                  ' operator files are located. This is also the directory where the '+
                  'output file will be saved.'
                 )
        )

class BandUpArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
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
                if(not castep_interface_available()):
                    arg_help = argparse.SUPPRESS
                parser_remaining_kwargs['action']=AssertCastepInterfaceAvailable
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
                    action=store_abs_path_action_gen(assert_existence=True), 
                    **parser_remaining_kwargs)
            elif('-out' in arg_name or '_out_' in arg_name):
                if('out_sckpts_file' in arg_name): continue # Pre-unf utility only
                output_files_args.add_argument(arg_name, help=arg_help, metavar='FILE',
                    action=store_abs_path_action_gen(assert_existence=False), 
                    **parser_remaining_kwargs)
            elif('-write' in arg_name):
                output_files_args.add_argument(arg_name, help=arg_help, 
                    **parser_remaining_kwargs)
            elif('file' in arg_name):
                input_files_args.add_argument(arg_name, help=arg_help, metavar='FILE',
                    action=store_abs_path_action_gen(assert_existence=True), 
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
                                help='Minimum energy in the grid (eV).')
        egrid_args.add_argument('-emax', type=float, default=None, 
                                help='Maximum energy in the grid (eV).')
        egrid_args.add_argument('-dE', type=float, default=0.05, metavar='dE',
                                help='Energy grid spacing (eV).')

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
            if(arg_passed(arg_name)):
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
        super(BandUpPlotArgumentParser, self).__init__(*args, **kwargs)

        # Non-argument definitions
        self.aux_settings = {}
        mpl_colors = get_matplotlib_color_names()
        linestyles = ['solid', 'dashed', 'dashdot', 'dotted', '-', '--', '-.', ':']
        # File format of the output figure
        # Change the following line if you want another format to be the default.
        self.aux_settings['default_fig_format'] = set_default_fig_format('tiff') 
        # Colormaps
        self.aux_settings['cmap_names'],self.aux_settings['cmaps']=get_available_cmaps() 
        # Change self.aux_settings['default_cmap'] if you want 
        # another colormap to be the default one. 
        # I used 'gist_ncar' and 'Greys' in my PRB(R) 2014 and 2015 papers, respectively 
        self.aux_settings['default_cmap'] = 'jet' 
        if('jet' not in self.aux_settings['cmap_names']):
            self.aux_settings['default_cmap'] = self.aux_settings['cmap_names'][0]

        # Argument groups
        obsolete_args = self.add_argument_group('Obsolete/discontinued options')
        kptgrid_args = self.add_argument_group('Options controlling the kpt grid')
        efermi_args = self.add_argument_group('Fermi energy line options')
        graph_args = self.add_argument_group('Graph appearance')
        spinor_args = self.add_argument_group('Spinor-related options')
        hsline_args = self.add_argument_group('High-symmetry line options')
        egrid_args = self.add_argument_group('Options controlling the energy grid')
        cbar_args = self.add_argument_group('Colorbar options')
        colormap_args = self.add_argument_group('Options controlling colormap settings')
        output_files_args = self.add_argument_group('Output file options')
        input_files_args = self.add_argument_group('Input file options')

        # Obsolete options
        # These are no longer accepted, but were kept here to warn those that 
        # normally use them and, when applicable, offer alternatives
        #obsolete_args.add_argument('positional_input_file',
        #    default=argparse.SUPPRESS, 
        #    action=obsolete_arg_action_gen('-input_file'), 
        #    nargs='?', help=argparse.SUPPRESS)
        #obsolete_args.add_argument('positional_output_file',
        #    default=argparse.SUPPRESS, 
        #    action=obsolete_arg_action_gen('-output_file'), 
        #    nargs='?', help=argparse.SUPPRESS)

        # Input files
        input_files_args.add_argument('-input_file', 
            default='unfolded_EBS_symmetry-averaged.dat',
            action=store_abs_path_action_gen(assert_existence=True), 
            help='Name of the input file.')
        input_files_args.add_argument('-kpts', '--kpoints_file', '-pckpts_file',
            default='KPOINTS_prim_cell.in', metavar='FILE',
            action=store_abs_path_action_gen(assert_existence=True), 
            help=('Name of the file containing information' + 
                  'about the primitive cell k-points. '
                 )
        )
        input_files_args.add_argument('-pc_file', '--prim_cell_file', metavar='FILE', 
            default='prim_cell_lattice.in', 
            action=store_abs_path_action_gen(assert_existence=True), 
            help=('Name of the file containing information'
                  'about the primitive cell lattice vectors. '
                 )
        )
        input_files_args.add_argument('-efile', '--energy_info_file', '-energy_file',
            default='energy_info.in', metavar='FILE', 
            action=store_abs_path_action_gen(assert_existence=True), 
            help=('Name of the file containing information about '+
                  'the energy grid and Fermi energy to be used. '+
                  'This file is optional for the plotting tool.'))

        # Deciding how to output results
        output_files_args.add_argument('-output_file', default=None, 
            action=store_abs_path_action_gen(assert_existence=False), 
            help='%s %s'%('Optional: Name of the output file.', 
                 'If not given, it will be based on the name of the input file.'))
        output_files_args.add_argument('--save', action='store_true', 
                                       default=False, 
                                       help='Saves the figue to a file')
        output_files_args.add_argument('--show', action='store_true', default=False, 
           help='Shows the figue. Default: False if --save is selected, True otherwise.')
        output_files_args.add_argument(
            '--saveshow', action='store_true', default=False, 
            help=('Saves the figue to a file and opens the saved file'
                  'instead of creating a pyplot window')
        )
        output_files_args.add_argument('-res', '--fig_resolution', default='m',
                                       choices=['l','m','h'],
                                       help='Resolution of the figure:' 
                                       'l = 100 dpi, m = 300 dpi, h = 600 dpi.'
                                      )
        output_files_args.add_argument(
            '-fmt', '--file_format', 
            default=self.aux_settings['default_fig_format'],
            help='File format of the figure.'
        )

        # Colormap stuff
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

        # High-symm lines
        hsline_args.add_argument('--high_symm_linecolor', default=None, 
                                 choices=mpl_colors, 
           help='Color of the lines marking the positions of the high-symmetry points.', 
           metavar='some_matplotlib_color_name')
        hsline_args.add_argument('--line_width_high_symm_points', type=float, 
                                 default=None)
        hsline_args.add_argument('--line_style_high_symm_points', default=None, 
                                 choices=linestyles)
        hsline_args.add_argument('--no_symm_lines', action='store_true', 
                                 help='Hides the high-symmetry lines.')
        hsline_args.add_argument('--no_symm_labels', action='store_true',
                                 help='Hides the high-symmetry labels.')
        # E-Fermi line
        efermi_args.add_argument('--e_fermi_linecolor', default=None, 
                                 choices=mpl_colors, 
                                 help='Color of the Fermi energy line.', 
                                 metavar='SOME_MATPLOTLIB_COLOR_NAME')
        efermi_args.add_argument('--line_width_E_f', type=float, default=None,
                                  metavar='WIDTH')
        efermi_args.add_argument('--line_style_E_f', default=None, choices=linestyles)
        efermi_args.add_argument('--no_ef', action='store_true',
                                 help='Hides the E-Fermi line.')

        # Appearance of figure
        graph_args.add_argument('--tick_marks_size', type=float, default=8.0,
                                metavar='SIZE')

        graph_args.add_argument('--disable_auto_round_vmin_and_vmax', 
                                action='store_true', 
                                help='Disable normalization of the vmin and vmax of '
                                     'the color scale to integer numbers.')
        graph_args.add_argument('-ar', '--aspect_ratio', 
                                type=type(Fraction('3/4')), 
                                default=Fraction('3/4'), 
                                help='Aspect ratio of the generated plot.')
        graph_args.add_argument('-interp', '--interpolation', 
                                default="nearest", 
                                choices=['nearest', 'linear', 'cubic'], 
                                help='Interpolation scheme used.')
        graph_args.add_argument(
            '-nlev', '--n_levels', type=int, default=101, 
            help='Number of different levels used in the contour plot.'
        )
        self.possible_fig_orientations = self.add_mutually_exclusive_group()
        self.possible_fig_orientations.add_argument('--landscape', action='store_true', 
                                                    default=False)
        self.possible_fig_orientations.add_argument('--portrait', action='store_true', 
                                                    default=False)

        # Spinor/spin stuff
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
        spinor_args.add_argument('--clip_spin', type=float, default=None,
                                 metavar='VALUE')
        spinor_args.add_argument('--spin_marker', default='_', choices=['_', 'o'])

        # Kpt grid
        kptgrid_args.add_argument('-kmin', type=float, default=None)
        kptgrid_args.add_argument('-kmax', type=float, default=None)
        kptgrid_args.add_argument('-shift_k', '--shift_kpts_coords', 
                                  type=float, default=0.0, metavar='SHIFT', 
                                  help='Shift in the k-points.')
        # Energy grid
        egrid_args.add_argument('-emin', type=float, default=None)
        egrid_args.add_argument('-emax', type=float, default=None)
        egrid_args.add_argument('-dE', type=float, default=None, metavar='dE_VAL')
        egrid_args.add_argument(
            '-shift_e', '--shift_energy', type=float, default=0.0,
            metavar='SHIFT',
            help='Shift in the energy grid.'
        )

        # Colorbar
        cbar_args.add_argument('--no_cb', action='store_true',help='Hides the colorbar.')
        cbar_args.add_argument('-vmin', '--minval_for_colorbar', type=float, 
                               metavar='VAL',
                               help='Value to which the first color of the colormap'
                               'will be normalized.')
        cbar_args.add_argument('-vmax', '--maxval_for_colorbar', type=float,
                               metavar='VAL',
                               help='Value to which the last color of the colormap '
                               'will be normalized.')
        cbar_args.add_argument('--round_cb', type=int, default=1, 
                               metavar='NDEC',
                               help='Number of decimal digits displayed'
                               'in the colobar ticks.')
        self.show_cb_label_options = self.add_mutually_exclusive_group()
        self.show_cb_label_options.add_argument('--cb_label', action='store_true', 
                                                help='Show the colorbar label.')
        self.show_cb_label_options.add_argument('--cb_label_full', action='store_true', 
                                                help='Show the colorbar label (full).')
        
        # Args not intended to be passed directly by users
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
        args.aux_settings['using_default_cmap'] = args.colormap is None
        if(args.aux_settings['using_default_cmap']):
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
        self._positionals.title = 'BandUP tasks (positional arguments)'

        self.bandup_parser = BandUpArgumentParser(add_help=False)
        self.bandup_plot_parser = BandUpPlotArgumentParser(add_help=False)
        self.bandup_pre_unf_parser = BandUpPreUnfoldingArgumentParser(add_help=False)
        self.bandup_at_proj_parser = (
            BandUpAtOrbProjUnfoldingArgumentParser(add_help=False)
        )

        # Print code version
        self.add_argument('-v', '--version', action='version',
                          version='%(prog)s {version}'.format(version=__version__()))
        # Defining available tasks and setting the default one to 'unfolding'
        self.allowed_tasks = OrderedDict()
        self.allowed_tasks['kpts-sc-get'] = {'subparser_name':'bandup_pre_unf',
                                            'help':"Runs BandUP's pre-unfolding tool "+
                                            "to get the SC-KPTs needed for unfolding.",
                                            'parents':[self.bandup_pre_unf_parser],
                                         'parent_class':BandUpPreUnfoldingArgumentParser}
        self.allowed_tasks['unfold'] = {'subparser_name':'bandup',
                                        'help':"Runs BandUP's main unfolding code.",
                                        'parents':[self.bandup_parser],
                                        'parent_class':BandUpArgumentParser}
        self.allowed_tasks['projected-unfold'] = {'subparser_name':'bandup_at_proj',
                                        'help':(
                                                "Reads SC atomic orbital projections, "+
                                                "creates SC projection operators for "+
                                                "the chosen atoms and/or atomic "+
                                                "orbitals, and then uses unfolding-"+
                                                "density operators to produce orbital- "+
                                                'and/or atom-decomposed unfolded band '+
                                                'structures. N.B.: For this to work, '+
                                                'please run BandUP first using the '+
                                                '"--orbitals" flag.'
                                               ),
                                   'parents':[self.bandup_at_proj_parser],
                                   'parent_class':BandUpAtOrbProjUnfoldingArgumentParser}
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
        for task, task_info in iteritems(self.allowed_tasks):
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
            action=store_abs_path_action_gen(assert_existence=False), 
            default=defaults['plot_dir'],
            help='Directory where plot will be saved')
        bandup_plot_extra_opts.add_argument('-results_dir', 
            action=store_abs_path_action_gen(assert_existence=True), 
            default=defaults['results_dir'],
            help='Dir where BandUP was run.')
        bandup_plot_extra_opts.add_argument('--gui', action='store_true',
            help='Open the GUI.')
        bandup_plot_extra_opts.add_argument('--overwrite', action='store_true',
            help='Overwrite files/directories if copy paths coincide.')
        # BandUP subparser
        bandup_extra_opts = self.bandup_subparser.add_argument_group('Options available'+
                                                 ' only through this Python interface')
        bandup_extra_opts.add_argument('-wavefunc_calc_dir', 
            action=store_abs_path_action_gen(assert_existence=True), 
            default=defaults['wavefunc_calc_dir'],
            help='Directory where the WFs to be unfolded are stored.')
        bandup_extra_opts.add_argument('-self_consist_calc_dir',
            action=store_abs_path_action_gen(assert_existence=True), 
            default=defaults['self_consist_calc_dir'],
            help='Dir containing the self-consistent calc files.')
        bandup_extra_opts.add_argument('-results_dir', 
            action=store_abs_path_action_gen(assert_existence=False), 
            default=defaults['results_dir'],
            help='Dir where the BandUP will be run.')
        bandup_extra_opts.add_argument('--overwrite', action='store_true',
            help='Overwrite files/directories if copy paths coincide.')
        bandup_extra_opts.add_argument('--orbitals', action='store_true', 
            help=('If set, BandUP will attempt to parse projections onto atomic '+
                  'orbitals. These need to have been previously calculated by your '+
                  'plane-wave code. All projections must be available, and they '+
                  'should contain both real and imaginary parts -- BandUP will need '+
                  'these to calculate their duals. This option will also request '+
                  'unfolding-density operators to be written out.'))
        # BandUP's pre-unfolding tool
        bandup_pre_unf_extra_opts = self.bandup_pre_unf_subparser.add_argument_group(
                                     'Options available'+
                                     ' only through this Python interface')
        bandup_pre_unf_extra_opts.add_argument('-inputs_dir', 
            default=defaults['pre_unfolding_inputs_dir'],
            help=('Dir where the input files for the pre-unfolding "kpts-sc-get" task '+
                  'are located. This is also the directory where the output files '+
                  'will be saved.'
                 )
        )

    def print_help(self, *args, **kwargs):
        calling_script = sys.argv[0]
        print(self.format_help())
        extra_help_msg = 'Task name abbreviations are allowed. For instance,\n'
        extra_help_msg += '"./bandup kpts" requests the same task as '
        extra_help_msg += '"./bandup kpts-sc-get".\n'
        extra_help_msg += '\n'
        extra_help_msg += "Each task has its own help as well, which can be requested\n"
        extra_help_msg+= 'by passing "-h" or "--help" after the name of the task.\n'
        extra_help_msg+= 'Eg.: %s unfold -h'%(calling_script)
        print(extra_help_msg)
        print('')
        sys.exit(0)

    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one

        # Positional args passed to program. This will be used to handle the 
        # "no subparser selected" and set task to self.default_main_task in this case
        positional_args = []
        first_positional_arg = None
        for iarg, arg in enumerate(sys.argv[1:]):
            if(arg.startswith('-')): break 
            positional_args.append(arg)
        if(positional_args): 
            first_positional_arg = positional_args[0]

        task_defined = first_positional_arg is not None
        help_requested = len(set(['-h', '--help']).intersection(sys.argv[1:])) > 0
        version_requested = len(set(['-v', '--version']).intersection(sys.argv[1:])) > 0
        if(task_defined):
            # To allow abbreviations in the task names:
            n_compat_choices = 0
            compatible_task = None
            for allowed_task in self.allowed_tasks:
                if(first_positional_arg==allowed_task[:len(first_positional_arg)]):
                    n_compat_choices += 1
                    compatible_task = allowed_task
            if(n_compat_choices!=1): compatible_task = None
            if(compatible_task is not None): 
                i_first_pos_arg = sys.argv.index(first_positional_arg)
                first_positional_arg = compatible_task
                sys.argv[i_first_pos_arg] = first_positional_arg
        elif(not help_requested and not version_requested):
            # Defining 'unfold' as default task
            sys.argv.insert(1, self.default_main_task) 

        # Making sure arguments are compatible
        needed_energy_opts = ['-emin', '-emax']
        # EFermi will be read from file if not passed
        individual_energy_opts = needed_energy_opts + ['-dE', '-efermi']

        efile_passed = arg_passed('-efile') or arg_passed('--energy_info_file')
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
            super(BandUpPythonArgumentParser, self).parse_known_args(
                                                        *fargs, **fkwargs))
        args = self.filter_args(args) 

        return args, unknown_args

    def filter_args(self, args, *fargs, **fkwargs):
        if(args.main_task == 'kpts-sc-get'):
            assert_exec_exists(BANDUP_PRE_UNFOLDING_BIN)
            subparser = self.bandup_pre_unf_parser
            # args.argv will be passed to Popen to run BandUP's Fortran core code
            args.argv = subparser.get_argv(args, run_dir=args.inputs_dir)
        elif(args.main_task == 'unfold'):
            assert_exec_exists(BANDUP_BIN)
            subparser = self.bandup_parser
            # args.argv will be passed to Popen to run BandUP's Fortran core code
            args.argv = subparser.get_argv(args, run_dir=args.results_dir)
            if(args.castep and not args.seed):
                args.seed = guess_castep_seed(args)
                args.argv.append('-seed')
                args.argv.append(args.seed)
                warnings.warn('Seed not passed as argument. Using "%s".'%(args.seed))
            # Asking BandUP to write unfolding-density operators if --orbitals is passed
            if(arg_passed('--orbitals')):
                args.argv.append('-write_unf_dens_op')
        elif(args.main_task == 'projected-unfold'):
            subparser = self.bandup_at_proj_parser
        elif(args.main_task == 'plot'):
            subparser = self.bandup_plot_parser
            args = subparser.filter_args_plot(args)

        if(hasattr(args, 'efermi')):
            # Using hasattr here generated cleaner code than another try->except
            try:
                if(args.efermi is None): 
                    args.efermi = get_efermi(args)
                    if(args.efermi is None):
                        warnings.warn('Could not get E-Fermi!') 
            except(AttributeError):
                if(args.efermi is None): 
                    warnings.warn('Could not get E-Fermi!') 
                else:
                    raise

        args.default_values = {}
        for arg in dir(args):
            if(arg.startswith('_')): continue
            try:
                args.default_values[arg] = subparser.get_default(arg)
            except(TypeError):
                pass
        return args
