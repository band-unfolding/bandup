import argparse
import os
from bandup_python_wrapper.environ import *
from bandup_python_wrapper.files import continuation_lines
from bandup_python_wrapper.figs import *

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

class BandUpArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if('remove_conflicts_with_bandup' in kwargs): 
            del kwargs['remove_conflicts_with_bandup']
        if('formatter_class') not in kwargs:
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        if('add_help' not in kwargs): kwargs['add_help'] = False
        super(BandUpArgumentParser, self).__init__(*args, **kwargs)

        # Adding BandUP's supported args to the parser
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
                #parser_remaining_kwargs['const'] = arg_name
            self.add_argument(arg_name, help=arg_help, **parser_remaining_kwargs)

    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one
        args, unknown_args = (
            super(BandUpArgumentParser, self).parse_known_args())
        args.argv = self.get_argv(args)
        return args, unknown_args
    def get_argv(self, args):
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

        self.indent = 4 * ' '
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
        self.cmap_names, self.cmaps = get_available_cmaps() 
        # Change self.default_cmap if you want another colormap to be the default one. 
        # I used 'gist_ncar' in my paper. 
        self.default_cmap = 'gist_ncar' 
        if('gist_ncar' not in self.cmap_names):
            self.default_cmap = self.cmap_names[0]
        self.possible_cmap_choices = self.add_mutually_exclusive_group()
        self.possible_cmap_choices.add_argument('-cmap', '--colormap', default=None, 
                                                choices=self.cmap_names, 
                        help='Choice of colormap for the plots (name of the colormap).'
                             'You might want to try a few.', metavar='colormap_name')
        self.possible_cmap_choices.add_argument('-icmap', '--icolormap', type=int, 
                                                default=None, 
                                                choices=range(len(self.cmap_names)), 
                             help='Choice of colormap for the plots (integer number).', 
                                   metavar='0~'+str(len(self.cmap_names) - 1))
        # E-Fermi and high-symm lines
        mpl_colors = sorted(mpl.colors.cnames.keys())
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
        self.default_fig_format = set_default_fig_format('tiff') 
        self.add_argument('-fmt', '--file_format', default=self.default_fig_format,
                          help='File format of the figure. Default: ' +
                               self.default_fig_format)

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
        if(args.icolormap is not None):
            args.colormap = self.cmap_names[args.icolormap]
        self.using_default_cmap = args.colormap is None
        if(self.using_default_cmap):
            args.colormap = self.default_cmap

        if(args.saveshow):
            args.save = True

        if(args.output_file is not None):
            args.output_file =(
                abspath(output_file_with_supported_extension(args.output_file, 
                                                             self.default_fig_format
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
        kwargs['add_help'] = False
        if('formatter_class' not in kwargs):
            kwargs['formatter_class'] = argparse.ArgumentDefaultsHelpFormatter
        super(BandUpPythonArgumentParser, self).__init__(*args, **kwargs)

        self.default_main_task = 'unfold'
        self.bandup_parser = BandUpArgumentParser(add_help=False)
        self.bandup_plot_parser = BandUpPlotArgumentParser(add_help=False)
        self.subparsers = self.add_subparsers(dest='main_task',
                              help='Task to be performed. (default: unfold)')
        # This is to prevent the .subparsers.add_parser method to fail,
        # attempting to init this class recursively
        self.subparsers._parser_class = argparse.ArgumentParser
        self.bandup_subparser = self.subparsers.add_parser('unfold', 
            help="Runs BandUP's main code", parents=[self.bandup_parser])
        self.bandup_plot_subparser = self.subparsers.add_parser('plot', 
            help="Plots BandUP's output files.", parents=[self.bandup_plot_parser])
        self.add_argument('-h', '--help', action='store_const', const='True',
                          help='show this help message and exit', 
                          default=argparse.SUPPRESS)

    def print_help(self, *args, **kwargs):
        print self.format_help()
        extra_help_msg = "Each task has its own help as well"
        sys.exit(0)

    def parse_known_args(self, *fargs, **fkwargs):
        # Argparse calls parse_known_args when we call parse_args, so we only need
        # to override this one

        # Finding positional args passed to program
        positional_args = []
        for arg in sys.argv[1:]:
            if(arg.startswith('-')): continue
            positional_args.append(arg)
        # Deciding whether to print the general help or not
        help_requested = len(set(['-h', '--help']).intersection(sys.argv[1:])) > 0
        pos_h_flag = float('Inf')
        if('-h' in sys.argv[1:]): pos_h_flag = sys.argv.index('-h')
        if('-help' in sys.argv[1:]): pos_h_flag=min(pos_h_flag, sys.argv.index('-help'))
        pos_subparser_choice = float('Inf')
        if('unfold' in sys.argv[1:]): pos_subparser_choice = sys.argv.index('unfold')
        if('plot' in sys.argv[1:]): 
            pos_subparser_choice = min(pos_subparser_choice, sys.argv.index('plot'))
        main_help_requested = help_requested and (pos_h_flag < pos_subparser_choice)
        if(main_help_requested):
            self.print_help()

        # Defining task (subparser) as 'unfold' if no specific task is requested
        if(fargs==(None, None) and fkwargs=={} and not help_requested and not 
           positional_args): 
            fargs=[[self.default_main_task]]
        args, unknown_args = (
            super(BandUpPythonArgumentParser, self).parse_known_args(*fargs, **fkwargs))
        args = self.filter_args(args) 

        return args, unknown_args

    def filter_args(self, args, *fargs, **fkwargs):
        if(args.main_task == 'unfold'):
            # args.argv will be passed to Popen to run BandUP's Fortran core code
            args.argv = self.bandup_parser.get_argv(args)
        elif(args.main_task == 'plot'):
            args = self.bandup_plot_parser.filter_args_plot(args)
        return args
