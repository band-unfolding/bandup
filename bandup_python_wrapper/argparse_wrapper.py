import argparse
import os
from bandup_python_wrapper.environ import *
from collections import OrderedDict
from bandup_python_wrapper.files import *
from plot_unfolded_EBS_BandUP import (
    main as bandup_plot,
    BandUpPlotOptions,
)


# Getting the options BandUP currently supports
cla_file = os.path.join(bandup_dir, 'src', 'cla_wrappers_mod.f90')
bandup_cla = OrderedDict()
with open(cla_file, 'r') as f:
    bandup_registered_clas = []
    cla_register_args = ['key', 'description', 'kkind', 'default']
    for line in continuation_lines(f):
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
                arg_k,arg_v=line[cla_arg_start:cla_arg_end].strip().strip(',').split('=')
                arg_v = arg_v.strip().strip('"').strip("'")
                arg_definitions_in_line[arg_k] = arg_v
            bandup_registered_clas.append(arg_definitions_in_line)

# Now using argparse to configure those BandUP supported options here
parser = argparse.ArgumentParser(
             formatter_class=argparse.ArgumentDefaultsHelpFormatter,
             add_help=False
         )

for cla_dict in bandup_registered_clas:
    arg_name = '%s'%(cla_dict['key'])
    arg_help = cla_dict['description']
    if(not arg_help): arg_help = ' '
    parser_remaining_kwargs = {}
    if(cla_dict['kkind'] in 
       ['cla_int', 'cla_float', 'cla_char', 'cla_xchar', 'cla_logical']):
        parser_remaining_kwargs['default'] = cla_dict['default']
    else:
        parser_remaining_kwargs['action'] = 'store_const'
        parser_remaining_kwargs['const'] = arg_name
        parser_remaining_kwargs['metavar'] = ' '

    parser.add_argument(arg_name, help=arg_help, **parser_remaining_kwargs)

# Defining task(s) (i)
parser.add_argument('-task', choices=['plot', 'unfold', 'p', 'u', 'pu', 'up'], 
                    default='up', 
                    help='Chooses whether to unfold, plot available data, or both.')

# Help options
parser.add_argument('-hb', '--helpbandup', action='store_true', help="Help for BandUP.")
parser.add_argument('-hp', '--helpplot', action='store_true', 
                    help="Help for BandUP's plotting tool.")
parser.add_argument('-h', '--help', action='store_true',
                    help="Help for both BandUP and BandUP's plotting tool.")
bandup_args, unknown_bandup_args = parser.parse_known_args()
print vars(bandup_args)
print unknown_bandup_args
# Defining task(s) (ii)
if('p' in bandup_args.task.lower()):
    bandup_args.plot = True
if('u' in bandup_args.task.lower()):
    bandup_args.unfold = True

# Figuring out the options supported by the plot script
sysargv_for_plotting_tool = [bandup_raw['symmetry_avgd']]
# Position of 1st non-positional args
for iopt, opt in enumerate(sys.argv):
    if(opt.startswith('-')): continue
sysargv_for_plotting_tool = [arg for arg in sys.argv[iopt:] if arg not in 
                             ['-h', '-hb', '-hp']]

if(unknown_bandup_args):
    print('The following arguments are not supported by BandUP. \n' 
          'They will be passed on to the plotting tool: %s'%(
          ', '.join(unknown_bandup_args))
         )

if(bandup_args.help):
    bandup_args.helpbandup = True
    bandup_args.helpplot = True
if(bandup_args.helpbandup): 
    print '##################################################################'
    print '                      Help for BandUP                             '
    print '##################################################################'
    parser.print_help()
    sys.exit(0)
if(bandup_args.helpplot): 
    if(bandup_args.helpbandup): print '\n \n'
    plot_opts = BandUpPlotOptions(sysargv_for_plotting_tool)
    print '##################################################################'
    print "                      Help for BandUP's Plotting tool             "
    print '##################################################################'
    plot_opts.print_help()
    sys.exit(0)

# Working out which of the BandUP options have actually been passed, and
# generating a sys.argv-like list to pass to Popen
sysargv_for_bandup = []
for cla_dict in bandup_registered_clas:
    opt = cla_dict['key']
    if(opt in sys.argv):
        argparse_var_name = '%s_%s'%(opt.strip('-').strip('-'), 'bandup_main')
        sysargv_for_bandup.append(opt)
        if(cla_dict['kkind'] in
           ['cla_int', 'cla_float', 'cla_char', 'cla_xchar', 'cla_logical']):
            val = getattr(bandup_args, argparse_var_name) 
            sysargv_for_bandup.append(val)


