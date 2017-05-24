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
import os
import shutil
import pickle
import sys
import errno
from collections import OrderedDict
from scipy.constants import physical_constants
import argparse
import time
# Imports from within the package
from .defaults import defaults
from .warnings_wrapper import warnings, WarningError
from .constants import WORKING_DIR
from .sysargv import arg_passed
from .version import __version__


def file_header(msgs=None, next_line=None):
    basic_header = '# File created by BandUP (%s) at '%(__version__())
    basic_header += '%s\n'%(time.strftime('%-H:%M UTC%z on %b %d, %Y'))
    basic_header += '# Copyright (C) 2013-2017 Paulo V. C. Medeiros\n'
    header = 85 * '#' + '\n'
    header += basic_header
    header += 85 * '#' + '\n'
    if(msgs is None):
        pass
    elif(type(msgs)==str):
        header += msgs.strip('\n') + '\n'
        header += 85 * '#' + '\n'
    else:
        for msg in msgs:
            header += msg.strip('\n') + '\n'
        header += 85 * '#' + '\n'
    if(next_line is not None):
        if(next_line=='comment'):
            header += '#\n'
        else:
            header += next_line.strip('\n') + '\n'
    return header

def mkdir(path, ignore_existing=False):
    try:
        os.makedirs(path)
    except(OSError):
        if(os.path.isdir(path) and not ignore_existing): raise
        if(not os.path.isdir(path)): raise
def rmdir(path):
    try:
        shutil.rmtree(path)
    except(OSError):
        if(os.path.isdir(path)):
            raise
def rmfile(path):
    try:
        os.remove(path)
    except(OSError):
        if(os.path.isfile(path)):
            raise

def valid_path(arg):
    arg = os.path.abspath(os.path.relpath(arg, WORKING_DIR))
    if not os.path.exists(arg):
        raise argparse.ArgumentTypeError('"%s" is not a valid path.'%(arg))
    else:
        return arg
def assert_valid_path(path):
    try:
        path_stat = os.stat(path)
        st_size = path_stat.st_size 
    except OSError as e:
        if e.errno == errno.ENOENT:
            raise IOError('"%s" is not a valid path!' % (path))
        else:
            raise e
def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
def assert_exec_exists(exec_path):
    if(not is_exe(exec_path)):
        msg = 'The path below does not point to a valid executable:\n%s\n'%(exec_path)
        msg += 'Have you compiled BandUP?'
        warnings.warn(msg, category=WarningError)

def pickle_load(filename):
    # Source: http://stackoverflow.com/questions/20716812/saving-and-loading-multiple-objects-in-pickle-file
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break

def continuation_lines(fname, marker='&'):
    """Reads data from an opened file object taking into account line continuation marks.

    This routine reads a the data from an input file object and returns its contents 
    taking into account the presence of line continuation marks. 

    The code in this routine is an adaptation of the code available at
    http://stackoverflow.com/questions/16480495/read-a-file-with-line-continuation-
    characters-in-python

    """
    continuation_flines = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.rstrip()
            while line.endswith(marker):
                line = line[:-1] + next(f).rstrip()
            continuation_flines.append(line + '\n')
    return continuation_flines


def get_efermi_fpath(args):
    # The Fermi energy is taken only from the self consistent calculation. 
    # Please ensure convergence w.r.t. k-point mesh for your self-consistent calcs!
    efermi_file = None
    if(args.qe):
        out_files = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                     os.listdir(args.self_consist_calc_dir) if fname.endswith('.out') and
                     not fname.startswith('.')]
        if(len(out_files)>1):
            msg = 'More than 1 ".out" found at "%s"!'%(args.self_consist_calc_dir)
            msg += '\nThese are: %s\n'%(', '.join(out_files))
            msg += 'The first one containing "the Fermi energy is" will be used.'
            warnings.warn(msg)
            found_efermi = False
            for out_file in out_files:
                with open(out_file, 'r') as f:
                    for line in f:
                        if('the Fermi energy is' in line):
                            found_efermi = True
                            efermi_file = out_file
                            break
                        if(found_efermi): break
        else: 
            efermi_file = out_files[0]
    elif(args.abinit):
        files_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                      os.listdir(args.self_consist_calc_dir) if 
                      fname.endswith('.files') and not fname.startswith('.')]
        if(len(files_file)>1):
            msg = 'More than 1 "files file" found at "%s"!'%(args.self_consist_calc_dir)
            msg += '\nThese are: %s'%(', '.join(files_file))
            warnings.warn(msg)
        files_file = files_file[0]
        with open(files_file, 'r') as f:
            files_file_lines = f.readlines()
        for line in files_file_lines:
            if('.out' in line):
                efermi_file = os.path.join(args.self_consist_calc_dir,line.strip())
                continue
    elif(args.castep):
        efermi_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                       os.listdir(args.self_consist_calc_dir) if 
                       fname.endswith('.bands')][0]
    else:
        efermi_file = os.path.join(args.self_consist_calc_dir,'OUTCAR')

    return efermi_file


def get_efermi(args):
    fpath = get_efermi_fpath(args)
    efermi = None

    if(os.path.basename(fpath)=='OUTCAR'):
        with open(fpath, 'r') as f:
            for line in f:
                if('E-fermi' in line):
                    efermi = float(line.split()[2])
    elif(args.castep):
        with open(fpath, 'r') as f:
            for line in f:
                if("Fermi energy" in line):
                    efermi_au = float(line.split()[-1])
        efermi = efermi_au / physical_constants["electron volt-hartree relationship"][0]
    elif(args.abinit):
        with open(fpath, 'r') as f:
            for line in f:
                if("Fermi (or HOMO) energy (hartree)" in line):
                    efermi_au = float(line.split('=')[1].split()[0])
        efermi = efermi_au / physical_constants["electron volt-hartree relationship"][0]
    elif(args.qe):
        with open(fpath, 'r') as f:
            for line in f:
                if("the Fermi energy is" in line):
                    efermi = float(line.split('is')[-1].split()[0])
    return efermi


def guess_castep_seed(args):
    guessed_seed = None
    try:
        bands_file = [os.path.join(args.self_consist_calc_dir,fname) for 
                      fname in
                      os.listdir(args.self_consist_calc_dir) if
                      fname.endswith('.bands')][0]
        guessed_seed = os.path.splitext(os.path.basename(bands_file))[0]
    except(IndexError):
        pass
    return guessed_seed


def create_bandup_input(args):
    origin2dest = {}
    # Energy file
    if(args.emin is None):
        # The option "emax" will have been guaranteed to be None when
        # parsing the arguments
        efile_passed = arg_passed('-efile') or arg_passed('--energy_info_file')
        if(efile_passed):
            fpath = getattr(args, 'energy_file')
            new_fpath = None
        else:
            new_fpath = os.path.join(args.results_dir, 'energy_info.in').strip()
            fpath = os.path.join(args.results_dir, 'energy_info.in')
            if(not os.path.isfile(fpath)):
                fpath = os.path.join(WORKING_DIR, 'energy_info.in')
        fpath = fpath.strip()
        origin2dest[fpath] = {'dest':new_fpath, 'copy':True}
    else:
        # The option "emax" will have been guaranteed not to be None when
        # parsing the arguments
        energy_info_file = os.path.join(args.results_dir, 'energy_info.in')
        energy_info_file_contents = OrderedDict([("E_Fermi",args.efermi),
                                                 ("emin",args.emin),
                                                 ("emax",args.emax),
                                                 ("dE",args.dE)])
        with open(energy_info_file, 'w') as f:
            f.write(file_header())
            for k,v in energy_info_file_contents.iteritems():
                f.write("%.5f  # %s \n"%(v,k))

    # PC, SC and PC-KPT files
    input_file_args_var_names = ['pc_file', 'sc_file', 'pckpts_file']
    for arg_name in input_file_args_var_names:
        using_default = not arg_passed('-%s'%(arg_name))
        if(not using_default):
            fpath = getattr(args, arg_name)
            new_fpath = None
        else:
            if(defaults['pre_unfolding_inputs_dir'] is None):
                raise ValueError('Default "pre_unfolding_inputs_dir" not defined') 
            default_fname = args.default_values[arg_name]
            fpath = os.path.join(defaults['pre_unfolding_inputs_dir'], default_fname)
            new_fpath = os.path.join(args.results_dir, default_fname) .strip()
        fpath = fpath.strip()
        origin2dest[fpath] = {'dest':new_fpath, 'copy':True}

    # Wavefunction symlinks.
    wf_fname = None
    wf_file = None
    if(args.castep):
        wf_fname = "%s.orbitals"%(args.seed)
    elif(args.abinit and (not arg_passed('-wf_file'))):
        with open(args.files_file, 'r') as f:
            flines = f.readlines()
        try:
            wf_fname = flines[3].strip() + '_WFK'
        except(IndexError):
            pass
    elif(args.qe):
        # For Quantum Espresso, the WFs will be searched in args.outdir
        pass
    elif(not arg_passed('-wf_file')):
        wf_fname = args.default_values['wf_file']
    if(wf_fname is not None): 
        wf_file = os.path.join(args.wavefunc_calc_dir, wf_fname).strip()
        wf_file_symlink = os.path.join(args.results_dir, wf_fname).strip()
        origin2dest[wf_file] = {'dest':wf_file_symlink, 'copy':False}

    # Creating/verifying files
    # Wavefunction files will not be copied. Symlinks will be made in this case
    for source, dest_properties in origin2dest.iteritems():
        assert_valid_path(source)
        dest = dest_properties['dest']
        to_be_copied = dest_properties['copy']
        try:
            if(not to_be_copied):
                os.symlink(source, dest)
            else:
                if(os.path.exists(dest)): raise OSError
                shutil.copy(source, dest) 
        except(TypeError):
            if(dest is None):
                # If the user has requested a file in a specific location, in which case
                # it will be accessed directly rather than through a symlink/copy 
                pass
            else:
                raise
        except(OSError):
            if(args.overwrite and (dest!=source)): 
                rmfile(source)
                os.symlink(source, dest)
                warnings.warn('The file "%s" has been overwritten!'%(
                               os.path.basename(dest))) 
            elif(args.overwrite):
                warnings.warn('File "%s" not overwritten: Source = dest!'%(
                               os.path.basename(dest)))
            else:
                warnings.warn('File "%s" already exists and will be used! '%(
                               os.path.basename(dest))+
                              'Please use "--overwrite" if you want to '+
                              'overwrite it.')



def create_bandup_plot_input(args):
    origin2dest = {}
    if(not args.running_from_GUI):
        # Energy file
        if(args.dE is None):
            # The other related options will have been garanteed to be None when
            # parsing the arguments
            efile_passed = (arg_passed('-efile') or arg_passed('--energy_info_file') or
                            arg_passed('-energy_file'))
            if(efile_passed):
                fpath = getattr(args, 'energy_file')
                new_fpath = None
            else:
                new_fpath = os.path.join(args.plotdir, 'energy_info.in').strip()
                fpath = os.path.join(args.results_dir, 'energy_info.in')
                if(not os.path.isfile(fpath)):
                    fpath = os.path.join(WORKING_DIR, 'energy_info.in')
            fpath = fpath.strip()
            origin2dest[fpath] = {'dest':new_fpath, 'copy':True}
        else:
            # The other related options will have been garanteed not to be None when
            # parsing the arguments
            energy_info_file = os.path.join(args.plotdir, 'energy_info.in')
            energy_info_file_contents = OrderedDict([("E_Fermi",args.efermi),
                                                     ("emin",args.emin),
                                                     ("emax",args.emax),
                                                     ("dE",args.dE)])
            with open(energy_info_file, 'w') as f:
                for k,v in energy_info_file_contents.iteritems():
                    f.write("%.5f  ! %s \n"%(v,k))

    # PC, SC, PC-KPT, and main input files
    input_file_args_var_names = ['pckpts_file', 'pc_file', 'input_file']
    var_names = {'pckpts_file':'kpoints_file', 'pc_file':'prim_cell_file',
                 'input_file':'input_file'}
    for arg_name in input_file_args_var_names:
        argval = getattr(args, var_names[arg_name]).strip()
        default = args.default_values[var_names[arg_name]].strip()
        using_default = argval==default
        if(not using_default):
            fpath = os.path.abspath(os.path.relpath(argval, WORKING_DIR))
            new_fpath = None
        else:
            if(not os.path.exists(args.results_dir)):
                raise ValueError('Results directory "%s" could not be found'%(
                                 args.results_dir)) 
            default_fname = args.default_values[var_names[arg_name]]
            fpath = os.path.join(args.results_dir, default)
            new_fpath = os.path.join(args.plotdir, default).strip()
        fpath = fpath.strip()
        origin2dest[fpath] = {'dest':new_fpath, 'copy':True}

    # Creating/verifying files
    # Wavefunction files will not be copied. Symlinks will be made in this case
    for source, dest_properties in origin2dest.iteritems():
        assert_valid_path(source)
        dest = dest_properties['dest']
        to_be_copied = dest_properties['copy']
        try:
            if(not to_be_copied):
                os.symlink(source, dest)
            else:
                if(os.path.exists(dest)): raise OSError
                shutil.copy(source, dest) 
        except(TypeError):
            if(dest is None):
                # If the user has requested a file in a specific location, in which case
                # it will be accessed directly rather than through a symlink/copy 
                pass
            else:
                raise
        except(OSError):
            if(args.overwrite and (dest!=source)): 
                rmfile(source)
                os.symlink(source, dest)
                warnings.warn('The file "%s" has been overwritten!'%(
                               os.path.basename(dest))) 
            elif(args.overwrite):
                warnings.warn('File "%s" not overwritten: Source = dest!'%(
                               os.path.basename(dest)))
            else:
                warnings.warn('File "%s" already exists and will be used! '%(
                               os.path.basename(dest))+
                              'Please use "--overwrite" if you want to '+
                              'overwrite it.')
