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
import sys
import errno
from collections import OrderedDict
from scipy.constants import physical_constants
# Imports from within the package
from .defaults import defaults
from .warnings_wrapper import warnings
from .environ import working_dir


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

def assert_valid_path(path):
    try:
        path_stat = os.stat(path)
        st_size = path_stat.st_size 
    except OSError, e:
        if e.errno == errno.ENOENT:
            raise IOError('"%s" is not a valid path!' % (path))
        else:
            raise e

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
    if(args.qe):
        in_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                   os.listdir(args.self_consist_calc_dir) if fname.endswith('.in')][0]
        with open(in_file, 'r') as f:
            in_file_lines = f.readlines()
        for line in in_file_lines:
            if('prefix' in line):
                qe_prefix = line.split('=').split()[0]
        efermi_file = os.path.join(args.self_consist_calc_dir, "%.out"%(qe_prefix))
    elif(args.abinit):
        files_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                      os.listdir(args.self_consist_calc_dir) if 
                      fname.endswith('.files')][0]
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

    if(args.castep):
        with open(fpath, 'r') as f:
            flines = f.readlines()
            for line in flines:
                if("Fermi energy" in line):
                    efermi_au = float(line.split()[-1])
                    continue
        efermi = efermi_au / physical_constants["electron volt-hartree relationship"][0]

    return efermi


def create_bandup_input(args):
    origin2dest = {}
    # Energy file
    if(args.dE is None):
        # The other related options will have been garanteed to be None when
        # parsing the arguments
        efile_passed = (
            len(set(['-efile', '--energy_info_file']).intersection(sys.argv[1:]))>0
        )
        if(efile_passed):
            fpath = getattr(args, 'energy_file')
            new_fpath = None
        else:
            new_fpath = os.path.join(args.results_dir, 'energy_info.in').strip()
            fpath = os.path.join(args.results_dir, 'energy_info.in')
            if(not os.path.isfile(fpath)):
                fpath = os.path.join(working_dir, 'energy_info.in')
        fpath = fpath.strip()
        origin2dest[fpath] = {'dest':new_fpath, 'copy':True}
    else:
        # The other related options will have been garanteed not to be None when
        # parsing the arguments
        energy_info_file = os.path.join(args.results_dir, 'energy_info.in')
        energy_info_file_contents = OrderedDict([("E_Fermi",args.efermi),
                                                 ("emin",args.emin),
                                                 ("emax",args.emax),
                                                 ("dE",args.dE)])
        with open(energy_info_file, 'w') as f:
            for k,v in energy_info_file_contents.iteritems():
                f.write("%.5f  ! %s \n"%(v,k))

    # PC, SC and PC-KPT files
    input_file_args_var_names = ['pc_file', 'sc_file', 'pckpts_file']
    for arg_name in input_file_args_var_names:
        using_default = False
        if(not '-%s'%(arg_name) in sys.argv[1:]):
            using_default = True
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

    # Wavefunction symlinks
    wf_fname = None
    wf_file = None
    if(args.castep):
        wf_fname = "%s.orbitals"%(args.seed)
    elif(args.abinit and (not '-wf_file' in sys.argv[1:])):
        with open(args.files_file, 'r') as f:
            flines = f.readlines()
        try:
            wf_fname = flines[3].strip()
        except(IndexError):
            pass
    elif(args.qe):
        raise Exception('Not yet implemmented!')
    elif(not '-wf_file' in sys.argv[1:]):
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



