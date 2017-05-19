# Copyright (C) 2017 Paulo V. C. Medeiros
# Routines used by the build script
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
from __future__ import print_function
from subprocess import Popen, PIPE, STDOUT
import shutil
import os
from .files import mkdir, rmdir
from .warnings_wrapper import warnings, WarningError

def system_has_compiler(compiler):
    try:
        compiler_call = Popen([compiler], stdout=PIPE, stderr=STDOUT)
        return True
    except(OSError):
        return False
def get_gcc():
    comp = 'gcc'
    for v in range(2,7)[::-1]:
        new_comp = 'gcc-%d'%(v)
        if(system_has_compiler(new_comp)):
            comp = new_comp
            break
    return comp
def compatible_c(fortcomp):
    if('ifort' in fortcomp): return 'icc'
    elif('gfortran' in fortcomp): return get_gcc()
    else: return None
def compatible_omp(fortcomp):
    if('ifort' in fortcomp): return '-openmp'
    elif('gfortran' in fortcomp): return '-fopenmp'
    else: return None
def doc_gen(bandup_dir, remove=False, verbose=False):
    os.environ['BANDUPDIR'] = bandup_dir
    main_doc_dir = os.path.join(bandup_dir, 'doc')
    outdir = os.path.join(main_doc_dir, 'source_code')
    os.environ['DOXIGENOUTDIR'] = outdir
    config_file = os.path.join(bandup_dir, 'src', 'Doxyfile')
    if(remove):
        if(verbose): print('Removing dir %s'%(main_doc_dir))
        rmdir(main_doc_dir)
    else:
        try:
            if(verbose): print('Creating source code documentation')
            mkdir(outdir, ignore_existing=True)
            doxygen_call = Popen(['doxygen', config_file], stdout=PIPE, stderr=STDOUT)
            out, err = doxygen_call.communicate()
        except(OSError):
            warnings.warn('Could not create source code documentation!')
def assert_valid_compiler(args, supported_fortran_compilers):
    user_chose_compiler = args.compiler is not None
    if(user_chose_compiler):
        if(not system_has_compiler(args.compiler)):
            msg = 'The compiler you chose (%s) is not available.'%(args.compiler)
            warnings.warn(msg, category=WarningError)
        print('BandUP will be compiled using %s (compiler chosen by you)'%(
              args.compiler))
    else:
        for comp in supported_fortran_compilers:
            if(system_has_compiler(comp)):
                args.compiler = comp
                break
        if(args.compiler is None):
            msg = 'Could not find a supported Fortran compiler.'
            warnings.warn(msg, category=WarningError)
        print('BandUP will be compiled using %s (compiler determined automatically)'%(
              args.compiler))
    return args
