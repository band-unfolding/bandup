# Copyright (C) 2017 Paulo V. C. Medeiros
# Parsing VASP output
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
import numpy as np
import time
# Imports from within the package
from .constants import WORKING_DIR
from .orbital_contributions import KptInfo
from .warnings_wrapper import warnings, WarningError

class RealAndImagPartsNotInSameLineProcar(IOError):
    pass

def procar2bandup(fpath=os.path.join(WORKING_DIR, 'PROCAR'),
                  ignore_off_diag=False, force_hermitian=False, use_dual=True):

    read_procar(
        fpath=fpath,
        mode='convert2bandup',
        ignore_off_diag=ignore_off_diag,
        force_hermitian=force_hermitian,
        use_dual=use_dual,
        only_phase=True,
    )

def try_to_read_procar(fpath=os.path.join(WORKING_DIR, 'PROCAR'),
                mode='accumulate', 
                ignore_off_diag=False,
                force_hermitian=False,
                use_dual=True,
                only_phase=False,
                same_line_real_and_imag_phase=True):

    allowed_modes = ['accumulate', 'convert2bandup']
    if(mode not in allowed_modes):
        raise ValueError('mode: Allowed values are %s'%(
                         ', '.join(allowed_modes)))

    # Getting orbitals present in file
    available_orbs = []
    with open(fpath, 'r') as f:
        for line in f:
            if('ion' in line and 'tot' in line):
                available_orbs = line.split()[1:-1]
                break
    kpts_info = []
    with open(fpath, 'r') as f:
        # Will not use readlines because PROCAR files can often be very big
        start_norm = float('Inf')
        start_phase = float('Inf')
        end_norm = float('-Inf')
        end_phase = float('-Inf')
        ispin = -1
        for iline, line in enumerate(f):
            lsplit = line.split()
            if('# of k-points' in line):
                ispin += 1
                nkps = int(lsplit[3])
                nbands = int(lsplit[7])
                nions = int(lsplit[11])
                if(ispin==0 and mode=='accumulate'):
                    kpts_info = [None for ik in range(nkps)]
            elif('k-point' in line and 'weight' in line): 
                # Some PROCAR files have a formatting problem in this line,
                # where there will be no space between two numbers if the
                # second one has a negative sign. Accounting for this and 
                # producing a proper lsplit for this case:
                index_collon = line.index(":")
                lsplit = (line[:index_collon] + 
                          line[index_collon:].replace("-", " -")).split()

                kpt_number = int(lsplit[1])
                msg = "Kpt #%d/%d, spin_channel=%d: "%(kpt_number,nkps,ispin+1)
                msg += "Reading <atOrb|PW> projections..."
                print(msg)
                if(mode=='convert2bandup' or ispin==0):
                    kpt  = KptInfo(number=kpt_number, 
                               frac_coords=map(float, lsplit[3:6]),
                               weight = float((lsplit[8])),
                               nbands=nbands, nions=nions,
                               orbitals=available_orbs,
                               nkpts_in_parent_file=nkps,
                               create_ion_proj_norms_dicts=not only_phase)
                else:
                    kpt = kpts_info[kpt_number-1]
                    kpts_info[kpt_number-1] = None 
                if(ispin==1): 
                    kpt.nspins = 2
                kpt.select_spin(ispin)
                iband = 0
            elif(("energy" in line) and ('occ' in line)):
                iband += 1
                kpt.bands[iband-1]['number'] = iband
                kpt.bands[iband-1]['ener'] = float(lsplit[4])
                kpt.bands[iband-1]['occ'] = float(lsplit[7])
                start_norm = iline + 3
                end_norm = start_norm + nions - 1
                start_phase = end_norm + 3
                if(same_line_real_and_imag_phase):
                    end_phase = start_phase + nions - 1
                else:
                    end_phase = start_phase + 2*nions - 1
                iatom = 0
            elif(not only_phase and (start_norm <= iline <= end_norm+1)):
                if(iline==end_norm+1):
                    for i, orb in enumerate(kpt.orbitals):
                        kpt.bands[iband-1]['tot_orb_proj_norms'][orb]= (
                            float(lsplit[i+1])
                        )
                else:
                    for i, orb in enumerate(kpt.orbitals):
                        kpt.bands[iband-1]['ion_proj_norms'][iatom][orb] = (
                            float(lsplit[i+1])
                        )
                    iatom += 1
            elif(start_phase <= iline <= end_phase):
                # Constructing orbital projection matrices
                if(iline==start_phase): 
                    iatom = 0
                    reading_real_part = True
                comb_iat_iorbs = [kpt.combined_atom_orb_index(iat=iatom,iorb=i)
                                  for i in xrange(kpt.n_orbs_per_atom)]
                if(same_line_real_and_imag_phase):
                    try:
                        for i, comb_iat_iorb in enumerate(comb_iat_iorbs):
                            cplx_val = (float(lsplit[i+1]) + 
                                        1.0j * float(lsplit[i+2]))
                            kpt.proj_matrix[comb_iat_iorb, iband-1] = cplx_val
                        iatom += 1
                    except(IndexError):
                        raise(RealAndImagPartsNotInSameLineProcar)
                else:
                    try:
                        if(reading_real_part):
                            for i, comb_iat_iorb in enumerate(comb_iat_iorbs):
                                re_part = float(lsplit[i+1])
                                kpt.proj_matrix[comb_iat_iorb, iband-1]=re_part
                        else:
                            for i, comb_iat_iorb in enumerate(comb_iat_iorbs):
                                im_part = 1.0j * float(lsplit[i+1])
                                kpt.proj_matrix[comb_iat_iorb,iband-1]+=im_part
                            iatom += 1
                        reading_real_part = not reading_real_part
                    except(IndexError):
                        msg = "You should use LORBIT=12 to generate your "
                        msg += "PROCAR file. If you've done so, "
                        msg += "plese contact us."
                        warnings.warn(msg, category=WarningError)
                if(iline==end_phase and iband==nbands):
                    if(mode=='convert2bandup'):
                        print('    * Calculating duals and saving results...')
                        kpt.saveinfo()
                    elif(mode=='accumulate'):
                        kpts_info[kpt.number-1] = kpt
                    del kpt
    if(mode=='convert2bandup'):
        return None
    elif(mode=='accumulate'):
        return kpts_info

def read_procar(fpath=os.path.join(WORKING_DIR, 'PROCAR'),
                mode='accumulate', 
                ignore_off_diag=False,
                force_hermitian=False,
                use_dual=True,
                only_phase=False):
    try:
        try_to_read_procar(
            fpath=fpath,
            mode=mode,
            ignore_off_diag=ignore_off_diag,
            force_hermitian=force_hermitian,
            use_dual=use_dual,
            only_phase=only_phase,
            same_line_real_and_imag_phase=True
        )
    except(RealAndImagPartsNotInSameLineProcar):
        msg = "    * Failed: PROCAR file format differs from the one "
        msg += "originally assumed."
        print(msg)
        print("      Trying again assuming a different format:")
        try_to_read_procar(
            fpath=fpath,
            mode=mode,
            ignore_off_diag=ignore_off_diag,
            force_hermitian=force_hermitian,
            use_dual=use_dual,
            only_phase=only_phase,
            same_line_real_and_imag_phase=False
        )
