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
import numpy as np
# Imports from within the package
from .constants import WORKING_DIR
from .orbital_contributions import KptInfo, write_orbital_contribution_matrix_file

def procar2bandup(fpath=os.path.join(WORKING_DIR, 'PROCAR'),
                  out_file='orbital_contribution_matrix.dat',
                  picked_orbitals='all', selected_ion_indices=None,
                  ignore_off_diag=False, force_hermitian=False, use_dual=True):
    read_procar(
        fpath=fpath,
        mode='convert2bandup',
        picked_orbitals=picked_orbitals,
        selected_ion_indices=selected_ion_indices,
        out_file=out_file,
        ignore_off_diag=ignore_off_diag,
        force_hermitian=force_hermitian,
        use_dual=use_dual,
    )

def read_procar(fpath=os.path.join(WORKING_DIR, 'PROCAR'),
                mode='accumulate', 
                picked_orbitals='all', 
                selected_ion_indices=None,
                out_file='orbital_contribution_matrix.dat',
                ignore_off_diag=False,
                force_hermitian=False,
                use_dual=True):

    allowed_modes = ['accumulate', 'convert2bandup']
    if(mode not in allowed_modes):
        raise ValueError('mode: Allowed values are %s'%(', '.join(allowed_modes)))

    kpts_info = []
    with open(fpath, 'r') as f:
        # Will not use readlines here because PROCAR files can often be very big
        start_norm = float('Inf')
        start_phase = float('Inf')
        end_norm = float('-Inf')
        end_phase = float('-Inf')
        for iline, line in enumerate(f):
            lsplit = line.split()
            if('# of k-points' in line):
                nkps = int(lsplit[3])
                nbands = int(lsplit[7])
                nions = int(lsplit[11])
            elif('k-point' in line and 'weight' in line): 
                # Some PROCAR files have a formatting problem in this line
                new_lsplit = []
                for i, item in enumerate(lsplit):
                    if('-' in item and item!='k-point'):
                        before, after = item.split('-')
                        if(before.strip()): new_lsplit.append(before)
                        new_lsplit.append('-%s'%(after))
                    else:
                        new_lsplit.append(item)
                lsplit = new_lsplit
                #if(int(lsplit[1])>1): break # TEST
                print 'Reading info for Kpt #%d/%d...'%(int(lsplit[1]), nkps) # TEST
                iband = 0
                kpt  = KptInfo(number=int(lsplit[1]), 
                           frac_coords=np.array(map(float, lsplit[3:6])),
                           weight = float((lsplit[8])),
                           nbands=nbands, nions=nions,
                           nkpts_in_parent_file=nkps)
            elif(("# energy" in line) and ('# occ' in line)):
                iband += 1
                kpt.bands[iband-1]['number'] = iband
                kpt.bands[iband-1]['ener'] = float(lsplit[4])
                kpt.bands[iband-1]['occ'] = float(lsplit[7])
                start_norm = iline + 3
                end_norm = start_norm + nions - 1
                start_phase = end_norm + 3
                end_phase = start_phase + 2*nions - 1
                iatom = 0
            elif(start_norm <= iline <= end_norm):
                for i, orb in enumerate(kpt.orbitals):
                    kpt.bands[iband-1]['ion_projs'][iatom][orb]['norm'] = (
                        float(lsplit[i+1])
                    )
                iatom += 1
            elif(iline==end_norm+1):
                iatom = 0
                for i, orb in enumerate(kpt.orbitals):
                    kpt.bands[iband-1]['tot_orb_projs'][orb]['norm'] = float(lsplit[i+1])
            elif(start_phase <= iline <= end_phase):
                # Constructing orbital projection matrices
                if(iline==start_phase): reading_real_part = True
                if(reading_real_part):
                    for i, orb in enumerate(kpt.orbitals[:-1]):
                        fline_value = float(lsplit[i+1])
                        #if(abs(fline_value) < 1E-3): fline_value = 0.0
                        alpha = kpt.combined_atom_orb_index(iat=iatom, iorb=i)
                        kpt.proj_matrix[alpha,iband-1] = fline_value
                else:
                    for i, orb in enumerate(kpt.orbitals[:-1]):
                        fline_value = float(lsplit[i+1])
                        #if(abs(fline_value) < 1E-3): fline_value = 0.0
                        alpha = kpt.combined_atom_orb_index(iat=iatom, iorb=i)
                        kpt.proj_matrix[alpha,iband-1] += 1.0j * fline_value
                    iatom += 1
                reading_real_part = not reading_real_part
                if(iline==end_phase and iband==nbands):
                    if(mode=='accumulate'):
                        kpts_info.append(kpt)
                    elif(mode=='convert2bandup'):
                        open_mode = 'append'
                        if(kpt.number==1): open_mode = 'write'
                        print 'Writing info for Kpt #%d/%d...'%(kpt.number, # TEST
                              kpt.nkpts_in_parent_file) # TEST
                        write_orbital_contribution_matrix_file(
                            kpt, picked_orbitals=picked_orbitals,
                            selected_ion_indices=selected_ion_indices,
                            out_file=out_file, 
                            open_mode=open_mode,
                            ignore_off_diag=ignore_off_diag,
                            force_hermitian=force_hermitian,
                            use_dual=use_dual)
                        del kpt
    if(mode=='accumulate'):
        return kpts_info
    elif(mode=='convert2bandup'):
        return None

