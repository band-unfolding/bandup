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
from .orbital_contributions import KptInfo


def read_procar(fpath=os.path.join(WORKING_DIR, 'PROCAR')):
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
                #if(int(lsplit[1])>1): break # TEST
                print 'Reading info for Kpt #%d...'%(int(lsplit[1])) # TEST
                iband = 0
                kpt  = KptInfo(number=int(lsplit[1]), 
                           frac_coords=np.array(map(float, lsplit[3:6])),
                           weight = float((lsplit[8])),
                           nbands=nbands, nions=nions)
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
                iline_from_start = iline - start_phase
                if(not iline_from_start%2):
                    for i, orb in enumerate(kpt.orbitals[:-1]):
                        fline_value = float(lsplit[i+1])
                        if(abs(fline_value) < 1E-3): fline_value = 0.0
                        kpt.bands[iband-1]['ion_projs'][iatom][orb]['phase'] = (
                            fline_value + 0.0j
                        )
                else:
                    for i, orb in enumerate(kpt.orbitals[:-1]):
                        fline_value = float(lsplit[i+1])
                        if(abs(fline_value) < 1E-3): fline_value = 0.0
                        kpt.bands[iband-1]['ion_projs'][iatom][orb]['phase'] += (
                            1.0j * fline_value
                        )
                    iatom += 1
                if(iline==end_phase and iband==nbands):
                    # Constructing orbital projection matrices and getting duals of projs
                    proj_matrix = np.zeros([kpt.nions*len(kpt.orbitals[:-1]), 
                                            kpt.nbands],
                                           dtype=complex)
                    for iband in range(kpt.nbands):
                        for iorb, orb in enumerate(kpt.orbitals[:-1]):
                            for iatom2 in range(kpt.nions):
                                alpha = iatom2*len(kpt.orbitals[:-1]) + iorb
                                proj_matrix[alpha,iband]+=(kpt.bands[iband]['ion_projs']
                                                          [iatom2][orb]['phase'])
                    dual_proj_matrix = np.linalg.pinv(proj_matrix)
                    for iband in range(kpt.nbands):
                        for iorb, orb in enumerate(kpt.orbitals[:-1]):
                            kpt.bands[iband]['tot_orb_projs'][orb]['phase'] = []
                            kpt.bands[iband]['tot_orb_projs'][orb]['phase_dual'] = []
                            for iatom2 in range(kpt.nions):
                                alpha = iatom2*len(kpt.orbitals[:-1]) + iorb
                                kpt.bands[iband]['tot_orb_projs'][orb]['phase'].append(
                                    proj_matrix[alpha,iband])
                                (kpt.bands[iband]['tot_orb_projs'][orb]
                                          ['phase_dual'].append(
                                    dual_proj_matrix[iband,alpha]))
                    kpts_info.append(kpt)
    return kpts_info



