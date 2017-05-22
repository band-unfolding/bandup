# Copyright (C) 2017 Paulo V. C. Medeiros
# Defines and operates the unfolding-density operator and related quantities
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
from scipy.sparse import coo_matrix
import numpy as np
import itertools
# Imports from within the package
from .warnings_wrapper import warnings

class UnfDensOp(object):
    def __init__(
            self, pckpt_number, 
            pckpt_cart_coords, pckpt_frac_coords_scrl, pckpt_frac_coords_pcrl,
            pckpt_coord_in_plot, folding_sckpt_number, folding_sckpt_cart_coords,
            folding_sckpt_frac_coords_scrl, nbands, nener_parent_grid, 
            emin_parent_grid, emax_parent_grid, iener, 
            energy, unfolded_N,
            row_indices, col_indices, entries, nspins=1
        ):
        self.nbands = nbands
        self.pckpt_number = pckpt_number
        self.pckpt_cart_coords = pckpt_cart_coords
        self.pckpt_frac_coords_scrl = pckpt_frac_coords_scrl
        self.pckpt_frac_coords_pcrl = pckpt_frac_coords_pcrl
        self.pckpt_coord_in_plot = pckpt_coord_in_plot
        self.folding_sckpt_number = folding_sckpt_number
        self.folding_sckpt_cart_coords = folding_sckpt_cart_coords
        self.folding_sckpt_frac_coords_scrl = folding_sckpt_frac_coords_scrl
        self.emin_parent_grid = emin_parent_grid
        self.emax_parent_grid = emax_parent_grid
        self.nener_parent_grid = nener_parent_grid
        self.dE_parent_grid = ((emax_parent_grid - emin_parent_grid) / 
                               float(nener_parent_grid - 1))
        self.iener = iener
        self.energy = energy
        self.unfolded_N = unfolded_N
        self._nspins = nspins
        self._current_ispin = 0
        # Loading sparse matrix elements
        extended_entries = []
        extended_row_indices = []
        extended_col_indices = []
        for irow, icol, val in zip(row_indices,col_indices,entries):
            extended_row_indices.append(irow)
            extended_col_indices.append(icol)
            if(irow==icol): 
                extended_entries.append(np.real(val) + 0.0j)
            else:
                extended_entries.append(val)
                # Enforcing hermitian operator
                extended_row_indices.append(icol)
                extended_col_indices.append(irow)
                extended_entries.append(np.conj(val))
        self.csr_matrix = coo_matrix((extended_entries, 
                                      (extended_row_indices, extended_col_indices)),
                                     shape=(nbands, nbands)).tocsr()
    def parent_energy_grid(self, iener):
        ener = None
        if(0<iener<=self.nener_parent_grid):
            ener = self.emin_parent_grid + float(iener-1)*self.dE_parent_grid
        return ener

    @property
    def nspins(self):
        return self._nspins
    @nspins.setter
    def nspins(self, nspins):
        if(nspins in [1,2]):
            self._nspins = nspins
        else:
            if(nspins < 1):
                msg = 'Cannot chose nspins={0}<1! Setting nspins=1.'.format(nspins)
                self._nspins = 1
            else:
                msg = 'Cannot chose nspins={0}>2! Setting nspins=2.'.format(nspins)
                self._nspins = 2
            warnings.warn(msg)
    @property
    def current_ispin(self):
        return self._current_ispin
    @current_ispin.setter
    def current_ispin(self, ispin):
        if(0 <= ispin < self._nspins):
            self._current_ispin = ispin
        else:
            if(ispin < 0):
                msg = 'Cannot chose ispin={0}<0! Setting ispin=0.'.format(ispin)
                self._current_ispin = 0 
            else:
                msg = 'Cannot chose ispin={0}>{1}! Setting ispin={1}.'.format(
                      ispin, self._nspins-1)
                self._current_ispin = self._nspins-1 
            warnings.warn(msg)
    def select_spin(self, ispin):
        self.current_ispin = ispin
    def switch_spin_channel(self):
        if(self._nspins==1):
            warnings.warn('There is only 1 spin channel! Nothing has been changed.')
            return
        self._current_ispin = (self._current_ispin + 1) % self._nspins
    @property
    def trace(self):
        return sum(self.csr_matrix.diagonal())
    @property
    def gen_spec_weights(self):
        '''Generalised spectral weights'''
        return self.unfolded_N * self.csr_matrix 
    def unfold(self, general_operator, multiply_by_N=False,
               discard_imag=False, clip_interval=None, verbose=False):
        if(np.shape(self.csr_matrix) != np.shape(general_operator)):
            raise ValueError('Incompatible matrix shapes!')
        if((clip_interval is not None) and (not discard_imag)):
            msg = 'clip_interval can only be used if discard_imag=True.'
            raise ValueError(msg)
        temp_matrix = self.csr_matrix * general_operator
        unfolded_op_val = sum(temp_matrix.diagonal())
        if(multiply_by_N): 
            unfolded_op_val *= self.unfolded_N
        if(discard_imag):
            #if(verbose and abs(unfolded_op_val)>=1E-2):
            #    if(abs(np.imag(unfolded_op_val)) / abs(unfolded_op_val) > 0.1):
            #        msg = 'Ignoring imaginary part of unfolded value', unfolded_op_val
            #        warnings.warn(msg)
            unfolded_op_val = np.real(unfolded_op_val)
            if(clip_interval is not None):
                unfolded_op_val = np.clip([unfolded_op_val], 
                                          clip_interval[0], clip_interval[1])[0]
        return unfolded_op_val

def read_unf_dens_ops(
    filename='unfolding_density_operator_symm_avgd.dat',
    verbose=False,
    ):

    unf_dens_ops = []
    with open(filename, 'r') as udof:
        # Combining itertools and the "next" file method to have both 
        # current and next lines yiedled simultaneously
        f, f_one_step_ahead = itertools.tee(udof)
        f_one_step_ahead.next()
        f_one_step_ahead = itertools.chain(f_one_step_ahead,[None])

        next_line_contains_mat_el = False
        for iline, (line, next_line) in enumerate(itertools.izip(f, f_one_step_ahead)):
            if('SpinChannel' in line): spin_channel = int(line.split('=')[-1])
            elif('nScBands' in line): nbands = int(line.split('=')[-1])
            elif('emin' in line): emin_parent_grid = float(line.split('=')[1].split()[0])
            elif('emax' in line): emax_parent_grid = float(line.split('=')[1].split()[0])
            elif('nEner' in line): nener_parent_grid = int(line.split('=')[1])
            elif('PcKptNumber' in line):
                pckpt_number = int(line.split('=')[-1])
                if(pckpt_number>1): 
                    unf_dens_ops.append(unf_dens_ops_at_a_pckpt)
                unf_dens_ops_at_a_pckpt = []
            elif('PcKptCartesianCoords' in line):
                pckpt_cart_coords = map(float, line.split('=')[-1].split()[:3]) 
            elif('PcKptFractionalCoords' in line and 'SCRL' in line):
                pckpt_frac_coords_scrl = map(float, line.split('=')[-1].split()[:3]) 
            elif('PcKptFractionalCoords' in line and 'PCRL' in line):
                pckpt_frac_coords_pcrl = map(float, line.split('=')[-1].split()[:3]) 
            elif('PcKptLinearCoordsInBandPlot' in line):
                pckpt_coord_in_plot = float(line.split('=')[-1].split()[0])
            elif('Folds into ScKptNumber' in line):
                folding_sckpt_number = int(line.split('=')[-1])
            elif('ScKptCartesianCoords' in line):
                folding_sckpt_cart_coords = map(float, line.split('=')[-1].split()[:3])
            elif('ScKptFractionalCoords' in line and 'SCRL' in line):
                folding_sckpt_frac_coords_scrl = map(float, 
                                                     line.split('=')[-1].split()[:3])
            elif('iEnerPC' in line):
                iener = int(line.split()[3])
                energy = float(line.split()[6])
                unfolded_N = float(line.split()[10])
                row_indices = []
                col_indices = []
                entries = []
            elif(('m1' in line) and ('m2' in line) and 
                 ('GenSpecWeight' in line) and ('UnfDensOp' in line)):
                next_line_contains_mat_el = True
            elif(next_line_contains_mat_el):
                lsplit = line.replace(')','').replace('(','').replace(',','').split()
                irow, icol = (i-1 for i in map(int, lsplit[:2]))
                val = float(lsplit[4]) + 0.0j
                if(icol!=irow):
                    val += 1.0j * float(lsplit[5])
                row_indices.append(irow)
                col_indices.append(icol)
                entries.append(val)
                if(next_line is None):
                    next_line_contains_mat_el = False
                elif(next_line.startswith('#') or not next_line.strip()): 
                    next_line_contains_mat_el = False
                if(not next_line_contains_mat_el):
                    unf_dens_op = UnfDensOp(
                                            pckpt_number=pckpt_number,
                                            pckpt_cart_coords=pckpt_cart_coords, 
                                            pckpt_frac_coords_scrl=pckpt_frac_coords_scrl, 
                                            pckpt_frac_coords_pcrl=pckpt_frac_coords_pcrl,
                                            pckpt_coord_in_plot=pckpt_coord_in_plot, 
                                            folding_sckpt_number=folding_sckpt_number, 
                                            folding_sckpt_cart_coords=(
                                                folding_sckpt_cart_coords
                                            ),
                                            folding_sckpt_frac_coords_scrl=(
                                                folding_sckpt_frac_coords_scrl
                                            ),
                                            nspins=spin_channel, 
                                            nbands=nbands, 
                                            nener_parent_grid=nener_parent_grid,
                                            emin_parent_grid=emin_parent_grid, 
                                            emax_parent_grid=emax_parent_grid, 
                                            iener=iener, 
                                            energy=energy, unfolded_N=unfolded_N,
                                            row_indices=row_indices, 
                                            col_indices=col_indices, 
                                            entries=entries
                                  )
                    unf_dens_op.select_spin(spin_channel-1) 
                    unf_dens_ops_at_a_pckpt.append(unf_dens_op)
        unf_dens_ops.append(unf_dens_ops_at_a_pckpt)

    return unf_dens_ops
