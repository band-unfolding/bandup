# Copyright (C) 2017 Paulo V. C. Medeiros
# Module used by the Python interface to BandUP to calc. orbital/atom contributions
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
import numpy as np
from scipy.sparse import coo_matrix
import time
import os
import pickle
import sys
# Imports from within the package
from .files import mkdir
from .warnings_wrapper import warnings
from .lists import int_list_to_str_range 
from .unfolding_density_op import read_unf_dens_ops
from .files import pickle_load, file_header

def formatted_orb_choice(orb_choices, supported):
    spdf = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
            'fy(3x2-y2)', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fz(x2-y2)', 'fx(x2-3y2)']
    if(type(orb_choices)==list):
        temp_list = []
        for orb in orb_choices:
            temp_list += formatted_orb_choice(orb, supported)
        temp_list = list(set(temp_list))
        rtn = [orb for orb in spdf if orb in temp_list]
        return rtn
    elif(orb_choices is None):
        return formatted_orb_choice('all', supported)

    orb_dict = {orb:[orb] for orb in supported}
    orb_dict['all'] = supported
    for general_orb in ['s', 'p', 'd', 'f']:
        orb_dict[general_orb] = [orb for orb in supported if general_orb in 
                                 orb.lower()]
    orb_dict['sp2'] = ['s', 'px', 'py']
    orb_dict['spxy'] = orb_dict['sp2']
    orb_dict['spyx'] = orb_dict['sp2']
    orb_dict['spxz'] = ['s', 'px', 'pz']
    orb_dict['spzx'] = orb_dict['spxz']
    orb_dict['spyz'] = ['s', 'py', 'pz']
    orb_dict['spzy'] = orb_dict['spyz']
    orb_dict['sp3'] = orb_dict['sp2'] + ['pz']

    orb = orb_choices.lower().strip()
    try:
        orblist = orb_dict[orb]
    except(KeyError):
       msg = 'ERROR: Invalid choice of orbital: "%s".\n'%(str(orb)) 
       msg += 'Cannot continue. Stopping now.'
       print(msg)
       sys.exit(1)
    return orblist

class KptInfo(object):
    def __init__(self, number, frac_coords, weight, nbands, nions, orbitals,
                 nkpts_in_parent_file=1,
                 nspins=1, create_ion_proj_norms_dicts=True):
        self.number = number
        self.frac_coords = frac_coords
        self.weight = weight
        self.nbands = nbands
        self.nions = nions
        self.nkpts_in_parent_file = nkpts_in_parent_file
        self.orbitals = orbitals
        self.n_orbs_per_atom = len(orbitals)
        self.n_atomic_orbs = self.nions * self.n_orbs_per_atom
        self._nspins = nspins
        self._current_ispin = 0
        self._proj_matrix = [None, None]
        self._proj_matrix_dual = [None, None]
        # orb_contr_matrices will be flexible to allow any combination of 
        # atom/orbitals to be stored. As in the proj matrices, each list item 
        # corresponds to a spin channel
        self._orb_contr_matrices = [None, None] 
        self._bands = [None, None]
        self._create_ion_proj_norms_dicts = create_ion_proj_norms_dicts

    def __getstate__(self):
        stuff_to_be_dumped_by_pickle = self.__dict__.copy()
        for key in self.__dict__:
            if('_proj_matrix' in key):
                stuff_to_be_dumped_by_pickle[key] = [None, None]
        return stuff_to_be_dumped_by_pickle
    def __setstate__(self, state):
        self.__dict__.update(state)

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
    def bands(self):
        if(self._bands[self._current_ispin] is None):
            self._bands[self._current_ispin] = (
                [{'number':None, 'ener':None, 'occ':None} for ib in xrange(self.nbands)]
            )
            if(self._create_ion_proj_norms_dicts):
                for ib in xrange(self.nbands):
                    self._bands[self._current_ispin][ib]['ion_proj_norms'] = (
                        [{orb:None for orb in self.orbitals} for
                         i in range(self.nions)]
                    )
                    self._bands[self._current_ispin][ib]['tot_orb_proj_norms'] = (
                        {orb:None for orb in self.orbitals}
                    )
        return self._bands[self._current_ispin]
    @bands.setter
    def bands(self, value):
        msg = 'Bands cannot be set manually!'
        raise AttributeError(msg)
    @property
    def proj_matrix(self):
        if(self._proj_matrix[self._current_ispin] is None):
            self._proj_matrix[self._current_ispin] = (
                np.zeros([self.n_atomic_orbs, self.nbands], dtype=complex)
            )
        return self._proj_matrix[self._current_ispin]
    #@proj_matrix.setter
    #def proj_matrix(self, value):
    #    msg = 'The projection matrix cannot be set manually!'
    #    raise AttributeError(msg)
    @proj_matrix.setter
    def proj_matrix(self, value):
        self._proj_matrix[self._current_ispin] = value
    @property
    def proj_matrix_dual(self):
        if(self.proj_matrix[self._current_ispin] is None):
            msg = 'Projection matrix has to be calculated in the first place!'
            raise ValueError(msg)
        elif(self._proj_matrix_dual[self._current_ispin] is None):
            self._proj_matrix_dual[self._current_ispin]=np.linalg.pinv(self.proj_matrix)
        return self._proj_matrix_dual[self._current_ispin]
    #@proj_matrix_dual.setter
    #def proj_matrix_dual(self, value):
    #    msg = 'The dual of the projection matrix cannot be set manually!'
    #    raise AttributeError(msg)
    @proj_matrix_dual.setter
    def proj_matrix_dual(self, value):
        self._proj_matrix_dual[self._current_ispin] = value
    @property
    def orb_contr_matrices(self):
        return self._orb_contr_matrices[self._current_ispin]
    @orb_contr_matrices.setter
    def orb_contr_matrices(self, value):
        self._orb_contr_matrices[self._current_ispin] = value

    def combined_atom_orb_index(self, iat, iorb):
        # The "alpha" variables (these alpha vars are to be renamed)
        return iat*self.n_orbs_per_atom + iorb
    def individual_atom_orb_indices(self, comb_index):
        iat = comb_index / self.n_orbs_per_atom
        iorb = comb_index % self.n_orbs_per_atom
        return iat, iorb
    def orb2iorb(self, orb):
        return self.orbitals.index(orb)
    def saveinfo(self):
        proj_mat_file_prefix = 'proj_matrices'
        pickle_file_prefix = 'KptsInfo'
        outdir='orbital_projections'
        if(self.proj_matrix is None):
            msg = 'The orbital projection matrix has to exist before it can be saved!'
            raise ValueError(msg)
        # I'll keep ne numbering starting from 1 for consistency with the Fortran output
        outfile = '%s_iKpt_%d_iSpin_%d.npz'%(proj_mat_file_prefix, 
                                             self.number, self.current_ispin+1)
        outfile = os.path.join(outdir, outfile)
        pickle_file_name = os.path.join(outdir, 
                                        '%s_iSpin_%d.pkl'%(pickle_file_prefix, 
                                                           self.current_ispin+1))
        if(self.number==1 and self.current_ispin==0):
            mkdir(outdir)
            pickle_file = open(pickle_file_name, 'w')
        else:
            pickle_file = open(pickle_file_name, 'a+b')
            mkdir(outdir, ignore_existing=True)
        np.savez(outfile, direct=self.proj_matrix, dual=self.proj_matrix_dual)
        pickle.dump(self, pickle_file, pickle.HIGHEST_PROTOCOL)
        pickle_file.close()
    def contrib_matrix_element(self, iband1, iband2, orb, 
                               selected_ion_indices=None,
                               use_dual=True, force_hermitian=False,
                              ):
        """Returns sum_{a}[<iband1| Proj(a) |iband2>], a->combined_index(selec_ion,orb)

        This routine calculates the matrix element of the following sum of projectors:
         sum_{a}[Proj(a)]; a in [combined_index(iAt,orb) for iAt in selected_ion_indices]
        evaluated between SC states with band numbers iband1 and iband2 at a given Kpt.
        We refer to the matrix representation of sum of projectors defined above as the 
        "orbital contribution matrix" for the selected atoms.

        The diagonal element (m,m) of this matrix gives the total contribution of the 
        selected atoms' "orb"-type orbitals to SC electronic band "m" at the given K.

        The off-diagonal terms are normally not used in the SC calculation itself. 
        They are needed, however, for the unfolding of this sum of projectors, which will
        then give the contribution of the selected atoms' "orb"-type orbitals to the 
        *unfolded* bands located at (k,E).

        The reason why the off-diagonal terms are needed for unfolding is explained in:
             Unfolding spinor wave functions and expectation values of general operators: 
                             Introducing the unfolding-density operator
              Paulo V. C. Medeiros, Stepan S. Tsirkin, Sven Stafstrom, and Jonas Bjork
                                      Phys. Rev. B 91, 041116(R)

        If use_dual is False, then the atomic orbital basis functions q(a) will be used
        to define the (pseudo-)projectors (they will not generally add up to identity):
                                   Proj(a) = |q(a)><q(a)|.
        The orbital contribution matrix is hermitian is this case. This is not the proper
        way of defining the projectors, and works reasonably only if the overlaps between
        the atomic orbital basis functions are negligible. Therefore, please be careful 
        when setting this option to False.

        If use_dual is True (which is the default), then the projectors are defined as:
                     Proj(a) = |Q(a)><q(a)|; <Q(a)|q(b)> = delta_{a,b}
        (i.e., {|Q(a)>} is the dual basis of {|q(a)>}). This is the proper way of doing 
        the projections. 
        Such a definition of the individual projectors will resuly in sums of projectors 
        that are, by construction, projectors as well. Mind, however, that they are 
        generally NOT hermitian. Only their sum over all values of "a" is guaranteed to
        to be a hermitian operator (because it equals identity).

        The option "force_hermitian", as the name suggests, forces the contribution 
        matrix to be hermitian. This is done as 
                       new_M[m1,m2] = 0.5*(old_M[m1,m2] + conj(old_M[m2,m1])).
        Given that the unfolding-density operator is hermitian, this guarantees that the
        unfolded contributions will be real numbers. On the other hand, however, mind
        that forcing the SC contribution matrix to be hermitian when the projectors have
        been calculated with "use_dual=True" means that the contribution matrix no longer
        represents such projectors -- neither will it generally represent any projector.
        It is easy to verify, for instance, that it will not generally be idempotent.
        """

        if(selected_ion_indices is None): 
            selected_ion_indices = xrange(self.nions)
        iorb = self.orb2iorb(orb)
        combined_at_orb_indices = (self.combined_atom_orb_index(iat, iorb=iorb) for 
                                   iat in selected_ion_indices)
        c_m_element = 0.0 + 0.0j
        if(use_dual):
            for comb_iat_iorb in combined_at_orb_indices:
                c_m_element += (self.proj_matrix_dual[iband1, comb_iat_iorb] *
                                self.proj_matrix[comb_iat_iorb, iband2])
        else:
            for comb_iat_iorb in combined_at_orb_indices:
                c_m_element += (np.conj(self.proj_matrix[comb_iat_iorb,iband1]) *
                                self.proj_matrix[comb_iat_iorb, iband2])
                
        if(force_hermitian):
            c_m_element += np.conj(self.contrib_matrix_element(iband2, iband1, orb,
                                       selected_ion_indices, use_dual, 
                                       force_hermitian=False
                                   )
                           ) 
            c_m_element *= 0.5

        return c_m_element
    def get_orbital_contribution_matrix(
        self,
        picked_orbitals='all',
        selected_ion_indices=None,
        ignore_off_diag=False,
        force_hermitian=False,
        use_dual=True,
        max_band_dE=0.1
    ): 

        picked_orbitals = formatted_orb_choice(picked_orbitals, supported=self.orbitals)
        str_selected_ion_indices = 'All'
        if(selected_ion_indices is not None):
            selec_at_indices_start_from_1 = [i+1 for i in selected_ion_indices]
            str_selected_ion_indices = (
                int_list_to_str_range(selec_at_indices_start_from_1)
            )

        row_indices = []
        col_indices = []
        entries = []
        for iband1 in range(self.nbands):
            if(force_hermitian): 
                min_iband2 = iband1
            else:
                min_iband2 = 0
            for iband2 in range(min_iband2,self.nbands):
                if(ignore_off_diag and (iband2!=iband1)): continue
                if(abs(self.bands[iband2]['ener'] - 
                       self.bands[iband1]['ener']) > max_band_dE):
                    continue
                contr = complex(0.0, 0.0)
                for orb in picked_orbitals:
                    this_orb_contr = self.contrib_matrix_element(
                                         iband1, iband2, orb, 
                                         force_hermitian=force_hermitian,
                                         use_dual=use_dual,
                                         selected_ion_indices=selected_ion_indices,
                                     )
                    if(iband1==iband2 and force_hermitian):
                        # No orbital should contribute with values outside this 
                        # interval
                        this_orb_contr = np.clip(np.real(this_orb_contr), 0.0, 1.0)
                        this_orb_contr += 0.0j
                    contr += this_orb_contr
                if(abs(contr) < 1E-2): continue
                row_indices.append(iband1)
                col_indices.append(iband2)
                entries.append(contr)
                if(force_hermitian):
                    row_indices.append(iband2)
                    col_indices.append(iband1)
                    entries.append(np.conj(contr))

        orb_contr_matrix = (
            coo_matrix((entries, (row_indices, col_indices)), 
                       shape=(self.nbands, self.nbands)).tocsr()
        )
        return orb_contr_matrix 


def get_unfolded_orb_projs(args, clip_contributions=False, verbose=False):

    spin_channel = args.spin_channel
    selected_ion_indices=args.atom_indices

    orb_projs_dir = os.path.join(args.results_dir, 'orbital_projections')
    unf_dens_ops_file = os.path.join(args.results_dir,
                                     'unfolding_density_operator_symm_avgd.dat')
    sckpts_info_file = os.path.join(orb_projs_dir, 
                                    'KptsInfo_iSpin_%d.pkl'%(spin_channel))
    orb_dependedt_unfolded_ebs_file = os.path.join(args.results_dir,
                                                   'unfolded_EBS_orb_dependent.dat')

    if(verbose):
        print('Reading unfolding-density operators...')
        print('    * File: %s'%(os.path.basename(unf_dens_ops_file)))
    unf_dens_ops = read_unf_dens_ops(filename=unf_dens_ops_file) 
    if(verbose):
        print('Done.')
        
    sckpts_info_without_projs = list(pickle_load(sckpts_info_file))

    sckpt_still_needed = [0 for i_sckpt in sckpts_info_without_projs]
    for ipckpt, unf_dens_ops_at_a_pckpt in enumerate(unf_dens_ops):
        i_sckpt = unf_dens_ops_at_a_pckpt[0].folding_sckpt_number-1
        sckpt_still_needed[i_sckpt] += 1

    sckpts_info_full = {}
    picked_orbitals = formatted_orb_choice(
                          args.orbs, 
                          supported=sckpts_info_without_projs[0].orbitals
                      )
    str_selected_ion_indices = 'All'
    if(selected_ion_indices is not None):
        # Internally we use python indexing, but externally the user sees fortran-style
        # indices. This is to keep consistency with the main BandUP (fortran) code output
        ion_indices_plus_one = [i+1 for i in selected_ion_indices]
        str_selected_ion_indices = int_list_to_str_range(ion_indices_plus_one)

    delta_N_times_orb_weights = []
    pkpt_indices = []
    energy_indices = []
    if(verbose):
        print('Parsing projections:')
        print('    * PickedAtomIndicesSC: %s'%(str_selected_ion_indices))
        print('    * OrbitalProjector: %s'%('+'.join(picked_orbitals)))
        print('')
    for ipckpt, unf_dens_ops_at_a_pckpt in enumerate(unf_dens_ops):
        i_sckpt = unf_dens_ops_at_a_pckpt[0].folding_sckpt_number-1
        if(verbose):
            msg = 'PcKpt #%d (--> ScKpt #%d), spin %d: '%(ipckpt+1, i_sckpt+1, 
                                                          spin_channel)
            print(msg)
        try:
            sckpt_info = sckpts_info_full[i_sckpt]
            sckpt_info.select_spin(spin_channel-1)
            if(verbose):
                print('    * Orb. contrib. matrix already calcd. Reusing.')
        except(KeyError):
            if(verbose):
                print('    * ScKpt requested for the first time. Loading info:')
            sckpt_info = sckpts_info_without_projs[i_sckpt]
            sckpt_info.select_spin(spin_channel-1)
            orb_projs_file = os.path.join(orb_projs_dir, 
                                          'proj_matrices_iKpt_%d_iSpin_%d.npz'%(
                                           sckpt_info.number, spin_channel))
            if(verbose):
                print('        > Reading direct and dual projection matrices...')
            with np.load(orb_projs_file) as data:
                sckpt_info.proj_matrix = data['direct']
                sckpt_info.proj_matrix_dual = data['dual']
            if(verbose):
                print('        > Calculating orbital contribution matrix...')
            sckpt_info.orb_contr_matrix = sckpt_info.get_orbital_contribution_matrix(
                                              picked_orbitals=picked_orbitals,
                                              selected_ion_indices=selected_ion_indices,
                                          )
            sckpts_info_full[i_sckpt] = sckpt_info

        for i_udo, unf_dens_op in enumerate(unf_dens_ops_at_a_pckpt):
            unf_dens_op.select_spin(sckpt_info.current_ispin)
            clip_interval = None
            if(clip_contributions): 
                clip_interval = (0.0, unf_dens_op.unfolded_N)
            unfolded_op_val = unf_dens_op.unfold(
                                  sckpt_info.orb_contr_matrix, 
                                  discard_imag=True,
                                  multiply_by_N=True, 
                                 clip_interval=clip_interval,
                                  verbose=verbose
                              )
            pkpt_indices.append(unf_dens_op.pckpt_number-1)
            energy_indices.append(unf_dens_op.iener-1)
            delta_N_times_orb_weights.append(unfolded_op_val)

        sckpt_still_needed[i_sckpt] -= 1
        if(not sckpt_still_needed[i_sckpt]):
            print('    * ScKpt #%d no longer needed. Unloading info.'%(i_sckpt+1))
            del sckpts_info_full[i_sckpt]

    n_pckpt = len(unf_dens_ops)
    nener = unf_dens_ops[0][0].nener_parent_grid
    unfolded_values_in_k_E_grid = coo_matrix((delta_N_times_orb_weights,
                                              (pkpt_indices, energy_indices)),
                                             shape=(n_pckpt,nener)).tocsr()

    if(verbose):
        print('Writing file "%s"...'%(orb_dependedt_unfolded_ebs_file))
    with open(orb_dependedt_unfolded_ebs_file, 'w') as f:
        msg = ['# Orbital/atom-decomposed unfolded band structure',
               '# nAtomsSC: %d'%(sckpts_info_without_projs[0].nions),
               '# PickedAtomIndicesSC: %s'%(str_selected_ion_indices),
               '# OrbitalProjector: %s\n'%('+'.join(picked_orbitals)),
               '# SpinChannel: %d\n'%(unf_dens_ops_at_a_pckpt[0].current_ispin+1)
        ]
        f.write(file_header(msg))
        msg = '#KCoord #E-E_Fermi #{N*Unf[<Proj>]}(k,E)\n'
        f.write(msg)
        for ikpt in range(n_pckpt):
            for iener in range(nener):
                kpt_coord = unf_dens_ops[ikpt][0].pckpt_coord_in_plot
                energy = unf_dens_ops[0][0].parent_energy_grid(iener+1)
                f.write('%.4f  %.4f   %.4f\n'%(
                        kpt_coord, energy, unfolded_values_in_k_E_grid[ikpt, iener]))
    if(verbose):
        print('Done.')
