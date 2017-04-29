import numpy as np
from scipy.sparse import coo_matrix
import time
# Imports from within the package
from .warnings_wrapper import warnings

def str_range_from_list(lst):
    lst = sorted(list(set(lst)))
    diffs = np.diff(lst)
    intervals = []
    interval = [lst[0], lst[0]]
    for idiff, diff in enumerate(diffs):
        if(diff!=1):
            interval[1] = lst[idiff]
            intervals.append(tuple(interval))
            interval[0] = lst[idiff+1] 
        interval[1] = lst[idiff+1] 
    intervals.append(tuple(interval))

    to_join = []
    for interval in intervals:
        if(interval[0]==interval[1]):
            to_join.append(str(interval[0]))
        elif(interval[1]==interval[0]+1):
            to_join.append(str(interval[0]))
            to_join.append(str(interval[1]))
        else:
            to_join.append('%d-%d'%(interval[0], interval[1]))
    str_intervals = ', '.join(to_join)
    return str_intervals

def formatted_orb_choice(orb_choices):
    if(type(orb_choices)==list): 
        return orb_choices
    SUPPORTED_ORBS = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
    orb = orb_choices.lower().strip()
    if(orb in ['all', 'spd']): orblist = SUPPORTED_ORBS
    elif(orb in ['p']): orblist = ['px', 'py', 'pz']   
    elif(orb in ['d']): orblist = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2']   
    elif(orb in ['sp3', 'sp']): orblist = ['s', 'px', 'py', 'pz']   
    elif(orb in ['sp2', 'spxy', 'spyx']): orblist = ['s', 'px', 'py']   
    elif(orb in ['spxz', 'spzx']): orblist = ['s', 'px', 'pz']   
    elif(orb in ['spxy', 'spyx']): orblist = ['s', 'px', 'py']   
    elif(orb in ['spzy', 'spyz']): orblist = ['s', 'py', 'pz']   
    else: 
        orblist = [orb_choices]
    return orblist

class KptInfo(object):
    def __init__(self, number, frac_coords, weight, nbands, nions, orbitals,
                 nkpts_in_parent_file,
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
        self._bands = [None, None]
        self._create_ion_proj_norms_dicts = create_ion_proj_norms_dicts
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
    @proj_matrix.setter
    def proj_matrix(self, value):
        msg = 'The projection matrix cannot be set manually!'
        raise AttributeError(msg)
    @property
    def proj_matrix_dual(self):
        if(self.proj_matrix[self._current_ispin] is None):
            msg = 'Projection matrix has to be calculated in the first place!'
            raise ValueError(msg)
        elif(self._proj_matrix_dual[self._current_ispin] is None):
            self._proj_matrix_dual[self._current_ispin]=np.linalg.pinv(self.proj_matrix)
        return self._proj_matrix_dual[self._current_ispin]
    @proj_matrix_dual.setter
    def proj_matrix_dual(self, value):
        msg = 'The dual of the projection matrix cannot be set manually!'
        raise AttributeError(msg)
    def combined_atom_orb_index(self, iat, iorb):
        # The "alpha" variables (these alpha vars are to be renamed)
        return iat*self.n_orbs_per_atom + iorb
    def individual_atom_orb_indices(self, comb_index):
        iat = comb_index / self.n_orbs_per_atom
        iorb = comb_index % self.n_orbs_per_atom
        return iat, iorb
    def orb2iorb(self, orb):
        return self.orbitals.index(orb)
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


def write_orbital_contribution_matrix_file(
    kpts_info_list, 
    picked_orbitals='all',
    selected_ion_indices=None,
    out_file='orbital_contribution_matrix.dat',
    open_mode='w',
    ignore_off_diag=False,
    force_hermitian=False,
    use_dual=True,
    max_band_dE=0.1
): 

    if(type(kpts_info_list)==list):
        kpts_info = kpts_info_list
    else:
        kpts_info = [kpts_info_list]
    picked_orbitals = formatted_orb_choice(picked_orbitals)

    str_selected_ion_indices = 'All'
    if(selected_ion_indices is not None):
        selec_at_indices_start_from_1 = [i+1 for i in selected_ion_indices]
        str_selected_ion_indices = str_range_from_list(selec_at_indices_start_from_1)
    msg = (
    '#################################################################################\n'
    '# File produced by BandUP at %s\n'%(time.strftime('%l:%M%p %z on %b %d, %Y'))+
    '# Copyright (C) 2017 Paulo V. C. Medeiros                                        \n'
    '#################################################################################\n'
        ) 
    if(force_hermitian): 
        msg += (
    '# Non-zero upper-triangular matrix elements of the orbital contribution matrix   \n'
    '# evaluated between SC states                                                    \n'
    '# This operator is hermitian                                                     \n'
        )
    else:
        msg += (
    '# Non-zero matrix elements of the orbital contribution matrix   \n'
    '# evaluated between SC states                                                    \n'
        )
    msg += (
    '# m1 and m2 refer to SC band indices, and the matrix elements ME are in the form \n'
    '#                         ME = (Re{ME}, Im{ME})                                  \n'
    '#                                                                                \n'
    '# OrbitalProjector = %s \n'%('+'.join(picked_orbitals))+
    '# nScBands = %d                                                                  \n'
    %(kpts_info[0].nbands)+
    '# nAtoms = %d                                                                  \n'
    %(kpts_info[0].nions)+
    '# PickedAtomIndices = %s                                                       \n'
    %(str_selected_ion_indices)+
    '# TotalNumberOfKpts = %d                                                    \n'
    %(kpts_info[0].nkpts_in_parent_file)+
    '#################################################################################\n'
    )

    with open(out_file, open_mode[0].lower()) as f:
        if(open_mode[0].lower()=='w'):
            f.write(msg)
            f.write('\n')
        for sc_kpt_info in kpts_info:
            f.write('# KptNumber = %d   SpinChannel = %d \n'%(
                     sc_kpt_info.number, sc_kpt_info.current_ispin+1))
            #f.write('#     KptCartesianCoords = %.8f  %.8f  %.8f'%(
            #         sc_kpt_info.cart_coords))
            coords = sc_kpt_info.frac_coords
            f.write('#     KptFractionalCoords = %.8f  %.8f  %.8f \n'%(
                     coords[0], coords[1], coords[2]))
            f.write('#        m1      m2      ME{m1,m2} \n')
            for iband1 in range(sc_kpt_info.nbands):
                if(force_hermitian): 
                    min_iband2 = iband1
                else:
                    min_iband2 = 0
                for iband2 in range(min_iband2,sc_kpt_info.nbands):
                    if(abs(sc_kpt_info.bands[iband2]['ener'] - 
                       sc_kpt_info.bands[iband1]['ener']) > max_band_dE):
                        continue
                    if(ignore_off_diag and (iband2!=iband1)): continue
                    contr = complex(0.0, 0.0)
                    for orb in picked_orbitals:
                        this_orb_contr = sc_kpt_info.contrib_matrix_element(
                                             iband1, iband2, orb, 
                                             force_hermitian=force_hermitian,
                                             use_dual=use_dual,
                                             selected_ion_indices=selected_ion_indices,
                                         )
                        if(force_hermitian and iband1==iband2):
                            # No orbital should contribute with values outside this 
                            # interval
                            this_orb_contr = np.clip(np.real(this_orb_contr), 0.0, 1.0)
                            this_orb_contr += 0.0j
                        contr += this_orb_contr
                    if(abs(contr) < 1E-2): continue
                    if(force_hermitian and iband1==iband2):
                        clipped_contr = np.clip(np.real(contr), 0.0, 1.0) + 0.0j
                        if((np.real(contr)<-0.009) or (np.real(contr)>1.009)):
                            msg = '|<ik=%d,m=%d|OrbProj|ik=%d,m=%d>|'%(
                                  sc_kpt_info.number,iband1,sc_kpt_info.number,iband2)
                            msg += ' clipped from %.2f to %.2f'%(
                                   np.real(contr), np.real(clipped_contr))
                        contr = clipped_contr
                    fline = 10*' ' + '%d       %d    (%.2f, %.2f) \n'%(iband1+1,iband2+1,
                            np.real(contr), np.imag(contr))
                    f.write(fline)


def read_orbital_contribution_matrix_file(
        orb_contr_file='orbital_contribution_matrix.dat'
    ):

    orb_contr_matrices = []
    with open(orb_contr_file, 'r') as orb_f:
        istart_values = float('Inf')
        is_hermitian = False
        for iline, line in enumerate(orb_f):
            if('hermitian' in line.lower()):
                is_hermitian = True
            elif('nScBands' in line):
                nbands = int(line.split('=')[-1])
            elif('KptNumber' in line):
                current_kpt_number = int(line.split('=')[1].split()[0])
                current_spin_channel = int(line.split('=')[2].split()[0])
                if(current_kpt_number>1):
                    # Appending data from previously parsed kpt
                    orb_contr_matrices.append(
                        coo_matrix((entries, (row_indices, col_indices)), 
                                   shape=(nbands, nbands)).tocsr()
                    )
                istart_values = float('Inf')
                row_indices = []
                col_indices = []
                entries = []
            elif('KptFractionalCoords' in line):
                current_kpt_frac_coords = map(float, line.split('=')[-1].split())
            elif('m1      m2' in line):
                istart_values = iline + 1
            elif(iline >= istart_values):
                lsplit = line.replace(')','').replace('(','').replace(',','').split()
                irow, icol = (i-1 for i in map(int, lsplit[:2]))
                val = float(lsplit[2]) + 1.0j * float(lsplit[3])
                if(is_hermitian):
                    if(icol==irow):
                        val = np.real(val)
                    else:
                        row_indices.append(icol)
                        col_indices.append(irow)
                        entries.append(np.conj(val))
                row_indices.append(irow)
                col_indices.append(icol)
                entries.append(val)
        # Appending last matrix read
        orb_contr_matrices.append(
            coo_matrix((entries, (row_indices, col_indices)), 
                       shape=(nbands, nbands)).tocsr()
        )

    return orb_contr_matrices
