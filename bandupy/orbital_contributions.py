import numpy as np
import copy
import time
# Imports from within the package
from .warnings_wrapper import warnings
  
SUPPORTED_ORBS = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']

def formatted_orb_choice(orb_choices):
    if(type(orb_choices)==list): 
        return orb_choices
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

class KptInfo():
    def __init__(self, number, frac_coords, weight, nbands, nions):
        self.number = number
        self.frac_coords = frac_coords
        self.weight = weight
        self.nbands = nbands
        self.nions = nions
        self.orbitals = SUPPORTED_ORBS + ['tot']
        __general_proj_info = {}
        for orb in self.orbitals:
            __general_proj_info[orb] = {'norm':None, 'phase':None, 'phase_dual':None}
        __band_info = {'number':None, 'ener':None, 'occ':None,
                       'ion_projs':[copy.deepcopy(__general_proj_info) for 
                                    i in range(nions)],
                       'tot_orb_projs':copy.deepcopy(__general_proj_info)}
        self.bands = [copy.deepcopy(__band_info) for ib in range(nbands)]

    def contrib_matrix_element(self, iband1, iband2, orb,
                               use_dual=True, force_hermitian=False):
        proj1 = self.bands[iband2]['tot_orb_projs'][orb]['phase']
        if(use_dual):
            proj2 = self.bands[iband1]['tot_orb_projs'][orb]['phase_dual']
        else:
            proj2 = map(np.conj, self.bands[iband1]['tot_orb_projs'][orb]['phase'])
        c_m_element = np.dot(proj2,proj1)

        if(force_hermitian):
            c_m_element += np.conj(self.contrib_matrix_element(iband2, iband1, orb,
                                   use_dual, force_hermitian=False))
            c_m_element *= 0.5

        return c_m_element


def write_orbital_contribution_matrix_file(
    kpts_info_list, 
    picked_orbitals='all',
    out_file='orbital_contribution_matrix.dat',
    open_mode='w'
): 

    if(type(kpts_info_list)==list):
        kpts_info = kpts_info_list
    else:
        kpts_info = [kpts_info_list]
    picked_orbitals = formatted_orb_choice(picked_orbitals)

    msg = (
    '#################################################################################\n'
    '# File produced by BandUP at %s\n'%(time.strftime('%l:%M%p %z on %b %d, %Y'))+
    '# Copyright (C) 2017 Paulo V. C. Medeiros                                        \n'
    '#################################################################################\n'
    '# Non-zero upper-triangular matrix elements of the orbital contribution matrix   \n'
    '# evaluated between SC states                                                    \n'
    '# This operator is hermitian                                                     \n'
    '# m1 and m2 refer to SC band indices, and the matrix elements ME are in the form \n'
    '#                         ME = (Re{ME}, Im{ME})                                  \n'
    '#                                                                                \n'
    '# Included orbitals: %s \n'%(' '.join(picked_orbitals))+
    '# nScBands = %d                                                                  \n'
    %(kpts_info[0].nbands)+
    '#################################################################################\n'
    )

    with open(out_file, open_mode[0].lower()) as f:
        if(open_mode[0].lower()=='w'):
            f.write(msg)
            f.write('\n')
        for sc_kpt_info in kpts_info:
            f.write('# ScKptNumber = %d \n'%(sc_kpt_info.number))
            #f.write('#     ScKptCartesianCoords = %.8f  %.8f  %.8f'%(
            #         sc_kpt_info.cart_coords))
            coords = sc_kpt_info.frac_coords
            f.write('#     ScKptFractionalCoords = %.8f  %.8f  %.8f \n'%(
                     coords[0], coords[1], coords[2]))
            f.write('#        m1      m2       ME{m1,m2} \n')
            for iband1 in range(sc_kpt_info.nbands):
                for iband2 in range(iband1,sc_kpt_info.nbands):
                    contr = complex(0.0, 0.0)
                    for orb in picked_orbitals:
                        this_orb_contr = sc_kpt_info.contrib_matrix_element(
                                             iband1, iband2, orb, force_hermitian=True
                                         )
                        contr += this_orb_contr
                        if((abs(this_orb_contr)>1.009) and (iband1==iband2)): 
                            warnings.warn(
                                '|<ik=%d,m=%d|Proj[%s]|ik=%d,m=%d>| = %.2f > 1'%(
                                   sc_kpt_info.number, iband1, orb, 
                                   sc_kpt_info.number, iband2, abs(contr)))
                    if(abs(contr) < 1E-3): continue
                    fline = 10*' ' + '%d       %d    (%.5f, %.5f) \n'%(iband1+1,iband2+1,
                            np.real(contr), np.imag(contr))
                    f.write(fline)
