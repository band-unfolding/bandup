#! /usr/bin/env python
import os
from scipy.sparse import coo_matrix
import numpy as np
import pickle
from bandupy.files import pickle_load
from bandupy.unfolding_density_op import read_unf_dens_ops
from bandupy.orbital_contributions import KptInfo
import time

spin_channel = 1
picked_orbitals = 'sp2'
#selected_ion_indices=None
selected_ion_indices = xrange(32) # Only C
#selected_ion_indices = range(32, 32+64) # Only Pt
#selected_ion_indices = range(32, 32+8) # Only O


unf_dens_ops_file = 'unfolding_density_operator_symm_avgd.dat'
print 'Reading unfolding-density operators...'
print '    * File: %s'%(unf_dens_ops_file)
unf_dens_ops = read_unf_dens_ops(filename=unf_dens_ops_file) 
print 'Done.' 

orb_projs_dir = '../step_3/orbital_projections'

sckpts_info_file = os.path.join(orb_projs_dir, 
                                'KptsInfo_iSpin_%d.pkl'%(spin_channel))
sckpts_info_without_projs = list(pickle_load(sckpts_info_file))
sckpts_info_full = {}

delta_N_times_orb_weights = []
pkpt_indices = []
energy_indices = []
for ipckpt, unf_dens_ops_at_a_pckpt in enumerate(unf_dens_ops):
    first_unf_dens_op = unf_dens_ops_at_a_pckpt[0]
    i_sckpt = first_unf_dens_op.folding_sckpt_number-1
    print 'PcKpt #%d (--> ScKpt #%d), spin %d: '%(ipckpt+1, i_sckpt+1, spin_channel),
    try:
        sckpt_info = sckpts_info_full[i_sckpt]
        sckpt_info.select_spin(spin_channel-1)
        print 'Orbital contrib. matrix already calculated. Reusing'
    except(KeyError):
        print '\n    * ScKpt requested for the first time. Loading info:'
        sckpt_info = sckpts_info_without_projs[i_sckpt]
        sckpt_info.select_spin(spin_channel-1)
        orb_projs_file = os.path.join(orb_projs_dir, 
                                      'proj_matrices_iKpt_%d_iSpin_%d.npz'%(
                                       sckpt_info.number, spin_channel))
        print '        > Reading direct and dual projection matrices...'
        with np.load(orb_projs_file) as data:
            sckpt_info.proj_matrix = data['direct']
            sckpt_info.proj_matrix_dual = data['dual']
        print '        > Calculating orbital contribution matrix...'
        sckpt_info.orb_contr_matrix = sckpt_info.get_orbital_contribution_matrix(
                                          picked_orbitals=picked_orbitals,
                                          selected_ion_indices = selected_ion_indices,
                                      )
        sckpts_info_full[i_sckpt] = sckpt_info

    for i_udo, unf_dens_op in enumerate(unf_dens_ops_at_a_pckpt):
        unf_dens_op.select_spin(sckpt_info.current_ispin)
        unfolded_op_val = unf_dens_op.unfold(sckpt_info.orb_contr_matrix, 
                                             discard_imag=True,
                                             multiply_by_N=True, 
                                             #clip_interval=(0, unf_dens_op.unfolded_N),
                                             verbose=True)
        pkpt_indices.append(unf_dens_op.pckpt_number-1)
        energy_indices.append(unf_dens_op.iener-1)
        delta_N_times_orb_weights.append(unfolded_op_val)

n_pckpt = len(unf_dens_ops)
nener = unf_dens_ops[0][0].nener_parent_grid
unfolded_values_in_k_E_grid = coo_matrix((delta_N_times_orb_weights,
                                          (pkpt_indices, energy_indices)),
                                         shape=(n_pckpt,nener)).tocsr()

orb_dependedt_unfolded_ebs_file = 'unfolded_EBS_orb_dependent.dat'
print 'Writing file "%s"...'%(orb_dependedt_unfolded_ebs_file)
with open(orb_dependedt_unfolded_ebs_file, 'w') as f:
    msg = '# File created by BandUP - '
    msg += 'Band Unfolding code for Plane-wave based calculations\n'
    f.write(msg)
    msg = '#KptCoord #E-E_Fermi #{N*Unfolded[<Proj>]}(k,E)\n'
    f.write(msg)
    for ikpt in range(n_pckpt):
        for iener in range(nener):
            kpt_coord = unf_dens_ops[ikpt][0].pckpt_coord_in_plot
            energy = unf_dens_ops[0][0].parent_energy_grid(iener+1)
            f.write('%.4f  %.4f   %.4f\n'%(
                    kpt_coord, energy, unfolded_values_in_k_E_grid[ikpt, iener]))
print 'Done.'
