#! /usr/bin/env python
from scipy.sparse import coo_matrix
import numpy as np
from bandupy.unfolding_density_op import read_unf_dens_ops

unf_dens_ops_file = 'unfolding_density_operator_symm_avgd.dat'
print 'Reading unfolding-density operators...'
print '    * File: %s'%(unf_dens_ops_file)
unf_dens_ops = read_unf_dens_ops(filename=unf_dens_ops_file) 
print 'Done.' 

orb_contr_file = 'orbital_contribution_matrix.dat'
print 'Reading orbital contribution matrix...'
print '    * File: %s'%(orb_contr_file)
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
            current_kpt_number = int(line.split('=')[-1])
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
print 'Done.'

delta_N_times_orb_weights = []
pkpt_indices = []
energy_indices = []
for ipckpt, unf_dens_ops_at_a_pckpt in enumerate(unf_dens_ops):
    for unf_dens_op in unf_dens_ops_at_a_pckpt:
        orb_contr_matrix = orb_contr_matrices[unf_dens_op.folding_sckpt_number-1]
        unfolded_op_val = unf_dens_op.unfold(orb_contr_matrix, discard_imag=True,
                                             multiply_by_N=True, verbose=True)
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
    msg = '#KptCoord #E-E_Fermi #OrbitalResolved_N(k,E)\n'
    f.write(msg)
    for ikpt in range(n_pckpt):
        for iener in range(nener):
            kpt_coord = unf_dens_ops[ikpt][0].pckpt_coord_in_plot
            energy = unf_dens_ops[0][0].parent_energy_grid(iener+1)
            f.write('%.4f  %.4f   %.4f\n'%(
                    kpt_coord, energy, unfolded_values_in_k_E_grid[ikpt, iener]))
print 'Done.'
