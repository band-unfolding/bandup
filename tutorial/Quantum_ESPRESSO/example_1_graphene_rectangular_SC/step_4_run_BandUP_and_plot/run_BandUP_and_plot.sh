#!/bin/bash
working_dir=`pwd`

# Modify or remove this to suit your needs
export ESPRESSO_TMPDIR="outdir" 
prefix='graphene_SC_exmpl1_BandUP'

# Your choice for the energy grid used for unfolding 
emin="-20.0"
emax="  5.0"
dE=0.050

# Path to folder where BandUP (and all related tools) are. 
# You may need to modify this to suit your folder tree.
bandup_folder=${bandup_folder} 


# Data from your self-consistent calc
sc_calc_folder="${working_dir}/../step_1_get_converged_charge_density"
E_Fermi=`grep 'the Fermi energy is' "${sc_calc_folder}/pwscf.out" | tail -1 | awk '{split($0,array," ")} END{print array[5]}'`

# Where the wavefunctions are
wavefunctions_calc_dir="${working_dir}/../step_3_get_SC_wavefunctions_to_be_used_for_unfolding/${ESPRESSO_TMPDIR}"

# Paths leading to BandUP inputs
pcbz_kpts_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files"
prim_cell_lattice_file_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files"
prim_cell_lattice_file="${prim_cell_lattice_file_folder}/prim_cell_lattice.in"
supercell_lattice_file="${prim_cell_lattice_file_folder}/supercell_lattice.in"
KPOINTS_prim_cell_file=${pcbz_kpts_folder}/KPOINTS_prim_cell.in

# Paths leading to BandUP executable and plotting tool script
BandUp_exe=${bandup_folder}/BandUP_bin/BandUP.x
plot_script=${bandup_folder}/utils/post_unfolding/plot/plot_unfolded_EBS_BandUP.py

# Folder where the unfolding is to be performed
band_unfolding_folder="band_unfolding_BandUP"
mkdir -p ${band_unfolding_folder}
cd ${band_unfolding_folder}

# Copying/creating needed files
rm -f KPOINTS_prim_cell.in && cp $KPOINTS_prim_cell_file KPOINTS_prim_cell.in
rm -f prim_cell_lattice.in && cp $prim_cell_lattice_file prim_cell_lattice.in
rm -f supercell_lattice.in && cp $supercell_lattice_file supercell_lattice.in
cat >energy_info.in <<!
${E_Fermi} # E-fermi
${emin}
${emax}
${dE}
!
ln -s ${wavefunctions_calc_dir} . &> /dev/null
ln -s $BandUp_exe bandup &> /dev/null

# Running BandUP
ulimit -s unlimited
# The option -qe tells BandUP that you're using Quantum ESPRESSO
./bandup -qe -prefix $prefix | tee out_BandUP_dir_${dir}.dat 

#exit # test

# Plotting results
ln -s ${plot_script} . &> /dev/null
./plot_unfolded_EBS_BandUP.py unfolded_EBS_symmetry-averaged.dat --round_cb 0 --show --save &

# Returning to the starting directory
cd ${working_dir}

# End of script
