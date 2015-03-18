#!/bin/bash
working_dir=`pwd`

wf_filename='out_graphene_rect_SC_WF_WFK'
filename_out_step_1='graphene_rect_SC_chg.out'
# Your choice for the energy grid used for unfolding 
emin="-20.0"
emax="  5.0"
dE=0.050

# Path to folder where BandUP (and all related tools) are. 
# You may need to modify this to suit your folder tree.
bandup_folder=${bandup_folder} 


# Data from your self-consistent calc. Check if the out file is indeed called "pwscf.out".
sc_calc_folder="${working_dir}/../step_1_get_converged_charge_density"
E_Fermi=`grep 'Fermi (or HOMO) energy (hartree)' "${sc_calc_folder}/${filename_out_step_1}" | tail -1 | \
         awk '{split($0,array," ")} END{print array[7]}'`
# Converting EF from Hartree to eV
E_Fermi=`echo "print $E_Fermi * 27.21138505" | python` 

# Where the wavefunctions are
wavefunctions_calc_dir="${working_dir}/../step_3_get_SC_wavefunctions_to_be_used_for_unfolding"
wf_file=${wavefunctions_calc_dir}/${wf_filename}

# Paths leading to BandUP inputs
pcbz_kpts_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files"
prim_cell_lattice_file_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files"
prim_cell_lattice_file="${prim_cell_lattice_file_folder}/prim_cell_lattice.in"
supercell_lattice_file="${prim_cell_lattice_file_folder}/supercell_lattice.in"
KPOINTS_prim_cell_file="${pcbz_kpts_folder}/KPOINTS_prim_cell.in"

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
ln -s ${wf_file} . &> /dev/null
ln -s $BandUp_exe bandup &> /dev/null

# Running BandUP
ulimit -s unlimited
# The option -abinit tells BandUP that you're using ABINIT
./bandup -abinit -wf_file $wf_filename | tee out_BandUP.dat 

# Plotting results
ln -s ${plot_script} . &> /dev/null
./plot_unfolded_EBS_BandUP.py unfolded_EBS_symmetry-averaged.dat --round_cb 0 --show --save &

# Returning to the starting directory
cd ${working_dir}

# End of script
