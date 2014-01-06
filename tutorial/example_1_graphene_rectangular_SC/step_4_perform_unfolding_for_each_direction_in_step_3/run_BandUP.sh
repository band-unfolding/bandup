#!/bin/bash

export OMP_NUM_THREADS=8 # Choose 1 if you do not want openmp parallelization. I normally use n_cores/2
E_Fermi=-2.2472 
dE=0.050
emin=-22.0
emax=7.0

working_dir=`pwd`
cd ../../../../../BandUP_code/BandUP_V1.0/BandUP_bin/
BandUp_exe=`pwd`/BandUP.x
cd ${working_dir}

wavecar_calc_dir="${working_dir}/../step_3_get_SC_wavefunctions_to_be_used_for_unfolding"
pcbz_kpts_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files/pcbz_kpts_files"
prim_cell_lattice_file_folder="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/input_files"
prim_cell_lattice_file="${prim_cell_lattice_file_folder}/prim_cell_lattice.in"

for dir in G-M1  G-M2  K1-G  K2-G  M1-K1  M2-K1  M2-K2 perp_to_K1-G_and_touching_K1 perp_to_K2-G_and_touching_K2
do

current_wavecar_folder="${wavecar_calc_dir}/to_unfold_onto_pcbz_direc_${dir}"
current_band_unfolding_folder="band_unfolding_onto_direc_${dir}"
KPOINTS_prim_cell_file=${pcbz_kpts_folder}/KPOINTS_${dir}

mkdir -p ${current_band_unfolding_folder}
cd ${current_band_unfolding_folder}

rm -f KPOINTS_prim_cell.in
cp $KPOINTS_prim_cell_file KPOINTS_prim_cell.in

rm -f prim_cell_lattice.in
cp $prim_cell_lattice_file prim_cell_lattice.in

cat >energy_info.in <<!
${E_Fermi} # E-fermi
${emin}
${emax}
${dE}
!

ln -s ${current_wavecar_folder}/WAVECAR WAVECAR
rm -f exe.x
ln -s $BandUp_exe exe.x
./exe.x > out_band_unfolding.dat 
rm -f exe.x

cd ..

done

cd ${working_dir}
