#!/bin/bash

working_dir=`pwd`
bandup_folder="${HOME}/work/band_unfolding/BandUP" # Modify this line to suit your folder tree
general_in_folder="${working_dir}/input_files"
exe="${bandup_folder}/utils/pre_unfolding/get_SCKPTS_pre_BandUP/get_SCKPTS_pre_BandUP.x"

cd ${general_in_folder}
    for direc in "L-G_G-X_X-U_K-G"
    do
        echo "Getting SCKPTS fo unfold onto pcbz direction ${direc}"
        cp KPOINTS_prim_cell_${direc}.in KPOINTS_prim_cell.in
        ${exe} -abinit
        mv -f KPOINTS_supercell.out ${working_dir}/KPOINTS_supercell_to_unfold_onto_pcbz_direc_${direc}.out
        rm KPOINTS_prim_cell.in
    done
cd ${working_dir}
