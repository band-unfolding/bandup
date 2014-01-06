#!/bin/bash

working_dir=`pwd`
general_in_folder="${working_dir}/input_files"

exe_path="${general_in_folder}"
exe="get_SCBZ_kpts_that_unfold_onto_the_selected_pcbz_kpts.x"

general_out_folder="SCBZ_kpts_that_unfold_onto_the_selected_pcbz_kpts"
rm -rf $general_out_folder && mkdir -p $general_out_folder

for dir in G-M1 G-M2 K1-G K2-G M1-K1 M2-K1 M2-K2 perp_to_K1-G_and_touching_K1 perp_to_K2-G_and_touching_K2
do
   folder=${dir}
   mkdir -p ${folder}
   cp ${general_in_folder}/pcbz_kpts_files/KPOINTS_${dir} ${folder}/KPOINTS_prim_cell.in
   cp ${general_in_folder}/prim_cell_lattice.in  ${folder}
   cp ${general_in_folder}/supercell_lattice.in  ${folder}
   cd ${folder}
   ln -s ${exe_path}/${exe} get_kpts.x
   ./get_kpts.x
   mv -f KPOINTS_supercell_not_reduced.out KPOINTS_to_be_unfolded_onto_pcbz_${dir}_not_reduced
   mv -f KPOINTS_supercell.out KPOINTS_to_be_unfolded_onto_pcbz_${dir}
   rm -f get_kpts.x
   cd ..
   rm -rf ${general_out_folder}/${folder}
   mv  ${folder} ${general_out_folder}
done
