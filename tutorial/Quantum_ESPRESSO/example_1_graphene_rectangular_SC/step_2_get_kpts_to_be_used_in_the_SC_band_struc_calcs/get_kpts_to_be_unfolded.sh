#!/bin/bash

working_dir=`pwd`
bandup_folder="${HOME}/work/band_unfolding/BandUP" # Modify this line to suit your folder tree
input_files_folder="${working_dir}/input_files"
exe="${bandup_folder}/utils/pre_unfolding/get_SCKPTS_pre_BandUP/get_SCKPTS_pre_BandUP.x"

cd ${input_files_folder}
    ln -s ${exe} get_SCKPTS_pre_BandUP.x
    ./get_SCKPTS_pre_BandUP.x
    mv -f KPOINTS_supercell.out ${working_dir}
    mv -f BandUP_suggestion_of_smaller_SC_based_on_your_input_SC.POSCAR ${working_dir} &> /dev/null
cd ${working_dir}
