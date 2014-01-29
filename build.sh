#!/bin/bash

working_dir=`pwd`
bin_folder="${working_dir}/BandUP_bin"
mkdir -p ${bin_folder}
rm -f ${bin_folder}/BandUP.x

cd ${working_dir}/src/external
    ./BandUP_configure_spglib.sh
cd ${working_dir}/utils/pre_unfolding/get_SCKPTS_pre_BandUP
    make
    make clean
cd ${working_dir}/src
    make
    make clean
    mv -f BandUP.x  ${bin_folder} 
cd ${working_dir}
