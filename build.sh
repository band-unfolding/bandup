#!/bin/bash


working_dir=`pwd`
bin_folder="${working_dir}/BandUP_bin"
mkdir -p ${bin_folder}
rm -f ${bin_folder}/BandUP.x
cd src
    make
    make clean
    mv -f BandUP.x  ${bin_folder} 
cd ${working_dir}
