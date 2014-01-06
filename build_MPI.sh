#!/bin/bash

mpi="impi/4.0.3.008"

working_dir=`pwd`
bin_folder="${working_dir}/BandUP_bin"
mkdir -p ${bin_folder}
rm -f ${bin_folder}/BandUP_MPI.x
cd src
    module load $mpi
    make -f Makefile_MPI
    make clean -f Makefile_MPI
    mv -f BandUP_MPI.x  ${bin_folder}
cd ${working_dir}
