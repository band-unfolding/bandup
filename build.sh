#!/bin/bash

# If you have both compilers, ifort will be automatically selected unless you run
# "build.sh gfortran". If You only have gfortran, it will be automatically used.
supported_compilers='ifort gfortran'

bandup_supports_chosen_compiler=false
system_has_compiler=false
if [ ! -z $1 ]; then
    user_chose_compiler=true
    FC=${1}
    if hash $FC 2>/dev/null; then
        system_has_compiler=true
    fi 
    if [[ $supported_compilers == *$FC* ]]; then
        bandup_supports_chosen_compiler=true
    fi
else
    user_chose_compiler=false
    for compiler in $supported_compilers; do
        if hash $compiler 2>/dev/null; then
            FC=$compiler
            system_has_compiler=true
            bandup_supports_chosen_compiler=true
            break
        fi
    done
fi

if [ $system_has_compiler = true ] && [ $bandup_supports_chosen_compiler = true ]; then
    if [ $user_chose_compiler = true ]; then
        echo "BandUP will be compiled using $FC (compiler chosen by you)"
    else
        echo "BandUP will be compiled using $FC (compiler determined automatically)"
    fi
else
    if [ $user_chose_compiler = true ]; then
        if [ $bandup_supports_chosen_compiler = true ]; then
            echo "BandUP supports $FC, but you don't seem to have it available."
        else
            echo "The Fortran compiler you've chosen ($FC) is not currently supported by BandUP."
        fi
    else
        echo 'Could not find a supported Fortran compiler.'
    fi

    if [ $bandup_supports_chosen_compiler = false ]; then
        echo 'The currently Fortran compilers supported by BandUP are:'
        for compiler in $supported_compilers; do
            echo "    * $compiler"
        done
    fi

    echo 'Cannot continue. Stopping now.'
    exit
fi

if [ ${FC} = 'ifort' ]; then
    CC=icc
    OMP_FLAG='-openmp'
elif [ ${FC} = 'gfortran' ]; then
    if hash gcc-6 2>/dev/null; then
        CC=gcc-6
    elif hash gcc-5 2>/dev/null; then
        CC=gcc-5
    else
        CC=gcc
    fi
    OMP_FLAG='-fopenmp'
fi

working_dir=`pwd`
BANDUP_EXE="${working_dir}/src/python_interface/bandup"
cd ${working_dir}/src/external
    ./BandUP_configure_espresso.sh $FC $CC
    tar -xvzf check2xsf2_modules_for_BandUP.tgz
    ./BandUP_configure_spglib.sh $FC $CC $OMP_FLAG
    ./BandUP_configure_cla.sh
cd ${working_dir}/src
    make FC=$FC CC=$CC
    make clean
cd ${working_dir}
rm -f bandup
ln -s ${BANDUP_EXE} bandup

BANDUPDIR="${working_dir}"
BANDUPCONFIGDIR="${HOME}/.bandup"
BANDUPCONFIGFILE="${BANDUPCONFIGDIR}/config"
mkdir -p ${BANDUPCONFIGDIR}
now="$(date +'%Y_%m_%d-%T')"
echo "# Created by BandUP in ${now}" >| ${BANDUPCONFIGFILE}
echo "BANDUPDIR=${BANDUPDIR}" >> ${BANDUPCONFIGFILE}
