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
    if hash gcc-5 2>/dev/null; then
        CC=gcc-5
    else
        CC=gcc
    fi
    OMP_FLAG='-fopenmp'
fi

working_dir=`pwd`
bin_folder="${working_dir}/BandUP_bin"
mkdir -p ${bin_folder}
rm -f ${bin_folder}/BandUP.x
rm -f ${bin_folder}/bandup

cd ${working_dir}/src/external
    ./BandUP_configure_espresso.sh $FC $CC
    tar -xvzf check2xsf2_modules_for_BandUP.tgz
    ./BandUP_configure_spglib.sh $FC $CC $OMP_FLAG
    ./BandUP_configure_cla.sh
cd ${working_dir}/utils/pre_unfolding/get_SCKPTS_pre_BandUP
    make FC=$FC CC=$CC
    make clean
    rm -f ${bin_folder}/get_SCKPTS_pre_BandUP.x
    ln -s `pwd`/get_SCKPTS_pre_BandUP.x ${bin_folder}/get_SCKPTS_pre_BandUP.x
cd ${working_dir}/src
    make FC=$FC CC=$CC
    make clean
    mv -f BandUP.x  ${bin_folder}
    ln -s ${working_dir}/src/python_interface/bandup  ${bin_folder}/bandup
cd ${working_dir}

BANDUP="${bin_folder}/BandUP.x"
BANDUPBINPATH="${bin_folder}"
BANDUPPLOTPATH="${working_dir}/utils/post_unfolding/plot/"
bandup_folder="${working_dir}"
BANDUPDIR="${working_dir}"

rm -f ${BANDUPBINPATH}/'bandup_plot'
rm -f ${BANDUPBINPATH}/'BandUP_plot_GUI.ui'
rm -f ${BANDUPBINPATH}/'plot_unfolded_EBS_BandUP.py'
ln -s ${BANDUPPLOTPATH}/plotting_tool_GUI/BandUP_plot_GUI.pyw ${BANDUPBINPATH}/'bandup_plot'
ln -s ${BANDUPPLOTPATH}/plotting_tool_GUI/BandUP_plot_GUI.ui ${BANDUPBINPATH}/'BandUP_plot_GUI.ui'
ln -s ${BANDUPPLOTPATH}/plot_unfolded_EBS_BandUP.py ${BANDUPBINPATH}/'plot_unfolded_EBS_BandUP.py'

now="$(date +'%Y_%m_%d-%T')"
BANDUPCONFIGDIR="${HOME}/.bandup"
BANDUPCONFIGFILE="${BANDUPCONFIGDIR}/config"
mkdir -p ${BANDUPCONFIGDIR}
echo "# Created by BandUP in ${now}" >| ${BANDUPCONFIGFILE}
echo "BANDUPDIR=${BANDUPDIR}" >> ${BANDUPCONFIGFILE}

# Comment the next line if you want BandUP to set some environment vars
exit 

for var_name in 'BANDUP' 'BANDUPBINPATH' 'BANDUPPLOTPATH' 'bandup_folder' 'BANDUPDIR'
do
    export_statement="export $var_name=${!var_name}"
    for filename in ".bashrc" ".bash_profile" ".profile" 
    do
        file="${HOME}/${filename}"
        if [ -e "${file}" ] && [ ! -L "${file}" ]
        then
            if ! grep -Fq "$export_statement" $file # if environment variable has not been exported before
            then
                # Making a backup of the file in case it needs to be changed
                backup_file="${file}.backup.before_installing_BandUP_${now}"
                if [ ! -e "${backup_file}" ]
                then
                    cp "${file}" "${backup_file}"
                    echo "" >> ${file}
                fi
                echo "${export_statement} # Set by BandUP in ${now}" >> ${file}
            fi
        fi
    done
done

# End of script
