#!/bin/bash

working_dir=`pwd`
bin_folder="${working_dir}/BandUP_bin"
mkdir -p ${bin_folder}
rm -f ${bin_folder}/BandUP.x
rm -f ${bin_folder}/bandup

cd ${working_dir}/src/external
    ./BandUP_configure_espresso.sh
    ./BandUP_configure_spglib.sh
    ./BandUP_configure_cla.sh
cd ${working_dir}/utils/pre_unfolding/get_SCKPTS_pre_BandUP
    make
    make clean
    rm -f ${bin_folder}/get_SCKPTS_pre_BandUP.x
    ln -s `pwd`/get_SCKPTS_pre_BandUP.x ${bin_folder}/get_SCKPTS_pre_BandUP.x
cd ${working_dir}/src
    make
    make clean
    mv -f BandUP.x  ${bin_folder}
    ln -s ${bin_folder}/BandUP.x  ${bin_folder}/bandup
cd ${working_dir}

BANDUP="${bin_folder}/BandUP.x"
BANDUPBINPATH="${bin_folder}"
BANDUPPLOTPATH="${working_dir}/utils/post_unfolding/plot/"
bandup_folder="${working_dir}"

rm -f ${BANDUPBINPATH}/'bandup_plot'
rm -f ${BANDUPBINPATH}/'BandUP_plot_GUI.ui'
rm -f ${BANDUPBINPATH}/'plot_unfolded_EBS_BandUP.py'
ln -s ${BANDUPPLOTPATH}/plotting_tool_GUI/BandUP_plot_GUI.pyw ${BANDUPBINPATH}/'bandup_plot'
ln -s ${BANDUPPLOTPATH}/plotting_tool_GUI/BandUP_plot_GUI.ui ${BANDUPBINPATH}/'BandUP_plot_GUI.ui'
ln -s ${BANDUPPLOTPATH}/plot_unfolded_EBS_BandUP.py ${BANDUPBINPATH}/'plot_unfolded_EBS_BandUP.py'

# Comment the next line if you want BandUP to set some environment vars
exit 

now="$(date +'%Y_%m_%d-%T')"
for var_name in 'BANDUP' 'BANDUPBINPATH' 'BANDUPPLOTPATH' 'bandup_folder'
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
