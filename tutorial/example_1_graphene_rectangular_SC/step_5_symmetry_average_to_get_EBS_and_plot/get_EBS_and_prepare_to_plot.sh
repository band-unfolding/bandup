#!/bin/bash

working_dir=`pwd`
cd ../../../utils/post_unfolding/get_EBS_graphene_from_rectangular_SCBZ/
EBS_scripts_folder=`pwd`
cd ${working_dir}
cd ../../../utils/post_unfolding/plot/
plot_script=`pwd`/plot_unfolded_band_structure.py
cd ${working_dir}


unfolded_bands_folder="${working_dir}/../step_4_perform_unfolding_for_each_direction_in_step_3"
EBS_script="${EBS_scripts_folder}/get_EBS_direcs_K-G_G-M_M-K_graphene_rectangular_SCBZ.py"
EBS_perp_to_KG_script="${EBS_scripts_folder}/get_EBS_direc_perp_to_K-G_and_touching_K_graphene_rectangular_SCBZ.py"

# Getting EBSs
EBS_folder="${working_dir}/EBS_out"
rm -rf ${EBS_folder} && mkdir -p ${EBS_folder}
cp $EBS_script ${EBS_folder}/ebs.py
cp $EBS_perp_to_KG_script ${EBS_folder}/ebs_perp_to_kg.py
for dir in G-M1  G-M2  K1-G  K2-G  M1-K1  M2-K1  M2-K2 perp_to_K1-G_and_touching_K1 perp_to_K2-G_and_touching_K2
do
    unfolded_bands_file="${unfolded_bands_folder}/band_unfolding_onto_direc_${dir}/unfolded_band_structure.dat"
    ln -s ${unfolded_bands_file} ${EBS_folder}/unfolded_band_structure_${dir}.dat
done
cd $EBS_folder
    ./ebs.py > /dev/null && rm -f ebs.py
    ./ebs_perp_to_kg.py > /dev/null && rm -f ebs_perp_to_kg.py
    rm -f unfolded_band_structure*.dat
    tar -cvzf EBSs_for_each_pcbz_direction.tgz G-M.dat M-K.dat K-G.dat > /dev/null
cd $working_dir

# Prparing a folder with all you need to plot the EBSs
inputs_for_plot_folder="${working_dir}/input_files_plot"
plot_folder="${working_dir}/out_plot"
rm -rf $plot_folder && mkdir -p $plot_folder
cp "${inputs_for_plot_folder}"/* ${plot_folder}
cp  $plot_script "${plot_folder}/plot_unfolded_band_structure.py"
ln -s ${EBS_folder}/EBS_K-G_G-M_M-K.dat ${plot_folder}/EBS_K-G_G-M_M-K.dat
ln -s ${EBS_folder}/EBS_perp_to_K-G_and_touching_K.dat ${plot_folder}/EBS_perp_to_K-G_and_touching_K.dat

# End of script
