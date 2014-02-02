#!/bin/bash
#SBATCH --nodes 34
#SBATCH -A liu5 -p green
#SBATCH -J  GCu111Band
#SBATCH -t 1-00:00:00
# Path to VASP EXE
vasp=$vasp_half

working_dir=`pwd`
vasp_files="${working_dir}/vasp_input_files"


# Path to where the CHGCAR, POSCAR and OUTCAR files for the self-consistent calculations are
sc_calc_folder="${working_dir}/../step_1_get_converged_CHGCAR"
CHGCAR=${sc_calc_folder}'/CHGCAR'
POSCAR=${sc_calc_folder}'/POSCAR'
NGX=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[5]}'`
NGY=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`
NGZ=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[11]}'`
NGXF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[4]}'`
NGYF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[6]}'`
NGZF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`

for direc in "K-G_G-M_M-K"  "perp_to_K-G_and_touching_K"
do
    folder_for_current_direc="to_unfold_onto_pcbz_direc_${direc}"
    mkdir -p ${folder_for_current_direc}

    kpts_file_for_current_direc="${working_dir}/../step_2_get_kpts_to_be_used_in_the_SC_band_struc_calcs/KPOINTS_supercell_to_unfold_onto_pcbz_direc_${direc}.out"

    cp $POSCAR ${folder_for_current_direc}/POSCAR
    ln -s $CHGCAR ${folder_for_current_direc}/CHGCAR
    cp $vasp_files/POTCAR ${folder_for_current_direc}/POTCAR
    cp ${kpts_file_for_current_direc} ${folder_for_current_direc}/KPOINTS
    cp $vasp_files/INCAR  ${folder_for_current_direc}

    cd ${folder_for_current_direc}
        echo "" >> INCAR
        echo " for consistency with static (total_energy) run:" >> INCAR
        echo "NGX = ${NGX}; NGY = ${NGY}; NGZ = ${NGZ}" >> INCAR
        echo "NGXF = ${NGXF}; NGYF = ${NGYF}; NGZF = ${NGZF} " >> INCAR

        mpprun  $vasp > std_out_vasp_direc_${direc}

        shopt -s extglob
        eval "rm -f !(WAVECAR|OUTCAR|POSCAR|KPOINTS|INCAR|POTCAR|std_out_vasp_direc_${direc})"
    cd ..
done

#End of script
