#!/bin/bash
#SBATCH --nodes 16
#SBATCH -A liu5 -p green_risk
#SBATCH -J  GBand
#SBATCH -t 1-00:00:00
# VASP version
vasp=$vasp_half

working_dir=`pwd`
vasp_files="${working_dir}/vasp_input_files"


# Path to where the CHGCAR, POSCAR and OUTCAR files for the self-consistent calculations are
sc_calc_folder="${HOME}/work/band_unfolding/BandUP_code/BandUP_V1.0/tutorial/example_1_graphene_rectangular_SC/step_1_get_converged_CHGCAR"
CHGCAR=${sc_calc_folder}'/CHGCAR'
POSCAR=${sc_calc_folder}'/POSCAR'
NGX=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[5]}'`
NGY=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`
NGZ=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[11]}'`
NGXF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[4]}'`
NGYF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[6]}'`
NGZF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`

for direc in G-M1  G-M2  K1-G  K2-G  M1-K1  M2-K1  M2-K2 perp_to_K1-G_and_touching_K1 perp_to_K2-G_and_touching_K2
do
    folder_for_current_direc="to_unfold_onto_pcbz_direc_${direc}"
    mkdir -p ${folder_for_current_direc}

    kpts_file_for_current_direc="$vasp_files/SCBZ_kpts_that_unfold_onto_the_selected_pcbz_kpts/${direc}/KPOINTS_to_be_unfolded_onto_pcbz_${direc}"

    cp $POSCAR ${folder_for_current_direc}/POSCAR
    ln -s $CHGCAR ${folder_for_current_direc}/CHGCAR
    cp $vasp_files/POTCAR ${folder_for_current_direc}/POTCAR
    cp ${kpts_file_for_current_direc} ${folder_for_current_direc}/KPOINTS
    cp $vasp_files/INCAR  ${folder_for_current_direc}

    cd ${folder_for_current_direc}
        # The following 3 lines are optional. You can delete them
        echo "" >> INCAR
        NBANDS=`echo "print (12)*2*4" | python` # Normally NBANDS = (VSC/Vpc) * n_electrons_pc gives a reasonable number of bands above EF
        echo "NBANDS = $NBANDS" >> INCAR
        # End of optional
        echo "" >> INCAR
        echo " for consistency with static (total_energy) run:" >> INCAR
        echo "NGX = ${NGX}; NGY = ${NGY}; NGZ = ${NGZ}" >> INCAR
        echo "NGXF = ${NGXF}; NGYF = ${NGYF}; NGZF = ${NGZF} " >> INCAR

        ulimit -s unlimited
        export OMP_NUM_THREADS=1
        mpprun  $vasp > std_out_vasp_direc_${direc}

        shopt -s extglob
        eval "rm -f !(WAVECAR|OUTCAR|POSCAR|KPOINTS|INCAR|POTCAR|std_out_vasp_direc_${direc})"
    cd ..
done

#End of script
