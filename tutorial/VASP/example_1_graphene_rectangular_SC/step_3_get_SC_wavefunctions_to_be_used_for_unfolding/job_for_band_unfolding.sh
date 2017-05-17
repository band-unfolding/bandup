#!/bin/bash

command_to_run='mpprun vasp > std_out_vasp.dat'
input_files_dir='vasp_input_files/'

# Path to where the CHGCAR, POSCAR and OUTCAR files for the self-cons. calculations are
sc_calc_folder='../step_1_*/'
CHGCAR=${sc_calc_folder}'/CHGCAR'
POSCAR=${sc_calc_folder}'/POSCAR'
# Grid stuff
NGX=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[5]}'`
NGY=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`
NGZ=`grep 'dimension x,y,z NGX =' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[11]}'`
NGXF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[4]}'`
NGYF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[6]}'`
NGZF=`grep 'dimension x,y,z NGXF=' "${sc_calc_folder}/OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`

# Preparing inputs
cp $POSCAR POSCAR
ln -s $CHGCAR CHGCAR
cp $input_files_dir/POTCAR POTCAR
cp ../step_2_*/KPOINTS_supercell.out KPOINTS
cp $input_files_dir/INCAR INCAR
# The following 5 lines are optional. You can delete them if you want
echo "" >> INCAR
# Normally NBANDS = (VSC/Vpc) * n_electrons_pc gives a reasonable number of bands 
# above EF
NBANDS=`echo "print (12)*2*4" | python` 
echo "NBANDS = $NBANDS" >> INCAR
# End of optional
echo "" >> INCAR
echo " for consistency with static (total_energy) run:" >> INCAR
echo "NGX = ${NGX}; NGY = ${NGY}; NGZ = ${NGZ}" >> INCAR
echo "NGXF = ${NGXF}; NGYF = ${NGYF}; NGZF = ${NGZF} " >> INCAR

# Running
eval $command_to_run

#End of script
