#!/bin/bash

sc_calc_folder="${working_dir}/../step_1_*/"
wavefunctions_calc_dir="${working_dir}/../step_3_*/"

wf_filename='out_graphene_rect_SC_WF_WFK'
filename_out_step_1='graphene_rect_SC_chg.out'

# Data from your self-consistent calc. Check if the out filename is correct.
#E_Fermi=`grep 'Fermi (or HOMO) energy (hartree)' "${sc_calc_folder}/${filename_out_step_1}" | tail -1 | \
#         awk '{split($0,array," ")} END{print array[7]}'`
# Converting EF from Hartree to eV
#E_Fermi=`echo "print $E_Fermi * 27.21138505" | python` 
# Where the wavefunctions are
#wf_file=${wavefunctions_calc_dir}/${wf_filename}


# If you compiled BandUP as recommended (using the build.sh script), then then following
# line will define the variable "BANDUPDIR". If you have not done so, then you will need
# to change this script so that "exe" points to the "bandup" executable.
# You really should compile BandUP as recommended, though.
source ~/.bandup/config

# This same executable will be used for all tasks (pre-unfolding, unfolding and plot)
exe="${BANDUPDIR}/bandup"

# Defining the task. Run the code with the "-h" option for a list of supported tasks.
task='unfold'

# Defining command line args. These are task-specific. To get a list of options
# available for a given task put "-h" after the task name when running BandUP. This will
# print the help for the specific task requested.
task_args="-abinit -wf_file $wf_file -emin -20 -emax 5 -dE 0.050 -efermi ${E_Fermi}"

# Command that will be run
command_to_run="${exe} ${task} ${task_args}"

# Running BandUP with the task and task-specific options requested
ulimit -s unlimited
eval $command_to_run

# Plotting results
task='plot'
task_args='-input_file unfolded_EBS_symmetry-averaged.dat --round_cb 0 --show --save -plotdir plot'
# Command that will be run
command_to_run="${exe} ${task} ${task_args}"
# Running BandUP with the task and task-specific options requested
eval $command_to_run

# Returning to the starting directory
cd ${working_dir}

# End of script
