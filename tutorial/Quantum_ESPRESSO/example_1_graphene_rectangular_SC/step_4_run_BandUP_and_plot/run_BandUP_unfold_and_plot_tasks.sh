#!/bin/bash

# Defining command line args. These are task-specific. To get a list of options
# available for a given task, put "-h" after the task name when running BandUP. This will
# print the help for the specific task requested.
# Run BandUP with the "-h" option only to get a list of supported tasks
# If the the task is not explicitly set, then "unfold" will be selected by default.
# BandUP tries to get the Fermi energy automatically from the self-consistent calc.
# If this fails, you can pass the Fermi energy by using the option "-efermi VALUE_IN_eV"
unfolding_task_args="-qe -prefix graphene_rectangular_SC -outdir ../step_3*/outdir"
unfolding_task_args="${unfolding_task_args} -emin -20 -emax 5 -dE 0.050"
plot_task_args='-input_file unfolded_EBS_symmetry-averaged.dat --show --save'
plot_task_args="${plot_task_args} -plotdir plot --round_cb 0"

# Choose OMP_NUM_THREADS=1 if you do not want openmp parallelization. 
# I normally use OMP_NUM_THREADS=n_cores/2
export OMP_NUM_THREADS=2

# You will probably not need to change anything below this line. I do recommend, however,
# that you read and understand what is happening in this whole file. Once you have fully
# understood what is done here, you'll realize that this script is not needed at all, as
# running BandUP directly is straightforward (particularly when you have it in your PATH)

# If you compiled BandUP as recommended (using the build.sh script), then then following
# line will define the variable "BANDUPDIR". If you have not done so, then you will need
# to change this script so that "exe" points to the "bandup" executable.
# You really should compile BandUP as recommended, though.
# If the "bandup" executable is in your PATH, then you can do just exe='bandup'
source ~/.bandup/config
# This same executable will be used for all tasks (pre-unfolding, unfolding and plot)
exe="${BANDUPDIR}/bandup"

# Preparing to run BandUP's "unfold" task
task='unfold'
task_args=${unfolding_task_args}
command_to_run="${exe} ${task} ${task_args}"
# Running BandUP with the task and task-specific options requested
ulimit -s unlimited
eval $command_to_run

# Preparing to run BandUP's "plot" task
task='plot'
task_args=${plot_task_args}
command_to_run="${exe} ${task} ${task_args}"
# Running BandUP with the task and task-specific options requested
eval $command_to_run

# End of script
