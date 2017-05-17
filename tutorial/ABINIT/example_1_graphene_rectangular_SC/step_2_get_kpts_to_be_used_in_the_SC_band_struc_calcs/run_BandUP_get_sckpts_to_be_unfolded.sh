#!/bin/bash

# You will probably not need to change anything in this script. I do recommend, however,
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

# Defining the task. Run BandUP with the "-h" option for a list of supported tasks.
# If the the task is not explicitly set, then "unfold" will be selected by default.
task='kpts-sc-get'

# Defining command line args. These are task-specific. To get a list of options
# available for a given task put "-h" after the task name when running BandUP. This will
# print the help for the specific task requested.
task_args=''

# Command that will be run. The syntax is always this. 
command_to_run="${exe} ${task} ${task_args}"

# Running BandUP with the task and task-specific options requested
eval $command_to_run

# End of script
