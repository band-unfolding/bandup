#!/bin/bash
#SBATCH --nodes 4
#SBATCH -A liu5 -p green
#SBATCH -J  GraphSCE0
#SBATCH -t 1-00:00:00

# Path to ABINIT EXE
if [ "$NSC_RESOURCE_NAME" == 'triolith' ]; then
    abinit=/software/apps/abinit/7.10.2/build01/bin/abinit
else
    abinit=/software/testing/abinit/6.12.1/bin/abinit
fi

rm -f out* *.out log.abinit tmp* *.nc *_DDB *_DEN
# Run ABINIT
mpprun $abinit < graphene_rect_SC_abinit.files >& log.abinit

# Cleaning up a bit
rm -f tmp_* *_DDB *.nc

#End of script
