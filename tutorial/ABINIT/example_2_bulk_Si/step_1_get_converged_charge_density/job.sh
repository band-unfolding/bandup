#!/bin/bash
#SBATCH --nodes 4
#SBATCH -A liu5 -p green
#SBATCH -J  BulkSiE0
#SBATCH -t 1-00:00:00

# Path to ABINIT EXE
if [ "$NSC_RESOURCE_NAME" == 'triolith' ]; then
    abinit=/software/apps/abinit/7.10.2/build01/bin/abinit
else
    abinit=/software/testing/abinit/6.12.1/bin/abinit
fi

rm -f out* *.out log.abinit tmp* *OUT.nc *_DDB *_DEN
# Run ABINIT
mpprun $abinit < bulk_Si_chg.files >& log.abinit

rm -f tmp* *OUT.nc *_DDB

#End of script
