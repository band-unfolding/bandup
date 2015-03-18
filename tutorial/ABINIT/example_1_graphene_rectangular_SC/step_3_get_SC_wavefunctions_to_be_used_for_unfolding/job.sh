#!/bin/bash
#SBATCH --nodes 4
#SBATCH -A liu5 -p green
#SBATCH -J  GraphWF
#SBATCH -t 0-05:00:00

# Path to ABINIT EXE
if [ "$NSC_RESOURCE_NAME" == 'triolith' ]; then
    abinit=/software/apps/abinit/7.10.2/build01/bin/abinit
else
    abinit=/software/testing/abinit/6.12.1/bin/abinit
fi

# We need the charge density from step #1
rm -f out* *.out log.abinit tmp* *.nc *_DDB *_DEN *_EIG
ln -s ../step_1_get_converged_charge_density/out_graphene_rect_SC_chg_DEN .

# Run ABINIT
mpprun $abinit < graphene_rect_SC_abinit_WF.files >& log.abinit

# Cleaning up a bit
rm -f tmp_* *_DDB *.nc

#End of script
