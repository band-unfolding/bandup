#!/bin/bash
#SBATCH --nodes 4 --exclusive
#SBATCH -J  BulkSiBands
#SBATCH -A liu5 -p green
#SBATCH -t 1-00:00:00


# Path to QE. You probably need to modify this!
pwscf=/software/testing/espresso/5.0.3/build01/bin/pw.x
# Modify or remove this to suit your needs
export ESPRESSO_TMPDIR="outdir" 
# Modify this to point to your pseudopotential folder
export ESPRESSO_PSEUDO=`pwd`/'../../upf_files'

ln -s ../step_1_get_converged_charge_density/${ESPRESSO_TMPDIR} .
mpprun $pwscf -input bulk_Si_pwscf_bands.in > pwscf_bands.out

#End of script
