#!/bin/bash
#SBATCH --nodes 4 --exclusive
#SBATCH -J Bulk_Si_SC
#SBATCH -A liu5 -p green
#SBATCH -t 1-00:00:00

# Path to QE. You probably need to modify this!
pwscf=/software/testing/espresso/5.0.3/build01/bin/pw.x
# Modify or remove this to suit your needs
export ESPRESSO_TMPDIR="./outdir" 
# Modify this to point to your pseudopotential folder
export ESPRESSO_PSEUDO=`pwd`/'../../upf_files'

mpprun $pwscf -input bulk_Si_pwscf.in > pwscf.out

#End of script
