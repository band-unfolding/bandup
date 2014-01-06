#!/bin/bash
#SBATCH -N 11
#SBATCH -t 0-00:05:00
#SBATCH -J  MPI
# Accounts
# For kappa
##SBATCH -A liu5 -p green
##SBATCH -A liu5 -p green_risk
#SBATCH -p kappa -A liu-2013-00109-30
##SBATCH -p huge -A liu-2013-00109-30-huge
# For triolith
##SBATCH -A snic2013-1-165
##SBATCH -A snic2013-1-99

# Only 1 MPI task per node
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
##SBATCH --mail-type=ALL --mail-user=paulovcmedeiros@gmail.com

export OMP_NUM_THREADS=`echo "print int($SLURM_CPUS_ON_NODE/1.0)" | python`
file_name_suffix="${NSC_RESOURCE_NAME}_mpi_${SLURM_NNODES}_nodes_omp_${OMP_NUM_THREADS}_threads"
std_out_file="out_band_unfolding_${file_name_suffix}.dat"
results_file="unfolded_band_structure_${file_name_suffix}.dat"

rm -f ${results_file} ${std_out_file}
echo "Job ID: ${SLURM_JOBID}" > $std_out_file
echo "Running on ${NSC_RESOURCE_NAME}. Node list: ${SLURM_NODELIST}" >> $std_out_file
mpprun  ./band_unfolding_V1.0_MPI.x >> $std_out_file 
mv unfolded_band_structure.dat $results_file
if [ ! -s slurm-${SLURM_JOBID}.out ]
then

    rm  slurm-${SLURM_JOBID}.out

fi
