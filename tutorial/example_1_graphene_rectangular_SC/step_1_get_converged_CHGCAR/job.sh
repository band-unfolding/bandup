#!/bin/bash
#SBATCH --nodes 5
#SBATCH -A liu5 -p green
#SBATCH -J  GE0
#SBATCH -t 1-00:00:00

# Path to VASP EXE
vasp=$vasp_half

vasp_files='vasp_input_to_get_CHGCAR'

cp $vasp_files/POSCAR_graphene_rect_zigzag_2x3 POSCAR
cp $vasp_files/POTCAR .
cp $vasp_files/KPOINTS .
cp $vasp_files/INCAR  .

ulimit -s unlimited
export OMP_NUM_THREADS=1
mpprun  $vasp  > vasp.out

rm -f WAVECAR CONTCAR DOSCAR EIGENVAL INCAR KPOINTS POTCAR PROCAR vasprun.xml XDATCAR OSZICAR IBZKPT EXHCAR CHG TMPCAR PCDAT ELFCAR PROOUT

#End of script
