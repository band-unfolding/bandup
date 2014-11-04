#!/bin/bash
#SBATCH --nodes 24
#SBATCH -A liu5 -p green
#SBATCH -J  BulkSiE0
#SBATCH -t 1-00:00:00

# Path to VASP EXE
vasp=$vasp_half

vasp_files='vasp_input_to_get_CHGCAR'

cp $vasp_files/POSCAR_Si_SC POSCAR
cp $vasp_files/POTCAR .
cp $vasp_files/KPOINTS .
cp $vasp_files/INCAR  .

mpprun  $vasp  > vasp.out

rm -f WAVECAR CONTCAR DOSCAR EIGENVAL INCAR KPOINTS POTCAR vasprun.xml XDATCAR OSZICAR IBZKPT EXHCAR CHG TMPCAR PCDAT ELFCAR PROOUT

#End of script
