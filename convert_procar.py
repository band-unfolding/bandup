#! /usr/bin/env python
from bandupy.vasp import procar2bandup

picked_orbitals = 'sp3'
selected_ion_indices = range(32) # Only C
#selected_ion_indices = range(32, 32+64) # Only Pt
procar2bandup(picked_orbitals=picked_orbitals, ignore_off_diag=False,
              force_hermitian=True,
              selected_ion_indices=selected_ion_indices, use_dual=True)
