#! /usr/bin/env python
from bandupy.vasp import procar2bandup

#selected_ion_indices = xrange(32) # Only C
#selected_ion_indices = range(32, 32+64) # Only Pt
#selected_ion_indices = range(32, 32+8) # Only O
procar2bandup()
