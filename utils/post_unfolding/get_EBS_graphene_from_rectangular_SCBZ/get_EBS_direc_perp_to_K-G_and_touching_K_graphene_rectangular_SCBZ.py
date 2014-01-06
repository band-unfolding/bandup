#! /usr/bin/env python
import numpy as np
from os import getcwd
from os.path import join, realpath, dirname
import sys
# Combines the files from unfolding calculations using a rectangular BZ and gives an EBS in a hexagonal BZ

direcs=['perp_to_K1-G_and_touching_K1', 'perp_to_K2-G_and_touching_K2']
mult_factor = [2.0, 1.0]
# Reading data from input file
mypath = realpath(join(getcwd(), dirname(__file__)))

data = []
for direc in direcs:
    myfile = 'unfolded_band_structure_'+direc+'.dat'
    print 'Reading file ', myfile 
    data_direc = np.loadtxt(join(mypath,myfile),
                 usecols=(0,1,2),
                 unpack=True)
    data.append(data_direc)

kpts_col = 0
ener_col = 1
weight_col = 2

print 'Parsing K-G direction.'
data_K_G = [[],[],[]]
data_K_G[kpts_col] = data[0][kpts_col]
data_K_G[ener_col] = data[0][ener_col]
data_K_G[weight_col] = [0.0 for weight in data[0][weight_col]]
for i in range(len(data[0][weight_col])):
    if ((data[0][kpts_col][i] !=  data[1][kpts_col][i])or(data[0][ener_col][i] !=  data[1][ener_col][i])):
        sys.exit("Error message")
    data_K_G[weight_col][i] = (2.0*data[0][weight_col][i] + 1.0*data[1][weight_col][i])/3.0


print 'Saving results.'
np.savetxt('EBS_perp_to_K-G_and_touching_K.dat',np.column_stack(data_K_G), fmt='%8.4f %8.4f %.3e')

