#! /usr/bin/env python
import numpy as np
from os import getcwd
from os.path import join, realpath, dirname
import sys
# Combines the files from unfolding calculations using a rectangular BZ and gives an EBS in a hexagonal BZ

direcs=['K1-G', 'K2-G', 'G-M1', 'G-M2', 'M1-K1', 'M2-K1', 'M2-K2']
mult_factor = [2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0]
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

 
print 'Parsing G-M direction.'
data_G_M = [[],[],[]]
data_G_M[kpts_col] = data[2][kpts_col]
data_G_M[ener_col] = data[2][ener_col]
data_G_M[weight_col] = [0.0 for weight in data[2][weight_col]]
for i in range(len(data[2][weight_col])):
    if ((data[2][kpts_col][i] !=  data[3][kpts_col][i])or(data[2][ener_col][i] !=  data[3][ener_col][i])):
        sys.exit("Error message")
    data_G_M[weight_col][i] = (1.0*data[2][weight_col][i] + 2.0*data[3][weight_col][i])/3.0


print 'Parsing M-K direction.'
data_M_K = [[],[],[]]
data_M_K[kpts_col] = data[4][kpts_col]
data_M_K[ener_col] = data[4][ener_col]
data_M_K[weight_col] = [0.0 for weight in data[4][weight_col]]
for i in range(len(data[4][weight_col])):
    if ( (data[4][kpts_col][i] != data[5][kpts_col][i])or(data[4][kpts_col][i] != data[6][kpts_col][i])or(data[5][kpts_col][i] != data[6][kpts_col][i])or
         (data[4][ener_col][i] != data[5][ener_col][i])or(data[4][ener_col][i] != data[6][ener_col][i])or(data[5][ener_col][i] != data[6][ener_col][i]) ):
        sys.exit("Error message")
    data_M_K[weight_col][i] = (1.0*data[4][weight_col][i] + 1.0*data[5][weight_col][i] + 1.0*data[6][weight_col][i])/3.0

print 'Saving results.'
np.savetxt('K-G.dat',np.column_stack(data_K_G), fmt='%8.4f %8.4f %.3e')
np.savetxt('G-M.dat',np.column_stack(data_G_M), fmt='%8.4f %8.4f %.3e')
np.savetxt('M-K.dat',np.column_stack(data_M_K), fmt='%8.4f %8.4f %.3e')


print 'Saving results in a single file...'
total_data = [[],[],[]]
total_data[0] = [kpt_coord for kpt_coord in data_K_G[0]] 
for kpt_coord in  data_G_M[0]:
    total_data[0].append(kpt_coord) 
for kpt_coord in data_M_K[0]:
    total_data[0].append(kpt_coord) 
total_data[1] = [energy for energy in data_K_G[1]] 
for energy in  data_G_M[1]:
    total_data[1].append(energy) 
for energy in data_M_K[1]:
    total_data[1].append(energy) 
total_data[2] = [weight for weight in data_K_G[2]] 
for weight in  data_G_M[2]:
    total_data[2].append(weight) 
for weight in data_M_K[2]:
    total_data[2].append(weight) 
total_data = np.array(total_data)

sorted_index_data = np.lexsort((total_data[0], total_data[1]))
total_data[0] = total_data[0][sorted_index_data]
total_data[1] = total_data[1][sorted_index_data]
total_data[2] = total_data[2][sorted_index_data]
np.savetxt('EBS_K-G_G-M_M-K.dat',np.column_stack(total_data), fmt='%8.4f %8.4f %.3e')

