#!/usr/bin/env python
#  Copyright (C) 2013 Paulo V. C. Medeiros
#
#  This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
#
#  BandUP is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BandUP is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
import argparse
import numpy as np
import math
import sys
import os
import heapq

print ('                                                                                      \n'
       '===================================================================================== \n'
       '             BandUP: Band Unfolding code for Plane-wave based calculations            \n'
       '==================================================================================== \n'
       'Copyright (C) 2013, 2014 Paulo V. C. Medeiros                                        \n'
       '                         paume@ifm.liu.se                                            \n'
       '                         Computational Physics Division                              \n'
       '                         Department of Physics, Chemistry and Biology - IFM          \n'
       '                         Linkoping University                                        \n'
       '                         Sweden                                                      \n'
       'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                          \n'
       '==================================================================================== \n'
       '                                                                                     \n'
       '            Post-processing utility "find_band_centers_and_broadenings.py"           \n'
       '   >>> Obtaining Band Centers and Smearing Widths from the calculated delta_Ns <<<   \n'
       '        (Phys. Rev. B 89, 041407(R) (2014), Supplemental Material, Section II.C)     \n')

min_sum_dNs_for_a_band = 5.0e-2 # 1e-1 is OK
# I used to set threshold_dN_2b_trial_band_center = 1e-1 or 2e-1 (work for most of the cases), but this is a complicated parameter...
threshold_dN_2b_trial_band_center = 5.0e-2
prec_pos_band_centers = 0.01e-3 # in eV


parser = argparse.ArgumentParser()
parser.add_argument("input_file", nargs='?', default="unfolded_EBS_symmetry-averaged.dat", help="The name of the input file.")
parser.add_argument("-kpt", default="g", help="")
args = parser.parse_args()

# Reading input file
min_dN = 1e-5  # 1e-4 should be OK
input_file=args.input_file
args.kpt = args.kpt.lower()[0]

print "Reading input file (",input_file,")"
# It is faster to read the file manually instead of using np.loadtxt, because we can skip vaues that won't be used
kpts = []
energies = []
dNs = []
with open(input_file,'r') as file:
    for line in file:
        try:
           line_dN = float(line.split()[2])
           if(line_dN >= min_dN):
               kpts.append(float(line.split()[0])) 
               energies.append(float(line.split()[1])) 
               dNs.append(line_dN) 
        except:
            pass
dNs = np.array(dNs)
kpts = np.array(kpts)
energies = np.array(energies) 

# Eliminating duplicated points
sorted_index_data = np.lexsort((energies, kpts))
kpts = kpts[sorted_index_data]
energies = energies[sorted_index_data]
dNs = dNs[sorted_index_data]
for line_number in range(len(kpts)-1):
    if(kpts[line_number]==kpts[line_number+1] and
       energies[line_number]==energies[line_number+1]):
       dNs[line_number + 1] = max(dNs[line_number], dNs[line_number + 1])
       dNs[line_number] = 0.0
# Done with eliminating duplicated points
# Notice that all data is now sorted first by K and then by energy

n_points_on_file = len(energies)

high_sym_kpts_dict = {'k':0, 'g':1,'m':2}
high_sym_kpts = [0.0, 1.6979287434219599, 3.1683781686735077, 4.0173425410212928]
chosen_kpt_coord = high_sym_kpts[high_sym_kpts_dict[args.kpt]]

available_kpt_closest_to_the_chosen = min(kpts, key=lambda x:abs(x-chosen_kpt_coord))
chosen_kpt = available_kpt_closest_to_the_chosen
print "Parsing results for kpt coord = ", chosen_kpt," A**-1."

unique_kpts_coords = np.unique(kpts).tolist()
min_d_kpts = min([abs(unique_kpts_coords[ikpt+1] - unique_kpts_coords[-1]) for ikpt in range(len(unique_kpts_coords) - 1)])
dNs_for_current_kpt = [dNs[ip] for ip in range(n_points_on_file) if abs(kpts[ip] - chosen_kpt) <= min_d_kpts and dNs[ip] >= min_dN]
dNs_for_current_kpt = [dNs[ip] for ip in range(n_points_on_file) if abs(kpts[ip] - chosen_kpt) <= min_d_kpts]
energies_for_current_kpt = [energies[ip] for ip in range(n_points_on_file) if abs(kpts[ip] - chosen_kpt) <= min_d_kpts]


# Guessing positions of band centers
sorted_dN_indices = np.argsort(dNs_for_current_kpt)
indices_major_band_centers = [index for index in sorted_dN_indices if dNs_for_current_kpt[index] >= threshold_dN_2b_trial_band_center]
band_centers_1st_guesses = sorted([energies_for_current_kpt[index] for index in indices_major_band_centers])
# Done with guessing positions of band centers

# Defining a band as a set of energies spread around the guessed band centers in an interval of size max_bandwidth
enegies_spread_in_band = []
weights_of_enegies_spread_in_band = []

possible_energy_range = [(min(energies_for_current_kpt), 0.5 * (band_centers_1st_guesses[0] + band_centers_1st_guesses[1]))]
for iband in range(1,len(band_centers_1st_guesses)-1):
    possible_energy_range.append((0.5 * (band_centers_1st_guesses[iband] + band_centers_1st_guesses[iband-1]), 
                                  0.5 * (band_centers_1st_guesses[iband] + band_centers_1st_guesses[iband+1])))
possible_energy_range.append( (0.5 * (band_centers_1st_guesses[-1] + band_centers_1st_guesses[-2]), max(energies_for_current_kpt) + 1.0) )

for iband in range(len(band_centers_1st_guesses)):
    indices_of_enegies_spread_in_band = [index for index in range(len(energies_for_current_kpt)) 
                                         if energies_for_current_kpt[index] >= possible_energy_range[iband][0] and 
                                            energies_for_current_kpt[index] < possible_energy_range[iband][1]]
    enegies_spread_in_band.append([energies_for_current_kpt[index] for index in indices_of_enegies_spread_in_band])
    weights_of_enegies_spread_in_band.append([dNs_for_current_kpt[index] for index in indices_of_enegies_spread_in_band])


def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights) 
    return (average, math.sqrt(variance))

# SC cycle to refine the guesses of the positions of the band centers
min_sum_dNs_for_a_band = abs(min_sum_dNs_for_a_band)
threshold_dN_2b_trial_band_center = abs(threshold_dN_2b_trial_band_center)
if(min_sum_dNs_for_a_band < threshold_dN_2b_trial_band_center):
    min_sum_dNs_for_a_band = threshold_dN_2b_trial_band_center
    print "WARNING: Resetting min_sum_dNs_for_a_band because it is smaller than threshold_dN_2b_trial_band_center."

guess_band_centers = band_centers_1st_guesses
n_guesses_bc_start = len(guess_band_centers)
refined_band_centers = []
count = 0
converged = False
while(not converged):
    count += 1; 
    # Reducing the number of too close energy values
    for iband in range(len(guess_band_centers)):
        b_weight = np.sum(weights_of_enegies_spread_in_band[iband])
        bc, bwidth = weighted_avg_and_std(values=enegies_spread_in_band[iband], weights=weights_of_enegies_spread_in_band[iband])
        try:
            previous_bc, previous_bwidth = weighted_avg_and_std(values=enegies_spread_in_band[iband-1], weights=weights_of_enegies_spread_in_band[iband-1])
        except:
            previous_bc = min(energies_for_current_kpt)
            previous_bwidth = 0.0

        valid_bc = False
        if(b_weight < min_sum_dNs_for_a_band or abs(bc - previous_bc) < 2.0 * max([bwidth, previous_bwidth])): 
            try:
                previous_b_weight = np.sum(weights_of_enegies_spread_in_band[iband - 1])
                if(abs(b_weight / previous_b_weight) > 1.0):
                    del refined_band_centers[-1]
                    valid_bc = True
            except:
                pass
        else:
            valid_bc = True

        if(valid_bc):
            refined_band_centers.append(bc)

    # Since we have a new set of trial band centers, we have to update the interval of energies that are spread around them
    possible_energy_range = [(min(energies_for_current_kpt), 0.5 * (refined_band_centers[0] + refined_band_centers[1]))]
    for iband in range(1,len(refined_band_centers)-1):
        possible_energy_range.append((0.5 * (refined_band_centers[iband] + refined_band_centers[iband-1]), 
                                      0.5 * (refined_band_centers[iband] + refined_band_centers[iband+1])))
    possible_energy_range.append( (0.5 * (refined_band_centers[-1] + refined_band_centers[-2]), max(energies_for_current_kpt) + 1.0) )

    enegies_spread_in_band = []
    weights_of_enegies_spread_in_band = []
    for iband in range(len(refined_band_centers)):
        indices_of_enegies_spread_in_band = [index for index in range(len(energies_for_current_kpt)) 
                                             if energies_for_current_kpt[index] >= possible_energy_range[iband][0] and 
                                                energies_for_current_kpt[index] < possible_energy_range[iband][1]]
        enegies_spread_in_band.append([energies_for_current_kpt[index] for index in indices_of_enegies_spread_in_band])
        weights_of_enegies_spread_in_band.append([dNs_for_current_kpt[index] for index in indices_of_enegies_spread_in_band])

    if(all(abs(guess_band_centers[i] - refined_band_centers[i]) <= prec_pos_band_centers for i in range(len(refined_band_centers)))):
        converged = True
 
        band_centers = [refined_band_centers[iband] for iband in range(len(refined_band_centers)) 
                        if np.sum(weights_of_enegies_spread_in_band[iband]) >= min_sum_dNs_for_a_band]
        E_around_band_centers = [enegies_spread_in_band[iband] for iband in range(len(enegies_spread_in_band))
                                 if np.sum(weights_of_enegies_spread_in_band[iband]) >= min_sum_dNs_for_a_band]
        weights_of_E_around_band_centers = [weights_of_enegies_spread_in_band[iband] for iband in range(len(weights_of_enegies_spread_in_band))
                                            if np.sum(weights_of_enegies_spread_in_band[iband]) >= min_sum_dNs_for_a_band]

        n_guesses_bc_end = len(band_centers)
        if(count>0):
            print "Positions of the band centers converged to a precision of",1000.0 * prec_pos_band_centers,"meV after",count,"self-consistent steps."
        if(n_guesses_bc_end != n_guesses_bc_start):
            print "The # of band centers has been reduced from",n_guesses_bc_start,"(starting guess based on the values of delta_N) to",n_guesses_bc_end,"."
    else:
        guess_band_centers = refined_band_centers
        refined_band_centers = []


print "#Band center #Band width (meV) #Sum of dN"
for iband in range(len(band_centers)):
    b_weight = np.sum(weights_of_E_around_band_centers[iband])
    bc, bwidth = weighted_avg_and_std(values=enegies_spread_in_band[iband], weights=weights_of_enegies_spread_in_band[iband])
    print " %0.2f        %3i             %0.1f" % (bc, 1000.0 * bwidth, round(b_weight,1))


sys.exit(0)

