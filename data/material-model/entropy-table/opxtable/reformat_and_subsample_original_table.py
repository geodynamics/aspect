#!/bin/python

# This interpolation script can be used on a HeFESTo table in its original format,
# subsamples it, reorders the columns, and adds a header so that it can be used
# by the AsciiDataLookup class inside ASPECT.

import numpy as np

# Change this to the original input file
input_filename = "opxtable2_s.aspect"
output_filename = "material_table.txt"

data = np.genfromtxt(input_filename)

n_pressure_points_original = 1401
n_entropy_points_original = 241

table = data.reshape(n_pressure_points_original,n_entropy_points_original,9)

pressure_subsampling = 10
entropy_subsampling = 5

modified = table[::pressure_subsampling,:,:]
modified = modified[:,::entropy_subsampling,:]
modified = modified.reshape(modified.shape[0]*modified.shape[1],9)

# Change column order by providing permutation from old column index to new
# column index.
# Original order: T(K) P(bar) rho,kg/m3 alpha,1/K cp,J/K/kg vp,km/s vs,km/s h,J/kg s,J/kg/K
# Required order: s,J/kg/K P(bar) T(K) rho,kg/m3 alpha,1/K cp,J/K/kg vp,km/s vs,km/s h,J/kg
permutation = [2, 1, 3, 4, 5, 6, 7, 8, 0]
idx = np.empty_like(permutation)
idx[permutation] = np.arange(len(permutation))
modified = modified[:, idx]

headerfile = open("header_table.txt","r")
headertext = headerfile.read()
headertext = headertext[:-1]

np.savetxt(output_filename, modified,fmt='%.6g', header=headertext, comments='')
