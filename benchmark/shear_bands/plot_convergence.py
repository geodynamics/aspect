#!/usr/bin/python

# This script can be used to plot data for the magmatic shear bands testcase
# for Newtonian rheology and analyze it using the relations given in Spiegelman (2003):
# Linear analysis of melt band formation by simple shear, Geochemistry, Geophysics,
# Geosystems, 4(9), 8615.

# It uses the output of the script 'run_plane_melt_bands.sh' and plots the growth rate
# error over the number of degrees of freedom in the model.


import numpy as np
import matplotlib.pyplot as plt
import colors

figsize=(7,5)
prop={'size':12}
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
figure=plt.figure(dpi=100,figsize=figsize)

# ndofs,u,p,p_f,p_c,phi,u_f
file_name=["plane_wave_melt_bands_8pi","plane_wave_melt_bands_16pi"]

data = []
for j in range(0,len(file_name)):
	data.append(np.genfromtxt(file_name[j],delimiter=' ', dtype = float))

dofs=[]
error=[]
resolution=[]

# data to plot
for i in range(0,len(file_name)):
	for j in range(0,6):
		dofs.append(abs(data[i][j][0]))
		error.append(abs(data[i][j][1]))
		resolution.append(3000./abs(data[i][j][0]))

plt.loglog(dofs[0:6],error[0:6],color=colors.color(1), marker=colors.marker(1),label='$k=8\pi \cdot 10^3$ m$^{-1}$')
plt.loglog(dofs[7:12],error[7:12],color=colors.color(3), marker=colors.marker(1),label='$k=16\pi \cdot 10^3$ m$^{-1}$')
plt.loglog(dofs[0:6],resolution[0:6],"--",color="black", label='$\mathcal{O}(h^2)$')

plt.xlim([2e7, 3000])
#plt.ylim([4e-3,0.4])
plt.xlabel("Degrees of freedom")
plt.ylabel("Relative error in $\dot s$")
plt.grid(True)
 

plt.legend(loc = "upper left",prop=prop)
plt.savefig('growth_rate_error.pdf', #bbox_extra_artists=(legend,), 
            bbox_inches='tight',dpi=200)

