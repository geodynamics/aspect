#!/usr/bin/python

# This script can be used to plot data for the magmatic shear bands testcase
# for Newtonian rheology and analyze it using the relations given in Spiegelman (2003):
# Linear analysis of melt band formation by simple shear, Geochemistry, Geophysics,
# Geosystems, 4(9), 8615.

# It uses the output of the script 'run_angle_melt_bands.sh' and plots the analytically 
# predicted growth rate and the numerically computed growth rate over the band angle.

import numpy as np
import matplotlib.pyplot as plt
import colors
import math

figsize=(7,5)
prop={'size':12}
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
figure=plt.figure(dpi=100,figsize=figsize)

# ndofs,u,p,p_f,p_c,phi,u_f
file_name="plane_wave_melt_bands_angle"

data = []
data.append(np.genfromtxt(file_name,delimiter=' ', dtype = float))

angle=[]
analytical=[]
numerical=[]

end=len(data[0])

# data to plot
for j in range(0,end):
	angle.append(data[0][j][0]*180/math.pi)
	analytical.append(data[0][j][1])
	numerical.append(data[0][j][2])

plt.plot(angle[0:end],numerical[0:end]," ",color=colors.color(3), marker=colors.marker(1), label='numerical')
plt.plot(angle[0:end],analytical[0:end],"--", color="black", marker="x", mew=1.5, ms=4, label='analytical')

plt.xlim([-5, 185])
#plt.ylim([4e-3,0.4])
plt.xlabel("Initial band angle in $^{\circ}$")
plt.ylabel("Melt band growth rate $\dot s$ in s$^{-1}$")
plt.xticks([0,45,90,135,180])
plt.grid(True)
 

plt.legend(loc = "lower left",prop=prop)
plt.savefig('growth_rate_angle.pdf', #bbox_extra_artists=(legend,), 
            bbox_inches='tight',dpi=200)

