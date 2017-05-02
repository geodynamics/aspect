#!/usr/bin/python

# Given one or several input files that contain the number of degrees
# of freedom of a model, the standard deviation and the mean of the 
# band angle distribution, this script outputs a convergence plot for
# the band angle.  

import numpy as np
import matplotlib.pyplot as plt
import csv as csv
import sys

filename=[]
if len(sys.argv)<2:
	print "usage: band_angle_uniform.csv"
	filename=["band_angle_uniform.csv"]
else:
	for i in range (1,len(sys.argv)):
		filename.append(sys.argv[i])

nplots=2

data = []
for i in range (0,len(filename)):
	data.append(np.genfromtxt(filename[i], dtype = float, names = True))

resolution=[]
width=[]
angle=[]

# data to plot
for i in range(0,len(data)):

#initialize empty lists we fill first and then append to the big list
	res_1_element = []
	wid_1_element = []
	ang_1_element = []

	for j in range(0,len(data[i])):
		res_1_element=np.append(res_1_element, data[i][j]['dofs'])
		wid_1_element=np.append(wid_1_element, data[i][j]['standard_dev'])
		ang_1_element=np.append(ang_1_element, data[i][j]['mean'])

	resolution.append(res_1_element)
	width.append(wid_1_element)
	angle.append(ang_1_element)

print resolution

# Two subplots, the axes array is 1-d
fig, ax = plt.subplots(1)
fig.set_size_inches(4,2.5)
colors=['k','r']
labels=['uniform mesh','adaptive mesh']

for i in range(0,len(data)):
	ax.plot(resolution[i],np.exp(np.log(angle[i])-width[i]*width[i]),marker='x',linestyle='None',color=colors[i],label=labels[i],markersize=10,mew=2)
ax.set_ylim([17.5,19.5])

# add labels and legend
ax.set_xlabel('Number of degrees of freedom')
ax.set_xscale('log')
ax.set_ylabel('Mode of the lognormal')

fig.savefig('band_angle', bbox_inches='tight',dpi=200)

