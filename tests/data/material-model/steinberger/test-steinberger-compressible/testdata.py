#!/bin/python

import numpy as np
import csv

def calc_adiabatic_conditions():
	npoints = 1500
	density = 3340
	alpha = 3e-5
	cp = 1200
	gravity = 9.81
	compressibility = 0.25e-11
	max_depth = (6371000 - 3481000) * npoints/1000
	delta_depth = max_depth / (npoints-1)
	temps = np.zeros(npoints)
	press = np.zeros(npoints)

	temps[0] = 1600
	press[0] = 0

	for i in range(npoints - 1):
		press[i+1] = press[i] + density * np.exp(compressibility*press[i]) * gravity * delta_depth
		temps[i+1] = temps[i] * (1 + alpha*gravity*delta_depth/cp)
	return temps, press

def get_Ta(temps,press,pa):
	for i in range(len(temps)-1):
		if press[i+1] > pa:
			x = (pa-press[i]) / (press[i+1]-press[i])
			return temps[i] + x * (temps[i+1]-temps[i])
	return -1


aditemps,adipress = calc_adiabatic_conditions()

# Number of temperature bands
ntemp = 101
# Number of depth bands
npress = 101

# Generate date for testdata.txt
data = np.zeros((ntemp*npress,8))

temps = np.linspace(250,4250,ntemp)
press = np.linspace(0,2e6,npress)

meshtemps,meshpress = np.meshgrid(temps,press)

data[:,0] = np.ravel(meshtemps)
data[:,1] = np.ravel(meshpress)

data[:,3] = np.ones((ntemp*npress)) * 3e-5

for i in range(data[:,2].size):
	adiabatic_temperature = get_Ta(aditemps,adipress,data[i,1]*1e5)
	data[i,2] = 3340.0 * (1 - (data[i,0]-adiabatic_temperature) * 3e-5) \
                       * np.exp(data[i,1] * 0.25e-6)


data[:,4] = 1200
data[:,5] = -1
data[:,6] = -1
data[:,7] = 0

np.savetxt("data.txt",data,fmt='%.12lg')

# Generate viscosity prefactor
ndepth = 21
vis_prefact = np.zeros((ndepth,2))
vis_prefact[:, 1] = np.linspace(0, 2889, ndepth)
np.savetxt("test-viscosity-prefactor.txt", vis_prefact, fmt='%.6lg')

# Generate over-resolved viscosity prefactor
ndepth = 101
vis_prefact = np.zeros((ndepth,2))
vis_prefact[:, 1] = np.linspace(0, 2889, ndepth)
np.savetxt("test-viscosity-prefactor-overres.txt", vis_prefact, fmt='%.6lg')
