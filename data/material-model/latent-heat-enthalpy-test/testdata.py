#!/bin/python

# This script creates a property table to recreate the latent heat benchmark
# with a self-consistent formulation for cp and alpha derived from an 
# enthalpy table in pressure-temperature space. Unfortunately the 
# benchmark itself assumes some thermodynamically inconsistent relations, so
# we need to use our knowledge about the expected result to set up the table

import numpy as np
import csv

gravity = 10.0
density = 3400.0
density_jump = 115.6
density_after_transition = density + density_jump
cp = 1000
alpha = 0.0
depth = 500000
zero_pressure = density*gravity*depth
zero_temperature = 1000
clapeyron_slope = 1e7

# note that the analytical latent heat model interprets the width to be
# a half_width, therefore 40000 m here is equivalent to the 20000 m of the
# cookbook
transition_width = 40000 
transition_width_pressure = transition_width * gravity * density
dS = -96.71

# This temperature change is usually the value to compare the result against, however
# since the enthalpy change is calculated as temperature times entropy change
# we need to know the temperature of the finished phase transition to calculate
# the enthalpy change. This is the only place we use this result.
dT = 109.08
TdS = dS * (zero_temperature+dT)

def phase_transition(temperature,pressure):
	equilibrium_transition_pressure = zero_pressure + (temperature-zero_temperature) * clapeyron_slope
        tanh_phase = 0.5 * (1 + np.tanh((pressure * 1e5 - equilibrium_transition_pressure) / transition_width_pressure))
	return tanh_phase

ntemp = 21
npress = 51

data = np.zeros((ntemp*npress,8))

temps = np.linspace(850,1250,ntemp)
press = np.linspace(0,4e5,npress)

meshtemps,meshpress = np.meshgrid(temps,press)

data[:,0] = np.ravel(meshtemps) # temperature
data[:,1] = np.ravel(meshpress) # pressure
for i in range(data[:,2].size): # density
	phase_proportion = phase_transition(data[i,0],data[i,1])
	data[i,2] = density + phase_proportion * density_jump
data[:,3] = 0 # thermal expansivity
data[:,4] = cp # specific heat
data[:,5] = -1 # vp
data[:,6] = -1 # vs

for i in range(data[:,2].size): # enthalpy
	phase_proportion = phase_transition(data[i,0],data[i,1])
	if phase_proportion == 0.0:
		data[i,7] = cp * data[i,0] + (1. / density_after_transition) * data[i,1] * 1e5
	elif phase_proportion == 1.0:
		data[i,7] = TdS + cp * data[i,0] + (1. / density_after_transition) * data[i,1] * 1e5
	else:
		data[i,7] = TdS * phase_proportion + cp * data[i,0] + (1. / density_after_transition) * data[i,1] * 1e5

np.savetxt("data.txt",data,fmt='%.12lg')
