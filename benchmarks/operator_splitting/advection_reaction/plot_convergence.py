#!/usr/bin/python

# This script can be used to plot errors for the operator
# splitting 'advection reaction' benchmark. 


import numpy as np
import matplotlib.pyplot as plt

figsize=(7,5)
prop={'size':12}
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']='\\usepackage{relsize}'
plt.rc('font', family='sanserif')
figure=plt.figure(dpi=100,figsize=figsize)

data = np.genfromtxt("output_reaction_advection",delimiter=' ', dtype = float)

dofs=[]
error=[]
errorT=[]
resolution=[]

# data to plot
for i in range(0,len(data)):
	dofs.append(abs(data[i][0]))
	error.append(abs(data[i][2]))
	errorT.append(abs(data[i][3]))
	resolution.append((data[0][2]-4e-4)*abs(data[i][0])/abs(data[0][0]))

fig = plt.figure()
#plt.grid(True)
plt.gca().invert_xaxis()
ax1 = fig.add_subplot(111)
ax1.set_yscale('log')
ax2 = ax1.twiny()

ax1.set_xscale('log')
ax1.set_xlim([2,1e-3])
ax1.plot(dofs[10:17],error[10:17],color='#ff0000', marker='D', markersize=8, label='Advection time step $=0.25=$ const.')
ax1.plot(dofs[10:17],errorT[10:17],color='#8b0000', marker='d', label='Temperature')
ax1.plot(0,0,color='#4169e1', marker='D', markersize=8, label='Reaction time step $=2\cdot 10^{-3}=$ const.')
ax1.plot(0,0,color='#00008b', marker='d',label='Temperature')
ax1.plot(0,0,"--",color="#00008b", label='$\mathcal{O}(\Delta t)$')

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_xlabel('Reaction time step $\Delta t_R$', color='r')
ax1.tick_params('x', colors='r')

def tick_function(X):
    V = X*1
    return ["%.3f" % z for z in V]

ax2.set_xlim([2,1e-3])
ax2.set_xticks([1, 0.1, 0.01])
ax2.set_xticklabels(tick_function([1, 0.1, 0.01]))
ax2.set_xscale('log')

ax2.plot(dofs[0:10],error[0:10],color='#4169e1', marker='D', markersize=8)
ax2.plot(dofs[0:10],errorT[0:10],color='#00008b', marker='d')
ax2.plot(dofs[0:6],resolution[0:6],"--",color="#00008b")

ax2.set_xlabel('Advection time step $\Delta t_A$', color='b')
ax2.tick_params('x', colors='b')

ax1.legend(loc='upper right')

ax1.set_ylabel("Relative error in $C$ and $T$")


plt.legend(loc = "upper left",prop=prop)
plt.savefig('composition_error.pdf', #bbox_extra_artists=(legend,), 
            bbox_inches='tight',dpi=200)

