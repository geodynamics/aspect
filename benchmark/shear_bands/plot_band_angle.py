#!/usr/bin/python

# this script can be used to generate plots from csv files

import numpy as np
import scipy.fftpack as fft
import scipy.interpolate
import matplotlib.pyplot as plt
import csv as csv
import math
import sys

if len(sys.argv)!=2:
	print "usage: shear_bands.csv"
	sys.exit(0)
filename = sys.argv[1]

nplots=3

#resolution=[2000]#[128000]
#label=['refinement 2','refinement 3','refinement 4','analytical']
#label=['resolution 2.5 m','resolution 1.25 m','resolution 0.63 m','analytical solution']

data = []

data.append(np.genfromtxt(filename, dtype = float, names = True))

coordsX=[]
coordsY=[]
porosity=[]

# data to plot
for i in range(0,len(data)):
	coordsX=np.append(coordsX, data[i]['x'])
	coordsY=np.append(coordsY, data[i]['y'])
	porosity=np.append(porosity, data[i]['porosity'])

x = np.asarray(coordsX)
y = np.asarray(coordsY)
p = np.asarray(porosity) - 0.05

print len(x), len(y), len(p)

xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

dx = 1e300
sortedx = np.sort(x)
for i in range(1,len(sortedx)):
	d = sortedx[i]-sortedx[i-1]
	if d>1e-10:
		dx = min(dx, d)
dy = 1e300
sortedy = np.sort(y)
for i in range(1,len(sortedy)):
	d = sortedy[i]-sortedy[i-1]
	if d>1e-10:
		dy = min(dy, d)

print "dx=",dx,"dy=",dy

nrows = np.round((ymax - ymin) / dy)+1
ncols = np.round((xmax - xmin) / dx)+1
aspect_ratio = nrows/ncols

# First we convert our x and y positions to indicies...
idx = np.round((x - xmin) / dx).astype(np.int)
idy = np.round((ymax - y) / dy).astype(np.int)

# Then we make an empty 2D grid...
grid = np.zeros((nrows, ncols), dtype=np.float)
xg = np.arange(xmin, xmax, dx)
yg = np.arange(ymin, ymax, dy)


#f = scipy.interpolate.griddata(x, y, p, kind='linear')
#grid = f(xg, yg)

xy=np.transpose(np.array([x,y]))
xyg = np.meshgrid(xg, yg)
grid = scipy.interpolate.griddata(xy, p, xyg, method='linear')

# make the fourier transformation
fourier = fft.fft2(grid)
fourier_centered = fft.fftshift(fourier)

fig, ax = plt.subplots(nplots)
ax[0] = plt.subplot2grid((2,2), (0, 0), colspan=2)
ax[1] = plt.subplot2grid((2,2), (1, 0))
ax[2] = plt.subplot2grid((2,2), (1, 1))

fouriermax = np.absolute(fourier_centered).max()
fourierabs = np.absolute(fourier_centered)/fouriermax

intensity_integral=[]
intensity_points=[]
angle=np.copy(grid)

# calculate the angle for each point
for i in range(0,int(ncols)):
	for j in range(0,int(nrows)):
		a = -math.atan2((i-(ncols-1)/2)*aspect_ratio, (j-(nrows-1)/2))
		if (a<0):
			a+=np.pi
		angle[j][i] = a/np.pi*180.0

# And now we plot it:
# input data
im0 = ax[0].imshow(grid, interpolation='nearest', 
        extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='RdBu_r')
# fourier transporm
im1 = ax[1].imshow((np.log(np.absolute(fourier_centered)/fouriermax)), interpolation='nearest', 
        extent=(xmin, xmax, ymax, ymin), aspect=1.0/aspect_ratio)
#plt.colorbar(im0)

# band angle histogram
angle_flat = angle.flatten()
fourier_flat = fourierabs.flatten()
bins=np.linspace(0,90,17)
ax[2].hist(angle_flat, bins, normed=True, weights=fourier_flat, rwidth=0.7)

fig.savefig('fft', dpi=200)


