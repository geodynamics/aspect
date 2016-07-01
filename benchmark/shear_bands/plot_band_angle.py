#!/usr/bin/python

# This script reads in a single csv file handed over as an argument
# and plots it. In addition, it computes a fourier transform of
# the field given in the input (which is intended to be the 
# porosity) and the band angle distribution for the given input,
# and plots both.

import numpy as np
import numpy.fft as fft
import scipy.interpolate
import matplotlib.pyplot as plt
import csv as csv
import math
import sys
from scipy import stats
from scipy.optimize import curve_fit

if len(sys.argv)!=2:
	print "usage: shear_bands.csv"
	sys.exit(0)
filename = sys.argv[1]

nplots=3

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
p = np.asarray(porosity)
p = p - p.mean()

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

# Then we make an empty 2D grid...
grid = np.zeros((nrows, ncols), dtype=np.float)
xg = np.arange(xmin, xmax, dx)
yg = np.arange(ymin, ymax, dy)

xy=np.transpose(np.array([x,y]))
xyg = np.meshgrid(xg, yg)
grid = scipy.interpolate.griddata(xy, p, tuple(xyg), method='nearest')

# apply a hamming window function to make porosity periodic
hx = np.ones(len(xg))
hy = np.hanning(len(yg))
ham2d = np.outer(hx,hy)
ham2d_trans = np.transpose(ham2d)
grid_window = np.multiply(grid,ham2d_trans)

# make the fourier transformation
fourier = fft.fft2(grid_window)
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
		a = math.atan2((j-(nrows)/2),(i-(ncols)/2)*aspect_ratio)
		if (a<0):
			a+=np.pi
		angle[j][i] = 90.0 - a/np.pi*180.0

# And now we plot it:
# input data
im0 = ax[0].imshow(grid_window, interpolation='nearest', 
        extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='RdBu_r')
# fourier transporm
fourier2plot = fourierabs[np.round(len(yg)*7/16):np.round(len(yg)*9/16),np.round(len(xg)*7/16):np.round(len(xg)*9/16)]
im1 = ax[1].imshow(fourier2plot, interpolation='nearest',
        extent=(xmin, xmax, ymin, ymax), origin='lower', aspect=1.0/aspect_ratio)
ax[1].axis('off')
#plt.colorbar(im0)

# band angle histogram
angle_flat = angle.flatten()
fourier_flat = fourierabs.flatten()
bins=np.linspace(0.1,90,18)
ax[2].hist(angle_flat, bins, normed=True, weights=fourier_flat, rwidth=0.7)

# fit lognormal
bins_fit=np.linspace(0.1,90,100)
y_hist, bin_edges = np.histogram(angle_flat, bins_fit, normed=True, weights=fourier_flat)
x_hist=bins_fit[1:]

(shape_out, scale_out), pcov = curve_fit(lambda xdata, shape, scale: stats.lognorm.pdf(xdata, shape, loc=0, scale=scale), x_hist, y_hist, p0=[0.25, 20])

ax[2].plot(x_hist,stats.lognorm.pdf(x_hist, shape_out, loc=0, scale=scale_out), 'r-', lw=3, label='Fitted distribution')
print shape_out, scale_out

fig.savefig('fft', dpi=200)


