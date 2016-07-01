#!/usr/bin/python

# This script reads in a number of csv file handed over as an argument
# and plots them. In addition, it computes a fourier transform of
# the field given in the input (which is intended to be the 
# porosity) and the band angle distribution for the given input,
# and plots plots the mode of a lognormal distribution fitted to the 
# distribution of band angles for each input file.
# It also generates output files 'band_angle_uniform.csv' and 
# 'band_angle_adaptive.csv' that contain the number of degrees of 
# freedom for a given model, the standard deviation and the mean band
# angle, which can be used to plot the convergence of the band angle
# unsing the script 'plot_band_angle_dofs.py'. 

import numpy as np
import numpy.fft as fft
import scipy.interpolate
import matplotlib.pyplot as plt
import csv as csv
import math
import sys
from scipy import stats
from scipy.optimize import curve_fit

files_regular=[ "output_uniform/shear_bands_" + str(n) + ".csv" for n in [7,8,9,10]]
files_adaptive=[ "output_adaptive/shear_bands_"+str(n)+"_adaptive.csv" for n in [7,8,9,10]]

labels_regular=['    298 118 dofs','  1 186 054 dofs',' 4 731 398 dofs','18 899 974 dofs']
labels_adaptive=['    158 648 dofs','    693 612 dofs',' 3 069 422 dofs','10 781 758 dofs']

# customise the plot
fig, ax = plt.subplots(len(files_regular),2,figsize=(13.0, 7.0))
plt.setp(ax, xticks=np.linspace(0.0,0.004,5),xticklabels=['0', '1', '2', '3', '4'])
plt.setp(ax, yticks=np.linspace(0.000005,0.000995,3),yticklabels=['0', '0.5', '1'])

plt.subplots_adjust(left=0.011, right=0.95, top=0.965, bottom=0.07, hspace = 0.3, wspace=0.0)

ax[0,0].set_title('Uniform mesh')
ax[0,1].set_title('Adaptive mesh')
ax[len(files_regular)-1,0].set_xlabel('x in mm')
ax[len(files_regular)-1,1].set_xlabel('x in mm')

plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
plt.setp([a.get_xticklabels() for a in ax[1, :]], visible=False)
plt.setp([a.get_xticklabels() for a in ax[2, :]], visible=False)
plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)

fig_hist, ax_hist = plt.subplots(2,len(files_regular), sharex='col', sharey='row', figsize=(15.0, 6.0))

ax_hist[0,0].set_title('Uniform mesh')
ax_hist[1,0].set_title('Adaptive mesh')
for n in range(0,len(files_regular)):
	ax_hist[1,n].set_xlabel('Band angle in degrees')

# write output to a file
f_uniform = open('band_angle_uniform.csv', 'w')
f_adaptive = open('band_angle_adaptive.csv', 'w')
f_uniform.write("dofs   standard_dev   mean \n")
f_adaptive.write("dofs   standard_dev   mean \n")

for n in range(0,len(files_regular)):

# read in data to plot
	data = []
	data_mesh = []

	data.append(np.genfromtxt(files_adaptive[n], dtype = float, names = True))
	data_mesh.append(np.genfromtxt(files_regular[n], dtype = float, names = True))

	coordsX_regular=[]
	coordsY_regular=[]
	porosity_regular=[]

	for i in range(0,len(data_mesh)):
		coordsX_regular=np.append(coordsX_regular, data_mesh[i]['x'])
		coordsY_regular=np.append(coordsY_regular, data_mesh[i]['y'])
		porosity_regular=np.append(porosity_regular, data_mesh[i]['porosity'])

	x_regular = np.asarray(coordsX_regular)
	y_regular = np.asarray(coordsY_regular)
	p_regular_plot = np.asarray(porosity_regular)
	p_regular = p_regular_plot - p_regular_plot.mean()

	coordsX_adaptive=[]
	coordsY_adaptive=[]
	porosity_adaptive=[]

	for i in range(0,len(data)):
		coordsX_adaptive=np.append(coordsX_adaptive, data[i]['x'])
		coordsY_adaptive=np.append(coordsY_adaptive, data[i]['y'])
		porosity_adaptive=np.append(porosity_adaptive, data[i]['porosity'])

	x_adaptive = np.asarray(coordsX_adaptive)
	y_adaptive = np.asarray(coordsY_adaptive)
	p_adaptive_plot = np.asarray(porosity_adaptive)
	p_adaptive = p_adaptive_plot - p_adaptive_plot.mean()

	xmin = x_regular.min()
	xmax = x_regular.max()
	ymin = y_regular.min()
	ymax = y_regular.max()
	pmin = min(p_regular_plot.min(),p_adaptive_plot.min())
	pmax = max(p_regular_plot.max(),p_adaptive_plot.max())

	dx = 1e300
	sortedx = np.sort(x_regular)
	for i in range(1,len(sortedx)):
		d = sortedx[i]-sortedx[i-1]
		if d>1e-10:
			dx = min(dx, d)
	dy = 1e300
	sortedy = np.sort(y_regular)
	for i in range(1,len(sortedy)):
		d = sortedy[i]-sortedy[i-1]
		if d>1e-10:
			dy = min(dy, d)

	nrows = np.round((ymax - ymin) / dy)+1
	ncols = np.round((xmax - xmin) / dx)+1
	aspect_ratio = nrows/ncols

	# Then we make an empty 2D grid...
	grid = np.zeros((nrows, ncols), dtype=np.float)
	xg = np.arange(xmin, xmax, dx)
	yg = np.arange(ymin, ymax, dy)

	xy_regular=np.transpose(np.array([x_regular,y_regular]))
	xy_adaptive=np.transpose(np.array([x_adaptive,y_adaptive]))
	xyg = np.meshgrid(xg, yg)
	grid_regular_plot = scipy.interpolate.griddata(xy_regular, p_regular_plot, tuple(xyg), method='nearest')
	grid_adaptive_plot = scipy.interpolate.griddata(xy_adaptive, p_adaptive_plot, tuple(xyg), method='nearest')
	grid_regular = scipy.interpolate.griddata(xy_regular, p_regular, tuple(xyg), method='nearest')
	grid_adaptive = scipy.interpolate.griddata(xy_adaptive, p_adaptive, tuple(xyg), method='nearest')

	# apply a hamming window function to make porosity periodic
	hx = np.ones(len(xg))
	hy = np.hanning(len(yg))
	ham2d = np.outer(hx,hy)
	ham2d_trans = np.transpose(ham2d)

	grid_window_regular = np.multiply(grid_regular,ham2d_trans)
	grid_window_adaptive = np.multiply(grid_adaptive,ham2d_trans)

	# make the fourier transformation
	fourier_regular = fft.fft2(grid_window_regular)
	fourier_centered_regular = fft.fftshift(fourier_regular)
	fourier_adaptive = fft.fft2(grid_window_adaptive)
	fourier_centered_adaptive = fft.fftshift(fourier_adaptive)

	fouriermax_regular = np.absolute(fourier_centered_regular).max()
	fourierabs_regular = np.absolute(fourier_centered_regular)/fouriermax_regular
	fouriermax_adaptive = np.absolute(fourier_centered_adaptive).max()
	fourierabs_adaptive = np.absolute(fourier_centered_adaptive)/fouriermax_adaptive

	# interpolate fourier transform to spherical grid
	# this allows us to have an equal number of points per bin 
	fourierx=np.linspace(-0.5*len(xg),0.5*len(xg),ncols)
	fouriery=np.linspace(-0.5*len(xg),0.5*len(xg),nrows)

	(fx,fy) = np.meshgrid(fourierx,fouriery)
	fourierGrid = np.transpose(np.array([fx.flatten(),fy.flatten()]))

	phi=np.linspace(0,2*math.pi,num=360,endpoint=False)
	rad=np.linspace(1e-5,0.1*len(xg),300)

	sphericalGrid = np.empty([len(phi)*len(rad),2])

	sphericalGrid[:,0] = (np.outer(rad,np.cos(phi))).flatten()
	sphericalGrid[:,1] = (np.outer(rad,np.sin(phi))).flatten()

	intensity_regular=scipy.interpolate.griddata(fourierGrid, fourierabs_regular.flatten(), sphericalGrid, method='nearest')
	intensity_adaptive=scipy.interpolate.griddata(fourierGrid, fourierabs_adaptive.flatten(), sphericalGrid, method='nearest')

	# calculate the angle for each point
	angle=[]

	for i in range(0,int(len(sphericalGrid[:,0]))):
		a = math.atan2(sphericalGrid[i,1],sphericalGrid[i,0])
		if (a<=0):
			a+=np.pi
		angle.append(90.0 - a/np.pi*180.0)

	# input data
	im0=ax[n,0].imshow(grid_regular_plot, interpolation='nearest', 
	        extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='RdBu_r')
	im1=ax[n,1].imshow(grid_adaptive_plot, interpolation='nearest', 
	        extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='RdBu_r')
	ax[n,0].text(0.73, 1.15, labels_regular[n], transform=ax[n,0].transAxes, va='top')
	ax[n,1].text(0.73, 1.15, labels_adaptive[n], transform=ax[n,1].transAxes, va='top')
	ax[n,0].set_ylabel('y in mm')

	cbar_ax = fig.add_axes([0.455, 0.7815-n*0.237, 0.015, 0.183])
	fig.colorbar(im0, cax=cbar_ax)
	cbar_ax = fig.add_axes([0.925, 0.7815-n*0.237, 0.015, 0.183])
	fig.colorbar(im1, cax=cbar_ax)

	# band angle histogram
	bins=np.linspace(1e-10,90,18)
	ax_hist[0,n].hist(angle, bins, normed=True, weights=intensity_regular, rwidth=0.7)
	ax_hist[1,n].hist(angle, bins, normed=True, weights=intensity_adaptive, rwidth=0.7)

	# fit lognormal
	bins_fit=np.linspace(1e-10,90,180)
	y_hist_regular, bin_edges_regular = np.histogram(angle, bins_fit, normed=True, weights=intensity_regular)
	y_hist_adaptive, bin_edges_adaptive = np.histogram(angle, bins_fit, normed=True, weights=intensity_adaptive)
	x_hist=bins_fit[1:]

	(shape_out, scale_out), pcov = curve_fit(lambda xdata, shape, scale: stats.lognorm.pdf(xdata, shape, loc=0, scale=scale), x_hist, y_hist_regular, p0=[0.25, 20])
	ax_hist[0,n].plot(x_hist,stats.lognorm.pdf(x_hist, shape_out, loc=0, scale=scale_out), 'r-', lw=3, label='Fitted distribution')
	ax_hist[0,n].text(0.44, 0.94, labels_regular[n], transform=ax_hist[0,n].transAxes, va='top')
	f_uniform.write(labels_regular[n].replace(" ", "").replace("dofs", "") + " " + str(shape_out) + " " + str(scale_out) + "\n")

	(shape_out, scale_out), pcov = curve_fit(lambda xdata, shape, scale: stats.lognorm.pdf(xdata, shape, loc=0, scale=scale), x_hist, y_hist_adaptive, p0=[0.25, 20])
	ax_hist[1,n].plot(x_hist,stats.lognorm.pdf(x_hist, shape_out, loc=0, scale=scale_out), 'r-', lw=3, label='Fitted distribution')
	ax_hist[1,n].text(0.44, 0.94, labels_adaptive[n], transform=ax_hist[1,n].transAxes, va='top')
	f_adaptive.write(labels_adaptive[n].replace(" ", "").replace("dofs", "") + " " + str(shape_out) + " " + str(scale_out) + "\n")

f_uniform.close()
f_adaptive.close()
fig.savefig('fft', dpi=300)
fig_hist.savefig('fft_hist', dpi=200)


