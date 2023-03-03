#!/usr/bin/env python
# coding: utf-8

# This script uses the slab2 model data, available to download from
# https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467
# (Choose the file "Slab2Distribute_Mar2018.tar.gz" under the
# "Attached files" heading and unpack it in this directory).
# 
# This script then formats the data in a way that it can be used by
# the slab model plugin in ASPECT.
# 
# Input files: Each input file contains attributes representing either depth
# "dep", dip ("dip"), strike ("str"), and thickness ("thk") information for
# each of the 27 slabs in the slab2 model in a structured ascii grid. The files
# contains columns corresponding to longitude, latitude, and attribute ("dep,
# "dip" etc). More information on the model can be found in the paper: 
# Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., Flamme, H., Furtney,
# M., & Smoczyk, G. M.  (2018). Slab2,a comprehensive subduction zone geometry
# model. Science, 362(6410), 58-61.
# 
# Output: A structured file containing all the slabs for the Earth with columns
# representing longitude, colatitude, slab depths, and slab thickness

import numpy as np
import os
from glob import glob
from scipy.interpolate import griddata

# Save the current working directory
cwd = os.getcwd()

# Set the path to the downloaded slab2 model here. This path is correct if you
# download the model from the address above and extract in the current folder.
os.chdir('./Slab2Distribute_Mar2018/Slab2_TXT')


def interpolate_data (name_key, type_key):
    """
    From initial data inspection, we found that most slabs have spacing of 0.05 degrees except for 5 slabs,   
    which have 0.02 degree spacing. Therefore, this function interpolates those slabs to 0.05 spacing for 
    consistency.
    """

    file_name = glob(name_key + '*' + type_key +'*.xyz')
    data_in   = np.loadtxt (file_name[0], delimiter=",")
    
    # We convert latitudes to co-latitudes for all the slab files.
    points_x  = data_in[:, 0] 
    points_y  = 90 - data_in[:, 1]
    field     = data_in[:, 2]
    
    x0 = np.min(points_x)
    x1 = np.max(points_x)
    y0 = np.min(points_y)
    y1 = np.max(points_y)
    
    # create the desired mesh
    x = np.arange (x0, x1 + 0.05, 0.05)
    y = np.arange (y0, y1 + 0.05, 0.05)
    
    x_grid, y_grid = np.meshgrid (x, y)
    
    depth_interp = griddata ((points_x, points_y), field, (x_grid, y_grid), method = 'linear')
    data_save    = np.column_stack ((x_grid.flatten(), 90 - y_grid.flatten(), depth_interp.flatten()))
    
    np.savetxt (file_name[0], data_save, fmt='%1.2f', delimiter=',')


# Slabs with 0.02 degree spacing
file_name_key = ['hin', 'man', 'mue', 'pam', 'puy']
file_type_key = ['dep', 'thk', 'dip']

for i in range(len(file_name_key)):
    for j in range(len(file_type_key)):
        interpolate_data (file_name_key[i], file_type_key[j])

# Create the grid for aspect, with global longitudes and colatitudes and
# resolution same as all the slabs (0.05 degrees)
longitudes_target = np.arange (0 , 360.05, 0.05)
latitudes_target  = np.arange (0,  180.05, 0.05)
longitude_grid, latitude_grid = np.meshgrid (longitudes_target, latitudes_target)

# We initialize the depth array with a very large number and the
# thickness value with 0 to make it easy to determine regions
# without slabs.
depths_target       = np.finfo(np.float64).max / 1e4 * np.ones (len(longitude_grid.flatten()),)
thickness_target    = np.zeros (len(longitude_grid.flatten()),)
dips_target         = np.zeros (len(longitude_grid.flatten()),)

nan_indx = 0 
def slab2_into_asciiboundary (slab_filename):
    """
    This function uses the above target mesh and fills it with the input slab depths, thicknesses
    and dip angles.
    Since the input data grid for each slab is a subset of this target grid, we can simply replace 
    the target grid at those locations with the input slab grid. 
    We have checked that in the input files longitude increases first and then latitude, suggesting that the 
    elements are ordered identically to the target mesh.
    """
    
    # load the input slab depth files
    slab_data = np.loadtxt (slab_filename, delimiter=",")
    points_x  = slab_data[:, 0] 
    points_y  = 90 - slab_data[:, 1]
    depths    = slab_data[:, 2]    
    
    x0 = np.min(points_x)
    x1 = np.max(points_x)
    y0 = np.min(points_y)
    y1 = np.max(points_y)
    
    # need this to compare two floats
    numerical_error = 1e-6
   
    indx = np.where ( (longitude_grid.flatten() >= x0 - numerical_error) & 
                      (longitude_grid.flatten() <= x1 + numerical_error) &
                      (latitude_grid.flatten()  >= y0 - numerical_error) & 
                      (latitude_grid.flatten()  <= y1 + numerical_error) )
    
    # We do not simply replace all the points with the slab grid. This is because there 
    # are some slabs in the western US with overlapping grid points. However, locations 
    # in the slab grid where depths and thickness are not defined are nan values, and we
    # want to avoid replacing existing depth and thickness values with nan values 
    # from another slab grid.
    # Therefore, we create a dummy variable that stores all the grid points in the slab 
    # and then copies only the non-nan values into the target array.
    depths_temp       = np.empty (len(longitude_grid.flatten()),)  
    depths_temp[:]    = np.nan
    depths_temp[indx] = slab_data[:, 2]
    
    # replace the target mesh with the slab data
    indx_not_a_nan                = np.where(~np.isnan (depths_temp))
    depths_target[indx_not_a_nan] = depths_temp[indx_not_a_nan]
    
    # do the same for the thickness file
    thk_file  = slab_filename.replace ('dep', 'thk')
    thk_data  = np.loadtxt (thk_file, delimiter=",", usecols=2)
    
    # do the same for the dip file
    dip_file  = slab_filename.replace ('dep', 'dip')
    dip_data  = np.loadtxt (dip_file, delimiter=",", usecols=2)

    # We checked that slab thickness and dips are defined at points where
    # slab depths are defined.
    target_indx_not_a_nan            = np.where(~np.isnan (depths))
    thickness_target[indx_not_a_nan] = thk_data[target_indx_not_a_nan]
    dips_target[indx_not_a_nan]      = dip_data[target_indx_not_a_nan]

    return (depths_target, thickness_target, dips_target)


input_filenames = glob('*_slab2*' + file_type_key[0]  + '*.xyz')

for i in range(len(input_filenames)):
    print (input_filenames[i])
    slab_depths, slab_thickness, slab_dips = slab2_into_asciiboundary (input_filenames[i])

# The thickness in the slab2 model is perpendicular to the slab surface, therefore, we need
# to convert it into equivalent distance vertically from the slab surface.
thickness_vertical = slab_thickness/np.cos(np.deg2rad(slab_dips))

# Save the arrays in the format for aspect ascii boundary :
# 1. use radians in longitudes and co-latitudes, 2. first longitudes increase then co-latitudes
# 3. convert depth and thickness from km to meters 4. depth is negative in the downward
# direction of the slab2 database, convert to positive values. 5. Add necessary header
# information.

km_to_m          = 1e3
output_dir       = cwd
output_filename  = '/slab2_depth_thickness_highres.txt'
output_data      = np.column_stack ((np.deg2rad(longitude_grid.flatten()), np.deg2rad(latitude_grid.flatten()),
                                     abs(depths_target.flatten())*km_to_m, thickness_vertical.flatten()*km_to_m))

aspect_header    = 'Test data for ascii data initial conditions.\n' + \
                   'Only next line is parsed in format: [nx] [ny] because of keyword "POINTS"\n' + \
                   'POINTS: %d %d\n' %(np.size(longitudes_target), np.size(latitudes_target)) + \
                   'Columns: phi theta depth(m) thickness(m)'

np.savetxt (output_dir + output_filename, output_data, header=aspect_header, fmt='%.4e')
