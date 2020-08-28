import os
import numpy as np

"""
This standalone python program reads in the four Crust1 files
obtainable from the Crust1 website https://igppweb.ucsd.edu/~gabi/crust1.html
(crust1.bnds, crust1.rho, crust1.vp, crust1.vs). It outputs ten
files containing the information stored for each boundary.
"""

boundary_topography = np.loadtxt('crust1.bnds')*1000.
r = 6371.e3 + boundary_topography
rho = np.loadtxt('crust1.rho')*1000.
vp = np.loadtxt('crust1.vp')*1000.
vs = np.loadtxt('crust1.vs')*1000.

lons = np.linspace(-179.5, 179.5, 360)
lats = np.linspace(89.5, -89.5, 180)

lons *= np.pi/180.
lats = (90. - lats)*np.pi/180.

lonlon, latlat = np.meshgrid(lons, lats)


lonlon = lonlon.flatten()
latlat = latlat.flatten()

layers = ['water',
          'ice',
          'upper_sediments',
          'middle_sediments',
          'lower_sediments',
          'upper_crust',
          'middle_crust',
          'lower_crust',
          'mantle']

for i in range(len(rho[0])):
    print('Writing file for the {0} layer...'.format(layers[i]))
    data = np.array([lonlon, latlat,
                     r[:, i],
                     [float(i+1)] * len(r[:, i]),
                     rho[:, i],
                     vp[:, i],
                     vs[:, i]]).T

    header = ('Crust1 data for the top of the {0} layer (Layer {1}).\n'
              'Only next line is parsed in format: [nx] [ny] '
              'because of keyword \"POINTS:\"\n'
              'POINTS: 360 180\n'
              'Columns: phi theta r layer_id density vp vs'.format(layers[i],
                                                                   i+1))
    fname = 'ascii_data_crust1_{0}_top_{1}.txt'.format(i+1, layers[i])
    np.savetxt(fname=fname, X=data, fmt='%.5f', header=header)


# Base air
i = 0
data = np.array([lonlon, latlat,
                 r[:, i],
                 [0.] * len(r[:, i]),
                 [1.] * len(r[:, i]),
                 [0.] * len(r[:, i]),
                 [0.] * len(r[:, i])]).T

header = ('Crust1 data for the base of the air layer (Layer 0).\n'
          'Only next line is parsed in format: [nx] [ny] '
          'because of keyword \"POINTS:\"\n'
          'POINTS: 360 180\n'
          'Columns: phi theta r layer_id density vp vs')

fname = 'ascii_data_crust1_0_base_air.txt'
np.savetxt(fname=fname, X=data, fmt='%.5f', header=header)

print('Files successfully written.')
