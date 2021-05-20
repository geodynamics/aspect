import numpy as np

# read in data files
upper_mantle = np.loadtxt('rent.d')
lower_mantle = np.loadtxt('tmelt.d')

# divide H by nR
upper_mantle[:,0] /= 3.5 * 8.31446

# convert from melting temperature to H/nR
lower_mantle[:,0] *= 12.0

# remove the last row from the upper mantle array to not have the 670 km depth value twice
# then attach the two together
viscosity = np.concatenate((upper_mantle[:-1, :], lower_mantle))

# change depth to positive values
viscosity[:,1] *= -1.0
 
np.savetxt('temp-viscosity-prefactor.txt', viscosity, fmt='%.7g', header='# H/nR [K], depth [km]', comments='')

