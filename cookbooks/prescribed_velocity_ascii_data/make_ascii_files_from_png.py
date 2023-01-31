# This script reads a grayscale image called aspect_name.png, and
# manipulates the pixel values to create tables of composition
# (scaled between 0 and 1) and temperature.
# The resulting composition and temperature tables are then plotted
# and saved to ASCII text files that can be read by ASPECT.

# This script can be called with the command
# python ./make_ascii_files_from_png.py

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# Open grayscale image
image = Image.open('aspect_name.png')
array = np.array(image.getdata())

# Get only the grayscale bit and flip along the y axis
nx, ny = image._size
grayscale = np.flip(array[:, 0].reshape((ny, nx)),
                    0).flatten()

# Create grids
x_spacing = 1000.  # m/pixel
y_spacing = 1000.  # m/pixel

xs = np.linspace(0., x_spacing*(nx-1), nx)
ys = np.linspace(0., y_spacing*(ny-1), ny)
xygrid = np.meshgrid(xs, ys)

xv = xygrid[0].flatten()
yv = xygrid[1].flatten()

# Create composition and temperature data from pixel values
compositions = 1. - (grayscale/255.)
composition_data = np.array([xv, yv, compositions]).T

temperatures = 1500. + 2.*compositions * (yv/1000. - 250.)
temperature_data = np.array([xv, yv, temperatures]).T

# Save ASCII files
header = f'POINTS: {nx} {ny}\nColumns: x y phase'
np.savetxt('aspect_name_initial_composition.txt',
           header=header, X=composition_data, fmt='%.1f')

np.savetxt('aspect_name_initial_temperature.txt',
           header=header, X=temperature_data, fmt='%.1f')

# Plot initial composition and temperature
xs_pad = np.linspace(0., x_spacing*(nx), nx+1)/1000
ys_pad = np.linspace(0., y_spacing*(ny), ny+1)/1000.

fig = plt.figure(figsize=(16, 5))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

for i in range(2):
    ax[i].set_aspect(1.)
    ax[i].set_xlim(0., 800.)
    ax[i].set_ylim(0., 500.)
    ax[i].set_xlabel('x (km)')
    ax[i].set_ylabel('y (km)')

cmap = plt.get_cmap('cividis')
msh = ax[0].pcolormesh(xs_pad, ys_pad,
                       compositions.reshape((ny, nx)),
                       cmap=cmap)
fig.colorbar(msh, ax=ax[0])
ax[0].set_title('Composition')

cmap = plt.get_cmap('RdBu_r')
msh = ax[1].pcolormesh(xs_pad, ys_pad,
                       temperatures.reshape((ny, nx)),
                       cmap=cmap)
fig.colorbar(msh, ax=ax[1])
ax[1].set_title('Initial temperature (K)')

fig.savefig('prescribed_velocities_ascii_data_initial_conditions.png')
plt.show()
