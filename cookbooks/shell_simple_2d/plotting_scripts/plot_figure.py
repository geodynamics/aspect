# This python script will reproduce the figures for the shell_simple_2d cookbook,
# which should be run prior to executing this code. The script assumes the results are located
# in the output folder (output-shell_simple_2d) one directory level above this folder. The script
# can be modified to plot different time steps by changing the pvtu file defined in the variable
# 'solution_file' below. It is recommended to run this script using the conda environment provided
# in contrib/python/env-py_aspect.yml

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

pv.set_plot_theme("document")

# inputs/config
solution_file = "../output-shell_simple_2d/solution/solution-01000.pvtu"
output_png    = "../doc/x-movie1000.png"

# fixed temperature range (inner/outer boundary temperatures from the .prm) so every
# time step shares the same colors and the panels can be compared side by side
clim = [973, 4273]

# show the mesh over the lower part of the annulus only
mesh_max_angle = 35   # degrees

# load the solution
mesh = pv.read(solution_file)
size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])

# pick the cells in the lower sector of the annulus for the mesh overlay
centers = mesh.cell_centers().points
angle = np.degrees(np.arctan2(centers[:, 1], centers[:, 0]))
sector = mesh.extract_cells(angle < mesh_max_angle)

# pyvista renders the bare temperature field; the shared colorbar lives in the svg
plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")
plotter.add_mesh(mesh, scalars="T", cmap=cm.lajolla, clim=clim,
                 lighting=False, show_scalar_bar=False)
plotter.add_mesh(sector, style="wireframe", color="ivory_black", line_width=0.3, opacity=0.1)

# straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
plotter.camera_position = "xy"
plotter.enable_parallel_projection()
plotter.camera.parallel_scale = size / 2

img = plotter.screenshot(return_img=True)
plotter.close()

# save the bare annulus with no axes or labels
fig, ax = plt.subplots(figsize=(9, 9))
ax.imshow(img)
ax.set_axis_off()
plt.savefig(output_png, dpi=100, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)
