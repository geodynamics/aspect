# It is recommended to run this script using the conda environment provided in
# contrib/python/env-py_aspect.yml
#
# This python script reproduces the figures for the
# platelike-boundary cookbook. It assumes the simulation results are located in
# the output folder (output-platelike-boundary) one directory level above this
# folder. To render a different time step, change the 'solution_file' below 
# and update 'output_png'.

import numpy as np
import pyvista as pv
from cmcrameri import cm  # Fabio Crameri's color maps

pv.set_plot_theme("document")

# inputs/config
solution_file = "../output-platelike-boundary/solution/solution-00000.pvtu"
output_png    = "../doc/visit0000.png"

# load the solution and compute the domain size used to fit the camera
mesh = pv.read(solution_file)
size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])

# calculate median velocity for arrow scaling
vel = mesh["velocity"]
vmag = np.linalg.norm(vel, axis=1)
median_v = float(np.median(vmag[vmag > 0])) if np.any(vmag > 0) else 1.0
arrows = mesh.glyph(scale="velocity", factor=0.03 / median_v, orient="velocity", tolerance=0.04)

# render the bare temperature field with velocity arrows (no scalar bar or axes)
plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")
plotter.add_mesh(mesh, scalars="T", cmap=cm.vik, lighting=False, show_scalar_bar=False)
glyph_actor = plotter.add_mesh(arrows, color="white")

# Setting use_bounds to false here allows the camera to cleanly crop out the arrows
# when using plotter.camera.tight()
glyph_actor.use_bounds = False

# straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
plotter.camera_position = "xy"
plotter.camera.tight(padding=0.0)

plotter.screenshot(output_png)
plotter.close()
