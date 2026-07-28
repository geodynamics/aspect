# It is recommended to run this script using the conda environment provided in 
# contrib/python/env-py_aspect.yml

import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from cmcrameri import cm # Fabio Crameri's color maps
import numpy as np
import os

pv.set_plot_theme("document")
plotter = pv.Plotter(off_screen=True)

# As is, reading in solution-00000.pvtu will plot the model at the first timestep.
# To plot other time steps, change the path to point at a different .pvtu file.
# For example, to generate an image at time step 1, change to path to: 
# "../../output-convection-box/solution/solution-00001.pvtu"
# Remember to change the name of the output file as well.
mesh = pv.read("../output-convection-box/solution/solution-00004.pvtu")
plot_spatial_bounds = [0, 1, 0, 1, 0, 0]
mesh = mesh.clip_box(bounds=plot_spatial_bounds, invert=False)
# calculate median velocity for arrow scaling
vel = mesh["velocity"]
vmag = np.linalg.norm(vel, axis=1)
median_v = float(np.median(vmag[vmag > 0])) if np.any(vmag > 0) else 1.0 # median of nonzero velocities, or 1.0 if all velocities are zero

arrows = mesh.glyph(scale="velocity", factor=0.025 / median_v,orient="velocity")

# scale the meshes down to put the scale bar on the side (like in the docs)
mesh = mesh.scale(0.75)
arrows = arrows.scale(0.75)

# Define properties for the scalar bar arguments
sargs = dict(
    width=0.2,
    height=0.6,
    title="T",
    label_font_size=24,
    title_font_size=32,
    color="black",
    position_x=0.0,
    position_y=0.20,
    vertical=True
)

mesh_actor = plotter.add_mesh(mesh, scalars="T", log_scale=False, scalar_bar_args=sargs,cmap=cm.vik)
arrows_actor = plotter.add_mesh(arrows,color='white')
plotter.view_xy()

# Displace the meshes to make room for the scalar bar
mesh_actor.position = (0.225,0.15,0)
arrows_actor.position = (0.225,0.15,0)

# shows the axes labels on the actual mesh (show_axes shows a widget in 3d space)
plotter.show_bounds(
use_3d_text=False,use_2d=True, # These parameters make it easier to plot 2D data
axes_ranges=plot_spatial_bounds, # We displaced mesh_actor earlier, so we need to reset the axes ranges
grid='front',ticks='outside',location='outer',
xtitle='X Axis',ytitle='Y Axis'
)

# Calculate Camera Position from Bounds
xmag = float(abs(plot_spatial_bounds[1] - plot_spatial_bounds[0]))
ymag = float(abs(plot_spatial_bounds[3] - plot_spatial_bounds[2]))

# xmag is multiplied by 1.2 here in order to change the aspect ratio of the 
# resulting image, in this case resulting in less white space on the 
# top and bottom of the image.
aspect_ratio = ymag / (xmag*1.2) 
plotter.window_size = (1024, int(1024 * aspect_ratio))
xmid = (plot_spatial_bounds[1] - plot_spatial_bounds[0]) /2
ymid = (plot_spatial_bounds[3] - plot_spatial_bounds[2]) /2
# zoom_level controls how close the camera should be to the mesh.
# A smaller value will move the camera closer resulting in a more
# zoomed in image.
zoom_level = 1.875
zoom = xmag * aspect_ratio * zoom_level  

position = (xmid, ymid, zoom)
focal_point = (xmid, ymid, 0)
viewup = (0, 1, 0)
camera = [position, focal_point, viewup]
plotter.camera_position = camera
plotter.camera_set = True
plotter.screenshot("../doc/visit0001.png")
