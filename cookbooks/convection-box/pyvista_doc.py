import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from cmcrameri import cm # Fabio Crameri's color maps
import numpy as np
import os

timereader = pv.get_reader("./output-convection-box/solution.pvd")
times = np.array(timereader.time_values) / 1e6  # convert yrs to Myrs

pv.set_plot_theme("document")
plotter = pv.Plotter(off_screen=True)

mesh = pv.read("output-convection-box/solution/solution-00049.pvtu")
plot_spatial_bounds = [0, 1, 0, 1, 0, 0]
mesh = mesh.clip_box(bounds=plot_spatial_bounds, invert=False)
arrows = mesh.glyph(scale="velocity", factor=0.0009,orient="velocity")

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
    position_y=0.35,
    vertical=True
)

mesh_actor = plotter.add_mesh(mesh, scalars="T", log_scale=False, scalar_bar_args=sargs,cmap=cm.vik)
arrows_actor = plotter.add_mesh(arrows,color='white')
plotter.view_xy()

# Displace the meshes to make room for the scalar bar
mesh_actor.position = (0.25,0.23,0)
arrows_actor.position = (0.25,0.23,0)

# shows the axes labels on the actual mesh (show_axes shows a widget in 3d space)
plotter.show_bounds(
use_3d_text=False,use_2d=True, # These parameters make it easier to plot 2D data
axes_ranges=[0,1,0,1,0,0], # We displaced mesh_actor earler, so we need to reset the axes ranges
grid='front',ticks='outside',location='outer',
xtitle='X Axis',ytitle='Y Axis'# Model bounds are in meters, kilometers make the labels more readable
)

# Calculate Camera Position from Bounds
bounds_array = np.array(plot_spatial_bounds)
xmag = float(abs(bounds_array[1] - bounds_array[0]))
ymag = float(abs(bounds_array[3] - bounds_array[2]))
aspect_ratio = ymag / xmag
plotter.window_size = (1024, int(1024 * aspect_ratio))
xmid = xmag / 2 + bounds_array[0]  # X midpoint
ymid = ymag / 2 + bounds_array[2]  # Y midpoint
zoom = xmag * aspect_ratio * 1.875  # Zoom level - not sure why 1.875 works

position = (xmid, ymid, zoom)
focal_point = (xmid, ymid, 0)
viewup = (0, 1, 0)
camera = [position, focal_point, viewup]
plotter.camera_position = camera
plotter.camera_set = True
plotter.screenshot("fig17.png")
