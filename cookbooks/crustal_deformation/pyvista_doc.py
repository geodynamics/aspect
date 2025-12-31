import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from cmcrameri import cm # Fabio Crameri's color maps
import numpy as np
import os

timereader = pv.get_reader("./output-crustal_model_2D/solution.pvd")
times = np.array(timereader.time_values) / 1e6  # convert yrs to Myrs

pv.set_plot_theme("document")
plotter = pv.Plotter(off_screen=True)

mesh = pv.read("output-crustal_model_2D/solution/solution-00020.pvtu")
plot_spatial_bounds = [0, 100e3, 0, 100e3, 0, 0]

# Define properties for the scalar bar arguments
sargs = dict(
    width=0.6,
    height=0.1,
    title="Viscosity (Pa s)",
    label_font_size=24,
    title_font_size=32,
    color="black",
    position_x=0.20,
    position_y=0.05,
    n_labels=2,
    vertical=False
)

mesh_actor = plotter.add_mesh(mesh, scalars="viscosity", log_scale=True, scalar_bar_args=sargs,cmap=cm.batlow_r)

# Displace the mesh to make room for the scalar bar
mesh_actor.position = (10e3,25e3,0)
plotter.view_xy()

plotter.show_bounds(
use_3d_text=False,use_2d=True, # These parameters make it easier to plot 2D data
axes_ranges=[0,80,0,16,0,0], # We displaced mesh_actor earler, so we need to reset the axes ranges
grid='front',ticks='outside',location='outer',
n_ylabels=3,n_xlabels=9,bold=False,
xtitle='X Axis (km)',ytitle='Y Axis (km)'# Model bounds are in meters, kilometers make the labels more readable
)

# Calculate Camera Position from Bounds
bounds_array = np.array(plot_spatial_bounds)
xmag = float(abs(bounds_array[1] - bounds_array[0]))
ymag = float(abs(bounds_array[3] - bounds_array[2]))
aspect_ratio = ymag / xmag
plotter.window_size = (1024, int(1024 * aspect_ratio))
xmid = xmag / 2 + bounds_array[0]  # X midpoint
ymid = ymag / 2 + bounds_array[2]  # Y midpoint
zoom = xmag * aspect_ratio * 1.875  # Zoom level

position = (xmid, ymid, zoom)
focal_point = (xmid, ymid, 0)
viewup = (0, 1, 0)
camera = [position, focal_point, viewup]
plotter.camera_position = camera
plotter.camera_set = True
plotter.screenshot("fig108.png")
