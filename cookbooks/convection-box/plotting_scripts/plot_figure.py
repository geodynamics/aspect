# It is recommended to run this script using the conda environment provided in 
# contrib/python/env-py_aspect.yml.

import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from cmcrameri import cm # Fabio Crameri's color maps
import numpy as np
import matplotlib.image as mpimg
from aspect_pyvista import render_mesh_2d
import os
import sys

pv.set_plot_theme("document")
plotter = pv.Plotter(off_screen=True)

# As is, reading in solution-00000.pvtu will plot the model at the first timestep.
# To plot other time steps, change the path to point at a different .pvtu file.
# For example, to generate an image at time step 1, change to path to: 
# "../../output-convection-box/solution/solution-00001.pvtu"
# Remember to change the name of the output file as well.
mesh = pv.read("../output-convection-box/solution/solution-00000.pvtu")

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

output_png = "../doc/visit0000.png"

# render_mesh_2d is a function from the aspect_pyvista module. The source code
# for the function can be found at contrib/python/plotting/src/aspect_pyvista.py
img = render_mesh_2d(mesh,scalars="T",grid_enabled=True,sargs=sargs,arrows_enabled=False,window_size=(1024,900),x_title="X Axis",y_title="Y Axis",n_x_labels=5,n_y_labels=5,axisFormat="%.1f")

mpimg.imsave(output_png, img)
