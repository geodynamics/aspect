# This python script will reproduce the figures for the geomio benchmark,
# which should be run prior to executing this code. The script assumes the results are located
# in the output folder (output-geomIO) one directory level above this folder. 
# It is recommended to run this script using the conda environment provided
# in contrib/python/env-py_aspect.yml

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

pv.set_plot_theme("document")

# inputs/config
solution_file = "../output-geomIO/solution/solution-00000.pvtu"
output_svg    = "../doc/jelly-paraview.svg"

mesh = pv.read(solution_file)
size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])
mesh_grid = mesh.extract_all_edges()

annotations = {
    0: "0",
    1: "1",
    2: "2",
    3: "3"
}


sargs = dict(
    title='Phase',
    vertical=True,
    position_x=0.03, 
    position_y=0.05,
    title_font_size=30,
    label_font_size=30,
    n_colors=4,
    n_labels=0,
    width=0.15,
    height=0.80,
)

plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")
plotter.add_mesh(mesh, scalars="C_1", cmap=cm.batlow,lighting=False, show_scalar_bar=True, scalar_bar_args=sargs, annotations=annotations)
plotter.add_mesh(mesh_grid, style="wireframe", line_width=1, opacity=0.25)

# straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
plotter.camera_position = "xy"
plotter.enable_parallel_projection()
plotter.camera.parallel_scale = size/2
plotter.save_graphic(output_svg)
plotter.close()
