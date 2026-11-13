# This python script will reproduce the figures for the rayleigh_taylor_instability benchmark,
# which should be run prior to executing this code. The script assumes the results are located
# in the output folder (output) one directory level above this folder. The script
# can be modified to plot different time steps by changing the pvtu file defined in the variable
# 'solution_file' below. It is recommended to run this script using the conda environment provided
# in contrib/python/env-py_aspect.yml

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cmcrameri import cm

# Helper function. This might be a good candidate to include in gdmate.
def camera_translate(plotter,translation_delta):
    old_position = plotter.camera_position
    # camera_position returns an array of tuples. 
    # Need to convert the tuple to a list to modify its contents.
    translation = list(old_position[0])
    translation[0] += translation_delta[0]
    translation[1] += translation_delta[1]
    translation[2] += translation_delta[2]
    translation_tuple = tuple(translation)

    focal = list(old_position[1])
    focal[0] += translation_delta[0]
    focal[1] += translation_delta[1]
    focal[2] += translation_delta[2]
    focal_tuple = tuple(focal)

    new_position = [translation_tuple,
    focal_tuple,
    old_position[2]]

    plotter.camera_position = new_position

# set this variable to False if you are rendering grid1.png
render_colorbar = True

pv.set_plot_theme("document")

# inputs/config
solution_file = "../output/solution/solution-00006.pvtu"
output_png    = "../doc/grid2.png"


mesh = pv.read(solution_file)
size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])
mesh_grid = mesh.extract_all_edges()

plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")
plotter.add_mesh(mesh, scalars="density", cmap=cm.vik,lighting=False, show_scalar_bar=False)
plotter.add_mesh(mesh_grid, style="wireframe", color="black", line_width=2, opacity=1)

# straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
plotter.camera_position = "xy"
plotter.enable_parallel_projection()
plotter.camera.parallel_scale = size/10
camera_translate(plotter,(-55000,0,0))
plotter.screenshot(output_png)
img = plotter.screenshot(return_img=True)
plotter.close()

# for the grid2 image, pyvista only renders the density surface.
# The whitespace to the right and the colorbar are handled with
# matplotlib. You can see additional examples of this pattern in
# cookbooks/composition-reaction/plotting_scripts/plot_figure.py and other
# plotting scripts.
if (render_colorbar):
    # Dividing image_height by right_offset produces an aspect ratio
    # which fits the rendered image neatly into the window, only 
    # leaving whitespace on the right side after adjusting the subplot.
    # the value of image_height here is arbitrary,
    # changing it will change the size of the final image.
    right_offset = 0.90
    image_height = 12
    fig, ax = plt.subplots(figsize=(image_height/right_offset, image_height))

    fig.subplots_adjust(left=0.00, right=right_offset, bottom=0.0, top=1)
    ax.set_axis_off()
    ax.imshow(img)

    # This normalizes the colormap used for the colorbar.
    # If you don't do this, you won't be able to use tick
    # values outside of 0.0-1.0.
    normalize = mcolors.Normalize(vmin=3000, vmax=3300)
    scalable = plt.cm.ScalarMappable(norm=normalize,cmap=cm.vik)
    
    cax = fig.add_axes([0.93, 0.25, 0.018, 0.50])
    cb = fig.colorbar(scalable, cax=cax,
                    ticks=[3.0e3, 3050, 3100, 3150, 3200, 3250, 3.3e3])
    cb.ax.set_title("Density\n (Kg/m^3)", fontsize=13, pad=18)
    cb.ax.tick_params(labelsize=13)

    plt.savefig(output_png, dpi=220,pad_inches=0.0000001, facecolor="white")
    plt.close(fig)
