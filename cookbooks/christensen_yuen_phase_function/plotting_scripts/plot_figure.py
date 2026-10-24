# This python script will reproduce the figures for the
# christensen_yuen_phase_function cookbook, which should
# be run before executing this code. The script assumes 
# the results are located in the
# output folders (christensen_yuen) and (christensen_yuen_2) 
# one level above this folder. 
# 
# In order to reproduce the figure
# from the documentation, you will need to run the model once 
# with a Clapeyron slope of -2700000 
# (which is the starting value in christensen_yuen_phase_function.prm)
# and then again edit the .prm to have a Clapeyron slope of +2700000.
# You should also edit the output directory in the .prm so as not to 
# overwrite the first run. This script assumes that the second run
# outputs to christensen_yuen_2.
# 
# It is recommended to run this script using the conda environment provided
# in contrib/python/env-py_aspect.yml.

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

pv.set_plot_theme("document")

# inputs/config
solution_file_negative_clapeyron = "../christensen_yuen/solution/solution-00538.pvtu"
solution_file_positive_clapeyron = "../christensen_yuen_2/solution/solution-00186.pvtu"
output_png    = "../doc/flow_field.png"

# This is the range of temperatures from the input file
t_clim = [0, 1000]

# Helper function
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

# Helper function. Returns a 2D image to compose into a composite image.
def render_2d(mesh, scalars, clim, cmap, draw_arrows=True, draw_scalarbar=True, sargs=None):
    # These variables are used to scale the arrow glyphs and
    # scalar bar ticks properly. In some examples,
    # like convection-box example, the mesh has a size of 1 so this isn't needed.
    bounds_original = mesh.bounds
    size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])

    plotter = pv.Plotter(off_screen=True)
    plotter.window_size = (1024, 1024)
    plotter.set_background("white")

    # Make room for the scalar bar if there is one
    if (draw_scalarbar == True):
        mesh = mesh.scale(0.80)

    if (draw_arrows == True):
        # median velocity for arrow scaling
        vel = mesh["velocity"]
        vmag = np.linalg.norm(vel, axis=1)
        median_v = float(np.median(vmag[vmag > 0])) if np.any(vmag > 0) else 1.0
        arrows = mesh.glyph(scale="velocity", factor=(0.06*size) / median_v, orient="velocity", tolerance=0.03)
        glyph_actor = plotter.add_mesh(arrows, cmap=cm.imola, show_scalar_bar=False)
        # Setting use_bounds to false here allows the camera to cleanly crop out the arrows
        # when using plotter.camera.tight(). This is noticeable when not drawing
        # scalar bars.
        glyph_actor.use_bounds = False

    # This enables you to pass in a different value 
    # for scalar bar arguments. By default the function
    # assumes you are plotting the temperature surface
    # with units of K.
    if (sargs is None):
        sargs = dict(
            title="Temperature (K)",
            height=0.13,
            width=0.75,
            position_x=0.125,
            position_y=0.01,
            title_font_size=40,
            label_font_size=32,
            fmt='%.2f'
        )

    mesh_actor = plotter.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=clim,
                    lighting=False, show_scalar_bar=draw_scalarbar,scalar_bar_args=sargs)
    
    # straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
    plotter.camera_position = "xy"
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = size / 2

    # shows the axes labels on the actual mesh (show_axes shows a widget in 3d space)
    plotter.show_bounds(
    use_3d_text=False,use_2d=True, # These parameters make it easier to plot 2D data
    axes_ranges=bounds_original, # We scaled the mesh earlier, so we need to set the axes ranges to their true values
    location='outer',
    ticks='outside',
    xtitle='X Axis (m)',ytitle='Y Axis (m)',
    grid=None,
    fmt='%.0f',
    n_xlabels=6,
    n_ylabels=6,
    )

    # If drawing the scalar bar, move the camera upwards to
    # make room for the scalar bar under it. Need to move the
    # camera after setting camera_position to "xy", not before.
    if (draw_scalarbar):
        camera_translate(plotter,[0,-0.09*size,0])

    img = plotter.screenshot(return_img=True)
    plotter.close()
    return img


# load the solutions
mesh = pv.read(solution_file_negative_clapeyron)
mesh2 = pv.read(solution_file_positive_clapeyron)

# -2.7 Mpa K^-1
negative_clapyron_img = render_2d(mesh, scalars="T", clim=t_clim, cmap=cm.vik, draw_arrows=True, draw_scalarbar=True)

# +2.7 Mpa K^-1
positive_clapyron_img = render_2d(mesh2, scalars="T", clim=t_clim, cmap=cm.vik, draw_arrows=True, draw_scalarbar=True)

# compose the two panels side by side and save
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
for ax, img, title in zip(axes, [negative_clapyron_img, positive_clapyron_img],
                          ["Clapeyron slope = -2.7 MPa $K^{-1}$", "Clapeyron slope = 2.7 MPa $K^{-1}$"]):
    ax.imshow(img)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24)

plt.tight_layout()
plt.savefig(output_png, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)
