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
import matplotlib.pyplot as plt
from cmcrameri import cm  # Fabio Crameri's color maps

from aspect_pyvista import render_mesh_2d

pv.set_plot_theme("document")

# inputs/config
timesteps = [0,6,19,31,51,178]
output_combined = f"../doc/platelike.png"
output_list = []

for timestep in timesteps:
    solution_file = f"../output-platelike-boundary/solution/solution-{timestep:05d}.pvtu"
    output_png    = f"../doc/temperature{timestep:05d}.png"

    mesh = pv.read(solution_file)
    img = render_mesh_2d(mesh,scalars="T",cmap=cm.vik,arrows_enabled=True,scalarbar_enable=False,arrows_normalize=mesh["velocity"],arrows_scale=0.02,x_title=" ",y_title=" ",bounds_enable=False)
    output_list.append(img)
    # Save the individual images
    plt.imsave(output_png, img)

# compose the six panels in a grid and save
fig, axes = plt.subplots(2, 3, figsize=(14, 5))
for ax, img, title in zip(axes.flat, [output_list[0], output_list[1], output_list[2], output_list[3], output_list[4], output_list[5]],
                          ["","","","","",""]):
    ax.imshow(img)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24)

plt.tight_layout()
plt.savefig(output_combined, dpi=180, bbox_inches="tight", facecolor="white")
plt.close(fig)
