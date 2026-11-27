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

from aspect_pyvista import render_mesh_2d

pv.set_plot_theme("document")

# inputs/config
solution_file_negative_clapeyron = "../christensen_yuen/solution/solution-00538.pvtu"
solution_file_positive_clapeyron = "../christensen_yuen_2/solution/solution-00186.pvtu"
output_png    = "../doc/flow_field.png"

# This is the range of temperatures from the input file
t_clim = [0, 1000]

# load the solutions
mesh = pv.read(solution_file_negative_clapeyron)
mesh2 = pv.read(solution_file_positive_clapeyron)


# -2.7 Mpa K^-1
negative_clapyron_img = render_mesh_2d(mesh, scalars="T", clim=t_clim, cmap=cm.vik, arrows_enabled=True, draw_scalarbar=True, arrows_normalize=mesh2["velocity"])

# +2.7 Mpa K^-1
positive_clapyron_img = render_mesh_2d(mesh2, scalars="T", clim=t_clim, cmap=cm.vik, arrows_enabled=True, draw_scalarbar=True)

# compose the two panels side by side and save
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
for ax, img, title in zip(axes, [negative_clapyron_img, positive_clapyron_img],
                          ["Clapeyron slope = -2.7 MPa $K^{-1}$", "Clapeyron slope = 2.7 MPa $K^{-1}$"]):
    ax.imshow(img)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24)

plt.tight_layout()
plt.savefig(output_png, dpi=90, bbox_inches="tight", facecolor="white")
plt.close(fig)
