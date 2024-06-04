# This python script uses the pyvista package to do math and plot
# data from a 2D cylindrical annulus output. The parameter file
# used to generate the output figures can be made by changing the
# onset-of-convection cookbook example in the ways described by the
# comments at the end of this file.

# This script uses TimeReader class.
# This class reads in .pvd files to extract the time information
# about a model output and reads .pvtu files for the T,v,p solution.
# The script also showcases simple mesh operations as well as
# the thresholding and contouring utilities of pyvista.

# For additional overview on the utility of pyvista with ASPECT
# please read the description from John Naliboff in
# contrib/python/scripts/pyvista_2d_example.py.

# Load modules
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

#   Include the path to the project directory you would like to
#   visualize.
dir = "/2D_annulus_simple_convection"  # [CHANGE THIS TO MATCH YOUR FILE SYSTEM]

#   Time step: initial and final values of the output file numbers
#    to traverse in this visualization
out_n_min = 0
out_n_max = 25

#   Set the 2D annulus dimensions in meters
r_inner = 3481000  # m
r_outer = 6371000  # m

#   Create a visualization directory in model directory for
#   saving the images.
os.makedirs(dir + "plots", exist_ok=True)  # make plot save directory if doesnt exist

#   Read in the time values associated with the solution output files.
#   More documentation about pyvista Time Reader at:
#   https://docs.pyvista.org/version/stable/api/readers/_autosummary/pyvista.timereader
timereader = pv.get_reader(dir + "output/solution.pvd")
times = np.array(timereader.time_values) / 1e6  # convert yrs to Myrs

#   Loop through the solution output files
for ts in range(out_n_min, out_n_max + 1):
    print("timestep number:", ts, "time =", times[ts], "Myrs")

    #   Read data from the pvtu file into a variable named mesh, which contains
    #   all information needed for the plot.
    num_str = str(ts).zfill(5)
    mesh = pv.read(dir + "output/solution/solution-" + num_str + ".pvtu")

    #   Extract cartesian coordinates of the mesh and x-y velocities
    mesh.point_data["x"] = mesh.points.T[0] # equivalently: [z[0] for z in mesh.points]
    mesh.point_data["y"] = mesh.points.T[1] # equivalently [z[1] for z in mesh.points]
    mesh.point_data["vx"] = mesh.point_data["velocity"].T[0]
    mesh.point_data["vy"] = mesh.point_data["velocity"].T[1]

    #   Compute polar coordinates of the mesh and radial
    #   and angular velocities
    mesh.point_data["r"] = np.sqrt(mesh["x"] ** 2 + mesh["y"] ** 2)
    mesh.point_data["phi"] = np.arctan2(mesh["y"], mesh["x"])
    mesh.point_data["v_r"] = (mesh["x"] * mesh["vx"] + mesh["y"] * mesh["vy"]) / mesh[
        "r"
    ]
    mesh.point_data["v_phi"] = (mesh["x"] * mesh["vy"] - mesh["y"] * mesh["vx"]) / mesh[
        "r"
    ]

    #   Define mid-mantle region, 400 km above the core-mantle boundary and
    #   400 km below the surface to avoid the temperature effects of the the non-
    #   mixed cold and hot thermal boundary layers on the boundaries.
    mesh_midmantle = mesh.threshold(
        scalars="r", value=r_inner + 400000, invert=False
    ).threshold(scalars="r", value=r_outer - 400000, invert=True)

    #   Compute the mean temperature value of the mid-mantle. This can be used
    #   along with the core and surface temperatures to define the thermal boundary
    #   layers
    mean_mantle_T = np.mean(mesh_midmantle["T"])

    #   Create a new point field in the mesh object which is the temperature
    #   deviation from the mean at all points in the domain.
    mesh.point_data["T_dev"] = mesh["T"] - mean_mantle_T

    #    Contours can be drawn and data can be extracted along curves of constant radius
    midmantle_contour = mesh.contour(scalars="r", isosurfaces=[(r_outer + r_inner) / 2])

    #   Plot the angular variations of the temperature field
    fig, ax = plt.subplots(1)
    fig.set_size_inches(12, 3)
    ax.scatter(midmantle_contour["phi"], midmantle_contour["T"])
    plt.show()
    plt.clf()
    plt.close()

    pv.set_plot_theme("document")
    #   Define formatting for the scalar color bar of the plot images
    sargs_vr = dict(
        width=0.6,
        fmt="%.1e",
        height=0.2,
        title="radial velocity",
        label_font_size=24,
        title_font_size=32,
        color="black",
        position_x=0.2,
        position_y=0.01,
    )

    #   Plot radial velocity
    pl = pv.Plotter(off_screen=True)
    pl.add_text(str(times[ts]) + " Myrs", position="upper_left", font_size=18)
    #   Add a mesh to the plot with scalar color bar limits, a colormap, and lighting=False since
    #   we are in 2D and don't care about directional lighting.
    pl.add_mesh(
        mesh,
        scalars="v_r",
        clim=[-1e-2, 1e-2],
        cmap="RdBu",
        lighting=False,
        scalar_bar_args=sargs_vr,
    )
    pl.view_xy()
    pl.screenshot(dir + "plots/vr_" + num_str + ".png")
    pl.clear()
    pl.close()

    sargs_vphi = dict(
        width=0.6,
        fmt="%.1e",
        height=0.2,
        title="angular velocity",
        label_font_size=24,
        title_font_size=32,
        color="black",
        position_x=0.2,
        position_y=0.01,
    )
    #   Plot angular velocity
    pl = pv.Plotter(off_screen=True)
    pl.add_text(str(times[ts]) + " Myrs", position="upper_left", font_size=18)
    pl.add_mesh(
        mesh,
        scalars="v_phi",
        clim=[-1e-2, 1e-2],
        cmap="RdBu",
        lighting=False,
        scalar_bar_args=sargs_vphi,
    )
    pl.view_xy()
    pl.screenshot(dir + "plots/vphi_" + num_str + ".png")
    pl.clear()
    pl.close()

    sargs_T = dict(
        width=0.6,
        fmt="%.1e",
        height=0.2,
        title="T-mean(T)",
        label_font_size=24,
        title_font_size=32,
        color="black",
        position_x=0.2,
        position_y=0.01,
    )
    #   Plot temperature deviation from the mean mantle temperature
    pl = pv.Plotter(off_screen=True)
    pl.add_text(str(times[ts]) + " Myrs", position="upper_left", font_size=18)
    pl.add_mesh(
        mesh,
        scalars="T_dev",
        clim=[-600, 600],
        cmap="RdBu_r",
        lighting=False,
        scalar_bar_args=sargs_T,
    )
    pl.add_mesh(midmantle_contour, color="red", line_width=2)
    pl.view_xy()
    pl.screenshot(dir + "plots/T_dev_" + num_str + ".png")
    pl.show
    pl.clear()
    pl.close()


#   ASPECT model output that can utilize this script can be created by making the following changes to the onset-of-convection cookbook:
#   https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/onset_of_convection/doc/onset_of_convection.html

# <   set Model name = box
# ---
# >   set Model name = spherical shell


# <   subsection Box
# <     set X extent = 9.424777e6 # pi * 3e6
# <     set Y extent = 3e6
# <     set X repetitions = 3
# <     set Y repetitions = 1
# ---
# >   subsection Spherical shell
# >     set Outer radius = 6371000
# >     set Inner radius = 3481000

# <     set Variable names      = x,z
# <     set Function constants  = p=1, L=9.424777e6, H=3.0e6, pi=3.1415926536, k=2
# <     set Function expression = 2500 * (1.0-z/H) - p*cos(k*pi*x/L)*sin(pi*z/H)
# ---
# >     set Coordinate system = spherical
# >     set Variable names = r,phi
# >     set Function expression = 1600 + 10*sin(10*phi)

# <   set List of model names = box
# ---
# >   set List of model names = spherical constant

# <   subsection Box
# <     set Bottom temperature = 2500
# <     set Left temperature   = 0
# <     set Right temperature  = 0
# <     set Top temperature    = 0
# ---
# >   subsection Spherical constant
# >     set Inner temperature = 2000
# >     set Outer temperature = 0

# <   set Tangential velocity boundary indicators = left, right, bottom, top
# ---
# >   set Tangential velocity boundary indicators = bottom, top

# <   set Model name = vertical
# ---
# >   set Model name = radial constant

# <   subsection Vertical
# ---
# >   subsection Radial constant
