# Plots stress output from chamber_elastic_expansion.prm with the
# results from an elastic solver (Cedric Thieulot, stone 123,
# https://cedrict.github.io/) and the 3D elastic analytical solution
# from Zhong et al. 2019 (below).

# Solutions are along two 5000m transects: one at the surface
# from (0, 4000) to (5000, 4000) and a second interior transect
# from (0, 2000) to (50000, 2000).

# 2D normal and shear stresses from this cookbook are plotted
# against a cross section of the spherical solution. There are expected
# differences between the 2D ad 3D solutions, and the 3D solutions can
# be used to benchmark future 3D runs.

# Xin Zhong, Marcin Dabrowski, Bjørn Jamtveit, Analytical solution for
# the stress field in elastic half-space with a spherical pressurized cavity
# or inclusion containing eigenstrain, Geophysical Journal International,
# Volume 216, Issue 2, February 2019, Pages 1100–1115, https://doi.org/10.1093/gji/ggy447

import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# input SOLUTION_FILE can be a single .pvtu, or a .pvd (in which case TIME_STEP selects
# which timestep index to load, 0-based).
SOLUTION_FILE = "../output_chamber_elastic_expansion/solution.pvd"
TIME_STEP = 1
OUTPUT_IMAGE = "./figures/stress_plot.png"

# Stone solution profile to overlay: an x-profile
# (columns: x y z u v w p xx yy zz xy xz yz) sampled at y=0, along the surface.
# Only the stress columns (xx, yy, zz, xy, xz, yz) are used, and only for
# x in [0, 5000] to match the ASPECT transect range.
XPROFILE_FILE = "./solutions/stone_surface.ascii"
XPROFILE_LABEL = "2D FieldStone"

XPROFILE_FILE_INTERIOR = "./solutions/stone_interior.ascii"
XPROFILE_LABEL_INTERIOR = "2D FieldStone"

# (label, transect start point, transect end point, solution CSV,
# xprofile stone solution file/label to overlay on this row (or None to skip
# the overlay))
TRANSECTS = [
    ("Interior transect", (-1, 2000, 0), (5000, 2000, 0), "./solutions/reference_soln_interior.csv",
     XPROFILE_FILE_INTERIOR, XPROFILE_LABEL_INTERIOR),
    ("Surface transect", (-1, 3999, 0), (5000, 3999, 0), "./solutions/reference_soln_surface.csv",
     XPROFILE_FILE, XPROFILE_LABEL),
]

# for line resolution
N_SAMPLES = 500

# Fields to plot, can also plot ve_stress_xx.. for elastic version
# Field, index in tensor, reference CSV column, xprofile stress column, label.
# ASPECT's y (vertical) corresponds to the reference/xprofile solutions' z
# (vertical), so ASPECT's yy and xy components are compared against the
# reference/xprofile zz and xz components.
PLOTS = [
    ("stress", "xx", "sigma_xx", "xx", r"$\sigma_{xx}$"),
    ("stress", "yy", "sigma_zz", "zz", r"$\sigma_{zz}$"),
    ("stress", "xy", "sigma_xz", "xz", r"$\sigma_{xz}$"),
]

TENSOR_COMPONENT_INDEX = {"xx": 0, "xy": 1, "xz": 2, "yx": 3, "yy": 4, "yz": 5, "zx": 6, "zy": 7, "zz": 8}

# xprofile stress columns (0-based) in the ascii file: x y z u v w p xx yy zz xy xz yz
XPROFILE_STRESS_COLUMN = {"xx": 7, "yy": 8, "zz": 9, "xy": 10, "xz": 11, "yz": 12}

if SOLUTION_FILE.endswith(".pvd"):
    reader = pv.get_reader(SOLUTION_FILE)
    reader.set_active_time_point(TIME_STEP)
    mesh = reader.read()[0]
else:
    mesh = pv.read(SOLUTION_FILE)


def resolve_scalar(line, field_name, tensor_component):
    """Resolve a field (which may be a scalar, vector, or tensor array) into a single scalar array."""
    field = line[field_name]
    if field.ndim == 1:
        return field
    elif field.shape[1] == 3:
        return np.linalg.norm(field, axis=1)
    elif field.shape[1] == 9:
        if tensor_component == "magnitude":
            return np.linalg.norm(field, axis=1)
        return field[:, TENSOR_COMPONENT_INDEX[tensor_component]]
    else:
        raise ValueError(f"Don't know how to plot '{field_name}' with shape {field.shape}")


def load_xprofile_stress(xprofile_file):
    """Load an xprofile ascii solution and return (distance, stress-by-component)
    for x in [0, 5000], to match the ASPECT transect range."""
    xprofile_data = np.loadtxt(xprofile_file, comments="#")
    xprofile_in_range = xprofile_data[:, 0] <= 5000
    xprofile_distance = xprofile_data[xprofile_in_range, 0]
    xprofile_stress = {
        component: xprofile_data[xprofile_in_range, column]
        for component, column in XPROFILE_STRESS_COLUMN.items()
    }
    return xprofile_distance, xprofile_stress


# increase title sizes
plt.rcParams.update({
    "axes.titlesize": 24,
    "figure.titlesize": 32,
})


fig, axes = plt.subplots(len(TRANSECTS), len(PLOTS), figsize=(5 * len(PLOTS), 5 * len(TRANSECTS)), sharex=True, sharey="col")
fig.suptitle("Elastic chamber stress")

# plotting loop for rows
for row_axes, (row_label, point_a, point_b, reference_file, xprofile_file, xprofile_label) in zip(axes, TRANSECTS):
    line = mesh.sample_over_line(point_a, point_b, resolution=N_SAMPLES - 1)

    # Distance along the transect, starting at 0 at point_a.
    distance = np.linalg.norm(line.points - line.points[0], axis=1)
    reference = pd.read_csv(reference_file) if reference_file is not None else None

    if xprofile_file is not None:
        xprofile_distance, xprofile_stress = load_xprofile_stress(xprofile_file)

    # for loop for each individual plot
    for ax, (field_name, tensor_component, ref_column, xprofile_component, title) in zip(row_axes, PLOTS):
        values = resolve_scalar(line, field_name, tensor_component)
        ax.plot(distance, values, linewidth=2, label="2D ASPECT")
        if xprofile_file is not None:
            # The stone solution uses the opposite sign convention from ASPECT.
            ax.plot(xprofile_distance, -xprofile_stress[xprofile_component], ":", linewidth=4, label=xprofile_label)
        if reference is not None and ref_column is not None:
            # The reference solution uses the opposite sign convention from ASPECT.
            ax.plot(reference["X"], -reference[ref_column], "--",linewidth=2, label="3D Solution")

        ax.set_xlabel("Distance along transect (m)", fontsize = 16)
        ax.set_ylabel(f"{title} (Pa)", fontsize = 16)
        ax.set_title(f"{row_label}: {title}")
        ax.legend(fontsize = 14)
        ax.grid(True)
        ax.axhline(0, color="gray", linewidth=2, zorder=1)

plt.tight_layout()
plt.savefig(OUTPUT_IMAGE, dpi=150)
