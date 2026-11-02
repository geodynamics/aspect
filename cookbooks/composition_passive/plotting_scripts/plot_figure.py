import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

# inputs/config
solution_file = "composition_passive/output-composition-passive/solution/solution-00178.pvtu"
output_png = "composition_passive/visit0017.png"

# sampling grid resolution and plot tuning
nx = 1200
n_arrows_x = 70
quiver_scale = 32
threshold = 0.5

# warm map for C_1, cool map for C_2
red_cmap = cm.lajolla
blue_cmap = cm.lapaz

# load solution and get domain bounds
mesh = pv.read(solution_file)
xmin, xmax, ymin, ymax = mesh.bounds[:4]
width = xmax - xmin
height = ymax - ymin
ny = int(nx * height / width)  # keep pixels square

# sample the mesh onto a regular grid
grid = pv.ImageData(
    dimensions=(nx, ny, 1),
    origin=(xmin, ymin, 0),
    spacing=(width / (nx - 1), height / (ny - 1), 1),
)
sampled = grid.sample(mesh)

# reshape data onto the grid (pyvista returns flattened Fortran-order arrays)
valid = np.asarray(sampled["vtkValidPointMask"]).reshape((nx, ny), order="F").T.astype(bool)
C1 = np.asarray(sampled["C_1"]).reshape((nx, ny), order="F").T
C2 = np.asarray(sampled["C_2"]).reshape((nx, ny), order="F").T
vel = np.asarray(sampled["velocity"]).reshape((nx, ny, 3), order="F").transpose(1, 0, 2)
U = vel[:, :, 0]
V = vel[:, :, 1]

# mask out points outside the original mesh
C1[~valid] = np.nan
C2[~valid] = np.nan
U[~valid] = np.nan
V[~valid] = np.nan

# color each field directly by its 0..1 value, transparent below the threshold
C1_color = red_cmap(C1)
C2_color = blue_cmap(C2)
C1_color[(C1 < threshold) | np.isnan(C1), 3] = 0
C2_color[(C2 < threshold) | np.isnan(C2), 3] = 0

# plotting
fig, ax = plt.subplots(figsize=(14, 14 * height / width), dpi=220)
extent = [xmin, xmax, ymin, ymax]

# draw compositions (C2 first so C1 sits on top where they overlap)
ax.imshow(C2_color, origin="lower", extent=extent, interpolation="bilinear")
ax.imshow(C1_color, origin="lower", extent=extent, interpolation="bilinear")

# uncomment to draw the element edges on top (used for the visit0017 mesh figure)
# from matplotlib.collections import LineCollection
# edges = mesh.extract_all_edges()
# pairs = edges.lines.reshape(-1, 3)[:, 1:]
# segs  = edges.points[pairs][:, :, :2]
# ax.add_collection(LineCollection(segs, colors="black", linewidths=0.25))

# draw velocity arrows
step = max(1, nx // n_arrows_x)
x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
ax.quiver(
    x[::step],
    y[::step, None],
    U[::step, ::step],
    V[::step, ::step],
    color="black",
    angles="xy",
    scale_units="xy",
    scale=quiver_scale,
    width=0.0012,
    pivot="mid",
)

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect("equal")

# axis labels and ticks
x_ticks = np.linspace(xmin, xmax, 5)
y_ticks = np.linspace(ymin, ymax, 5)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels([f"{v:g}" for v in x_ticks])
ax.set_yticklabels([f"{v:g}" for v in y_ticks])
ax.set_xlabel("X")
ax.set_ylabel("Y")

# make space for the two small colorbars
fig.subplots_adjust(left=0.055, right=0.86, bottom=0.085, top=0.99)

# colorbar for C_1 (top right)
cax1 = fig.add_axes([0.885, 0.56, 0.018, 0.28])
cb1 = fig.colorbar(plt.cm.ScalarMappable(cmap=red_cmap), cax=cax1)
cb1.ax.set_title("C_1", fontsize=9, pad=6)
cb1.ax.tick_params(labelsize=8)

# colorbar for C_2 (bottom right)
cax2 = fig.add_axes([0.885, 0.20, 0.018, 0.28])
cb2 = fig.colorbar(plt.cm.ScalarMappable(cmap=blue_cmap), cax=cax2)
cb2.ax.set_title("C_2", fontsize=9, pad=6)
cb2.ax.tick_params(labelsize=8)

# save final figure
# pick a dpi so the tight-cropped image comes out about 1024 px wide
fig.canvas.draw()
dpi = 1024 / (fig.get_tightbbox().width + 2 * 0.01)
plt.savefig(output_png, dpi=dpi, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)
