import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

# inputs/config
solution_file = "composition-reaction/output-composition-reaction/solution/solution-00000.pvtu"
output_png = "composition-reaction/doc/0.png"

c1_name = "C_1"
c2_name = "C_2"
temp_name = "T"

# domain
x0 = 0.0
x1 = 2.0
y0 = 0.0
y1 = 1.0

# grid resolution
nx = 1600
n_levels = 160

# turns long sampled data into a 2d grid
def reshape_data(sampled, name):
    return np.asarray(sampled[name]).reshape((nx, ny), order="F").T

# load solution and trim to the region of interest
mesh = pv.read(solution_file)
mesh = mesh.clip_box(bounds=[x0, x1, y0, y1, 0, 0], invert=False)

width = x1 - x0
height = y1 - y0
ny = int(nx * height / width)  # keep pixels square

# sample onto a regular grid
grid = pv.ImageData()
grid.dimensions = (nx, ny, 1)
grid.origin = (x0, y0, 0.0)
grid.spacing = (width / (nx - 1), height / (ny - 1), 1.0)
sampled = grid.sample(mesh)

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
X, Y = np.meshgrid(x, y)

# fields
T = reshape_data(sampled, temp_name)
C1 = reshape_data(sampled, c1_name)
C2 = reshape_data(sampled, c2_name)
valid = reshape_data(sampled, "vtkValidPointMask").astype(bool)

# remove bad points
T = np.where(valid, T, np.nan)
C1 = np.where(valid, C1, np.nan)
C2 = np.where(valid, C2, np.nan)

# plotting
fig_w = 14
fig_h = fig_w * height / width
fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=220)

# temperature field
levels = np.linspace(np.nanmin(T), np.nanmax(T), n_levels)
cf = ax.contourf(
    X,
    Y,
    T,
    levels=levels,
    cmap=cm.vik,
    zorder=1,
)

# black line: material from the bottom
ax.contour(
    X,
    Y,
    C1,
    levels=[0.5],
    colors="black",
    linewidths=1.3,
    zorder=2,
)

# white line: material made by reaction
ax.contour(
    X,
    Y,
    C2,
    levels=[0.5],
    colors="white",
    linewidths=1.1,
    zorder=3,
)

ax.set_xlim(x0, x1)
ax.set_ylim(y0, y1)
ax.set_aspect("equal")

# axis ticks and labels
x_ticks = np.linspace(x0, x1, 5)
y_ticks = np.linspace(y0, y1, 5)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels([f"{v:g}" for v in x_ticks], fontsize=9)
ax.set_yticklabels([f"{v:g}" for v in y_ticks], fontsize=9)
ax.set_xlabel("X", fontsize=11)
ax.set_ylabel("Y", fontsize=11)
ax.tick_params(direction="out", length=4, width=0.8)

# leave room for scalar bar
fig.subplots_adjust(left=0.06, right=0.90, bottom=0.10, top=0.99)

# scalar bar for temperature
ticks = np.linspace(np.nanmin(T), np.nanmax(T), 6)
cax = fig.add_axes([0.92, 0.25, 0.015, 0.50])
cb = fig.colorbar(cf, cax=cax, ticks=ticks)
cb.ax.set_title("T", fontsize=9, pad=6)
cb.ax.set_yticklabels([f"{v:g}" for v in ticks])
cb.ax.tick_params(labelsize=8)

# save final figure
# pick a dpi so the tight-cropped image comes out about 1024 px wide
fig.canvas.draw()
dpi = 1024 / (fig.get_tightbbox().width + 2 * 0.01)
plt.savefig(output_png, dpi=dpi, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)
