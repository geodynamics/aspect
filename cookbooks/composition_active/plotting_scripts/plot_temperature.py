import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

# inputs/config
solution_file = "composition_active/output-composition-active/solution/solution-00079.pvtu"
output_png    = "composition_active/doc/visit0003.png"

# sampling grid resolution and plot tuning
nx         = 1200
n_arrows_x = 70
threshold  = 0.5

def to_2d(sampled, name):
    # pyvista returns flattened Fortran-order arrays; reshape back to (ny, nx) image layout
    return np.asarray(sampled[name]).reshape((nx, ny), order="F").T

# load solution and get domain bounds
mesh = pv.read(solution_file)
xmin, xmax, ymin, ymax = mesh.bounds[:4]
width, height = xmax - xmin, ymax - ymin
ny = int(nx * height / width)  # keep pixels square

# resample the unstructured mesh onto a regular image grid
grid = pv.ImageData(
    dimensions=(nx, ny, 1),
    origin=(xmin, ymin, 0.0),
    spacing=(width / (nx-1), height / (ny-1), 1.0),
)
sampled = grid.sample(mesh)

# mask out points outside the original mesh
valid = to_2d(sampled, "vtkValidPointMask").astype(bool)

def get_field(name):
    return np.where(valid, to_2d(sampled, name), np.nan)

# extract fields needed for the plot
T   = get_field("T")
C1  = get_field("C_1")
vel = np.asarray(sampled["velocity"]).reshape((nx, ny, 3), order="F").transpose(1, 0, 2)
U   = np.where(valid, vel[:, :, 0], np.nan)
V   = np.where(valid, vel[:, :, 1], np.nan)
speed = np.sqrt(U**2 + V**2)

# plotting
fig, ax = plt.subplots(figsize=(16, 16 * height / width), dpi=220)
extent = [xmin, xmax, ymin, ymax]

# temperature background
t_img = ax.imshow(T, origin="lower", extent=extent,
                   cmap=cm.vik, interpolation="bilinear", zorder=1)

# velocity arrows colored by speed
step = max(1, nx // n_arrows_x)
q = ax.quiver(
    np.linspace(xmin, xmax, nx)[::step],
    np.linspace(ymin, ymax, ny)[::step, None],
    U[::step, ::step], V[::step, ::step], speed[::step, ::step],
    cmap=cm.batlow,
    angles="xy", scale_units="xy", scale=32,
    width=0.0012, pivot="mid", zorder=3)

# overlay C_1 = 0.5 contour as the compositional boundary
ax.contour(C1, levels=[threshold], colors="black", linewidths=1.5,
           extent=extent, origin="lower", zorder=4)

# axes formatting
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect="equal")
ax.set_xlabel("X", fontsize=11)
ax.set_ylabel("Y", fontsize=11)
ax.tick_params(direction="out", length=4, width=0.8, labelsize=9)
fig.subplots_adjust(left=0.06, right=0.86, bottom=0.08, top=0.99)

# temperature colorbar (top right)
cb1 = fig.colorbar(t_img, cax=fig.add_axes([0.885, 0.56, 0.018, 0.28]))
cb1.ax.set_title("T", fontsize=9, pad=6)
cb1.ax.tick_params(labelsize=8)

# speed colorbar (bottom right)
cb2 = fig.colorbar(q, cax=fig.add_axes([0.885, 0.20, 0.018, 0.28]))
cb2.ax.set_title("speed", fontsize=9, pad=6)
cb2.ax.tick_params(labelsize=8)

# save final figure
# pick a dpi so the tight-cropped image comes out about 1024 px wide
fig.canvas.draw()
dpi = 1024 / (fig.get_tightbbox().width + 2 * 0.01)
plt.savefig(output_png, dpi=dpi, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)
