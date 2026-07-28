import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

# inputs/config
solution_file = "composition_active/output-composition-active/solution/solution-00000.pvtu"
output_png    = "composition_active/doc/visit0007.png"
comp_name     = "C_1"

# domain bounds (model coordinates)
x0, x1 = 0.0, 2.0
y0, y1 = 0.0, 1.0

# sampling grid resolution and plot tuning
nx           = 1600
n_arrows_x   = 70
quiver_scale = 25

def reshape_data(sampled, name):
    # pyvista returns flattened Fortran-order arrays; reshape back to (ny, nx) image layout
    return np.asarray(sampled[name]).reshape((nx, ny), order="F").T

# load solution and trim to the region of interest
mesh = pv.read(solution_file)
mesh = mesh.clip_box(bounds=[x0, x1, y0, y1, 0, 0], invert=False)

width, height = x1 - x0, y1 - y0
ny = int(nx * height / width)  # keep pixels square

# resample the unstructured mesh onto a regular image grid
grid = pv.ImageData(
    dimensions=(nx, ny, 1),
    origin=(x0, y0, 0.0),
    spacing=(width / (nx-1), height / (ny-1), 1.0),
)
sampled = grid.sample(mesh)

# coordinate arrays for plotting/quiver
x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
X, Y = np.meshgrid(x, y)

# extract composition field, validity mask, and velocity components
C1    = reshape_data(sampled, comp_name)
valid = reshape_data(sampled, "vtkValidPointMask").astype(bool)
vel   = np.asarray(sampled["velocity"]).reshape((nx, ny, 3), order="F").transpose(1, 0, 2)

# mask out points outside the original mesh
C1 = np.where(valid, C1, np.nan)
U  = np.where(valid, vel[:, :, 0], np.nan)
V  = np.where(valid, vel[:, :, 1], np.nan)

# plotting
fig, ax = plt.subplots(figsize=(16, 16 * height / width), dpi=220)
ax.set_facecolor("white")

# composition field as an image
im = ax.imshow(C1, origin="lower", extent=[x0, x1, y0, y1], cmap=cm.lapaz,
               vmin=0, vmax=1, interpolation="bilinear", zorder=1)

# velocity arrows, subsampled on a regular grid so they're not too dense
step = max(1, nx // n_arrows_x)
ax.quiver(
    x[::step], y[::step, None],
    U[::step, ::step], V[::step, ::step],
    color="black", angles="xy", scale_units="xy",
    scale=quiver_scale, width=0.0010, pivot="mid", zorder=2)

# axes formatting
ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect="equal")
ax.set_xlabel("X", fontsize=11)
ax.set_ylabel("Y", fontsize=11)
ax.tick_params(direction="out", length=4, width=0.8, labelsize=9)
fig.subplots_adjust(left=0.06, right=0.91, bottom=0.10, top=0.99)

# colorbar on the right side
cax = fig.add_axes([0.925, 0.25, 0.015, 0.50])
cb  = fig.colorbar(im, cax=cax, ticks=np.linspace(0, 1, 6))
cb.ax.set_title("C_1", fontsize=9, pad=6)
cb.ax.tick_params(labelsize=8)

# save final figure
# pick a dpi so the tight-cropped image comes out about 1024 px wide
fig.canvas.draw()
dpi = 1024 / (fig.get_tightbbox().width + 2 * 0.01)
plt.savefig(output_png, dpi=dpi, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)
