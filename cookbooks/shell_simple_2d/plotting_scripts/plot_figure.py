import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

# tick values in scientific notation
def sci_split(vals):
    m = max(abs(min(vals)), abs(max(vals)))
    e = int(np.floor(np.log10(m))) if m > 0 else 0
    if -2 < e < 2:
        return [f"{v:g}" for v in vals], ""
    return [f"{v / 10**e:g}" for v in vals], f" ×10$^{{{e}}}$"

pv.set_plot_theme("document")

# inputs/config
solution_file = "shell_simple_2d/output-shell_simple_2d/solution/solution-01000.pvtu"
output_png    = "shell_simple_2d/doc/x-movie1000.png"

# show the mesh over the lower part of the annulus only
mesh_max_angle = 35   # degrees

# load the solution
mesh = pv.read(solution_file)
xmin, xmax, ymin, ymax = mesh.bounds[:4]
size = max(xmax - xmin, ymax - ymin)
T = np.asarray(mesh["T"], dtype=float)

# pick the cells in the lower sector of the annulus for the mesh overlay
centers = mesh.cell_centers().points
angle = np.degrees(np.arctan2(centers[:, 1], centers[:, 0]))
sector = mesh.extract_cells(angle < mesh_max_angle)

# pyvista renders the scene only; axes and scalar bar are added with matplotlib below
plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")

# temperature field
plotter.add_mesh(mesh, scalars="T", cmap=cm.lajolla, lighting=False, show_scalar_bar=False)

# wireframe over the lower sector
plotter.add_mesh(sector, style="wireframe", color="ivory_black", line_width=0.3, opacity=0.1)

# straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
plotter.camera_position = "xy"
plotter.enable_parallel_projection()
plotter.camera.parallel_scale = size / 2

img = plotter.screenshot(return_img=True)
plotter.close()

# plotting
fig, ax = plt.subplots(figsize=(9, 9))
ax.imshow(img, extent=[xmin / 1e3, xmax / 1e3, ymin / 1e3, ymax / 1e3])
ax.set_aspect("equal")

# axis ticks and labels in km
x_ticks = np.linspace(xmin / 1e3, xmax / 1e3, 5)
y_ticks = np.linspace(ymin / 1e3, ymax / 1e3, 5)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
x_labels, x_suf = sci_split(x_ticks)
ax.set_xticklabels(x_labels)
y_labels, y_suf = sci_split(y_ticks)
ax.set_yticklabels(y_labels)
ax.set_xlabel(f"X (km{x_suf})")
ax.set_ylabel(f"Y (km{y_suf})")
ax.tick_params(direction="out", length=4, width=0.8, labelsize=9)

fig.subplots_adjust(left=0.10, right=0.88, bottom=0.08, top=0.99)

# scalar bar
cax = fig.add_axes([0.90, 0.25, 0.018, 0.50])
sm = plt.cm.ScalarMappable(cmap=cm.lajolla)
sm.set_clim(T.min(), T.max())
cb_ticks = np.linspace(T.min(), T.max(), 6)
cb_labels, cb_suf = sci_split(cb_ticks)
cb = fig.colorbar(sm, cax=cax, ticks=cb_ticks)
cb.ax.set_title(f"T (K{cb_suf})", fontsize=9, pad=6, loc="left")
cb.ax.set_yticklabels(cb_labels)
cb.ax.tick_params(labelsize=8)

plt.savefig(output_png, dpi=100, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)