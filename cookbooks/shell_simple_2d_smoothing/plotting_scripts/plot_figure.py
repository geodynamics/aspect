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
solution_file = "shell_simple_2d_smoothing/output-shell_simple_2d_smoothing/solution/solution-00185.pvtu"
output_png    = "shell_simple_2d_smoothing/doc/smoothing.png"
warp_factor   = 6e6   # m of lift per unit of C_1 — keeps the bounding box roughly cubic

# optional interactive html file output (uncomment the export html line at the end too)
# html_file = "shell_simple_2d_smoothing/doc/with_smoothing.html"

# lift the annulus surface upward by the C_1 value
mesh = pv.read(solution_file)
warped = mesh.warp_by_scalar(scalars="C_1", factor=warp_factor)

# pyvista renders the scene only; axes and scalar bar are added with matplotlib below
plotter = pv.Plotter(off_screen=True)
plotter.window_size = (1024, 1024)
plotter.set_background("white")
plotter.add_mesh(warped, scalars="C_1", cmap=cm.lapaz, clim=[0, 1],
                 smooth_shading=True, lighting=False, show_scalar_bar=False)

# camera looking down from the north-east
xmin, xmax, ymin, ymax = mesh.bounds[:4]
scale = max(xmax - xmin, ymax - ymin)
cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
plotter.camera_position = [
    (cx + 1.5 * scale, cy + 1.5 * scale, warp_factor + 0.8 * scale),
    (cx, cy, warp_factor / 2),
    (0, 0, 1),
]
plotter.camera.zoom(0.8)

img = plotter.screenshot(return_img=True)
H = img.shape[0]

# turns a 3d point into a pixel position in the screenshot
def to_px(p):
    plotter.renderer.SetWorldPoint(p[0], p[1], p[2], 1.0)
    plotter.renderer.WorldToDisplay()
    x, y, _ = plotter.renderer.GetDisplayPoint()
    return x, H - y

# axis lines along the domain edges
x_ticks = np.linspace(0, xmax, 5)
y_ticks = np.linspace(0, ymax, 5)
z_ticks = np.linspace(0, warp_factor, 4)
tick = 0.025 * scale   # tick length in meters
x_labels, x_suf = sci_split(x_ticks / 1e3)
y_labels, y_suf = sci_split(y_ticks / 1e3)
z_labels, z_suf = sci_split(z_ticks / 1e3)

# axes go on the box edges facing the camera so they never cross the surface
axis_lines = [
    (to_px((0, ymax, 0)), to_px((xmax, ymax, 0))),
    (to_px((xmax, 0, 0)), to_px((xmax, ymax, 0))),
    (to_px((xmax, 0, 0)), to_px((xmax, 0, warp_factor))),
]
x_tk = [(to_px((t, ymax, 0)), to_px((t, ymax + tick, 0)), to_px((t, ymax + 3 * tick, 0))) for t in x_ticks]
y_tk = [(to_px((xmax, t, 0)), to_px((xmax + tick, t, 0)), to_px((xmax + 3 * tick, t, 0))) for t in y_ticks]
z_tk = [(to_px((xmax, 0, t)), to_px((xmax + 0.7 * tick, -0.7 * tick, t)),
         to_px((xmax + 2.2 * tick, -2.2 * tick, t))) for t in z_ticks]
titles = [
    (to_px((xmax / 2, ymax + 6.5 * tick, 0)), f"X (km{x_suf})"),
    (to_px((xmax + 6.5 * tick, ymax / 2, 0)), f"Y (km{y_suf})"),
    (to_px((xmax + 6 * tick, -6 * tick, warp_factor / 2)), f"Z (km{z_suf})"),
]

plotter.show_bounds(xtitle="X (m)", ytitle="Y (m)", ztitle="Z (m)")
plotter.add_scalar_bar(title="C_1")
plotter.camera.zoom(0.75)  # pull back so the axes fit in the html view
plotter.render()

# optional interactive html file (needs pip install "pyvista[jupyter]")
# plotter.export_html(html_file)
plotter.close()

# plotting
fig, ax = plt.subplots(figsize=(9, 9))
ax.imshow(img)
ax.set_axis_off()

for p0, p1 in axis_lines:
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="black", linewidth=0.8)

for ticks_px, tick_labels in [(x_tk, x_labels), (y_tk, y_labels), (z_tk, z_labels)]:
    for (base, tip, anchor), lab in zip(ticks_px, tick_labels):
        ax.plot([base[0], tip[0]], [base[1], tip[1]], color="black", linewidth=0.8)
        ax.text(anchor[0], anchor[1], lab, fontsize=9, ha="center", va="center")

for (px_, py_), label in titles:
    ax.text(px_, py_, label, fontsize=11, ha="center", va="center")

# crop to the drawn content so there is as little whitespace as possible
pts = np.array([p for line in axis_lines for p in line]
               + [p for tk in (x_tk + y_tk + z_tk) for p in tk]
               + [p for p, _ in titles])
rows, cols = np.where(img[:, :, :3].min(axis=2) < 250)
left   = min(pts[:, 0].min(), cols.min()) - 15
right  = max(pts[:, 0].max(), cols.max()) + 15
top    = min(pts[:, 1].min(), rows.min()) - 15
bottom = max(pts[:, 1].max(), rows.max()) + 15
ax.set_xlim(left, right)
ax.set_ylim(bottom, top)

fig.subplots_adjust(left=0.0, right=0.90, bottom=0.02, top=0.98)

# scalar bar
cax = fig.add_axes([0.92, 0.25, 0.018, 0.50])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm.lapaz), cax=cax,
                  ticks=[0, 0.25, 0.5, 0.75, 1])
cb.ax.set_title("C_1", fontsize=9, pad=6)
cb.ax.tick_params(labelsize=8)

plt.savefig(output_png, dpi=100, bbox_inches="tight", pad_inches=0.01, facecolor="white")
plt.close(fig)

# free the vtk objects before python shuts down to prevent errors on exit
del plotter, warped, mesh
