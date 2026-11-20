#!/usr/bin/env python3
"""
This script compares the original ASPECT hydrous melt fraction postprocessor results
with the extension of the water saturation limit.

For a fair comparison, the two Gnuplot files must come from the same benchmark setup.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Files relative to the /aspect/build: 
ORIGINAL = Path("../benchmarks/melt_fraction_hydrous/original_hydrous_melt_fraction.gnuplot")
EXTENDED = Path("../tests/melt_fraction_postprocessor_hydrous/solution/solution-00000.0000.gnuplot")
OUTPUT = Path("water_sat_validation")


# Primitive-mantle parameters (must agree with the benchmark .prm):
M_CPX = 0.17127832
R1 = 0.99133219
R2 = -1.2357699e-10
BULK_PPM = 280.0
D_H2O = 0.01
K_H2O = 43.0
GAMMA = 0.75
K1, K2, LAMBDA = 12.0, 1.0, 0.6
A1, B1, BETA1 = 1085.7, 1520.0, 1.5

#Read rows x, y, vx, vy, pressure (Pa), temperature (K), F:
def read_gnuplot(path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row = tuple(map(float, line.split()))
        if len(row) != 7:
            raise ValueError(f"Expected seven columns in {path}")
        rows.append(row)
    if not rows:
        raise ValueError(f"No data in {path}")
    return rows

#Average surface vertices with the same temperature:
def surface_curve(rows, y_top):
    values = {}
    for row in rows:
        if row[1] == y_top:
            values.setdefault(row[5], []).append(row[6])
    return {t: sum(fractions) / len(fractions) for t, fractions in values.items()}


# Global comparison:
original = read_gnuplot(ORIGINAL)
extended = read_gnuplot(EXTENDED)
if len(original) != len(extended):
    raise ValueError("Gnuplot files have different row counts")
if any(old[:6] != new[:6] for old, new in zip(original, extended)):
    raise ValueError("Gnuplot files have different coordinates, pressures or temperatures")

changed = [(old, new) for old, new in zip(original, extended) if old[6] != new[6]]
decreased = sum(new[6] < old[6] for old, new in changed)
increased = sum(new[6] > old[6] for old, new in changed)
maximum_reduction = max([0.0] + [old[6] - new[6] for old, new in changed])
old_opx = oversaturated = opx_to_cpx = 0

for old, new in changed:
    pressure = max(0.0, old[4])  # Pa.
    pressure_gpa = pressure / 1e9
    r_cpx = R1 + R2 * pressure
    f_cpx_out = 1.0 if r_cpx <= 0 else min(1.0, M_CPX / r_cpx)

    originally_in_opx = old[6] > f_cpx_out
    old_opx += originally_in_opx
    opx_to_cpx += originally_in_opx and new[6] <= f_cpx_out

    x_water = (BULK_PPM / 10000.0) / (D_H2O + old[6] * (1.0 - D_H2O))
    x_saturation = K1 * pressure_gpa**LAMBDA + K2 * pressure_gpa
    oversaturated += x_water > x_saturation

print(f"Gnuplot rows compared: {len(original)}")
print(f"Unchanged rows: {len(original) - len(changed)}")
print(f"Rows with decreased F: {decreased}")
print(f"Rows with increased F: {increased}")
print(f"Maximum reduction in F: {maximum_reduction:.6f}")
print(f"Changed rows previously in opx: {old_opx}")
print(f"Changed rows with calculated oversaturation: {oversaturated}")
print(f"Changed rows crossing opx to cpx: {opx_to_cpx}")
print(f"Changed y coordinates: {sorted({old[1] for old, _ in changed})}")
if changed:
    changed_pressures = [max(0.0, old[4]) / 1e9 for old, _ in changed]
    print(f"Changed pressure range: {min(changed_pressures):.6g} to {max(changed_pressures):.6g} GPa")
print("Note: Gnuplot rows may repeat mesh vertices.")


# Surface curves and reference temperatures at P=0:
y_top = max(row[1] for row in original)
old_surface = surface_curve(original, y_top)
new_surface = surface_curve(extended, y_top)
if old_surface.keys() != new_surface.keys():
    raise ValueError("The surface temperature grids differ")
temperatures = sorted(old_surface)

f_cpx_out_zero = 1.0 if R1 <= 0 else min(1.0, M_CPX / R1)
t_cpx_out_zero = A1 + 273.15 + f_cpx_out_zero**(1 / BETA1) * (B1 - A1)
x_water_at_cpx_out = (BULK_PPM / 10000.0) / (D_H2O + f_cpx_out_zero * (1 - D_H2O))
t_original_switch = t_cpx_out_zero - K_H2O * x_water_at_cpx_out**GAMMA


# Save PNG plots:
OUTPUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 4.8))
ax.plot(temperatures, [old_surface[t] for t in temperatures], "o-", markersize=3, color="#4C6A92", label="Original ASPECT (based on Ball et al. (2022))")
ax.plot(temperatures, [new_surface[t] for t in temperatures], "o-", markersize=3, color="#B85C5C", label="Extended water saturation")
ax.axvline(t_original_switch, color="#4C6A92", linestyle="--", linewidth=1.4, label=f"Original cpx/opx boundary: {t_original_switch:.2f} K")
ax.axvline(t_cpx_out_zero, color="#B85C5C", linestyle="--", linewidth=1.4, label=f"Extended cpx/opx boundary: {t_cpx_out_zero:.2f} K")
ax.set(xlabel="Temperature (K)", ylabel="Melt fraction F", title=f"Surface comparison (y = {y_top:g})")
ax.grid(alpha=0.25)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTPUT / "surface_melt_fraction_comparison.png", dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.8))
ax.plot(temperatures, [new_surface[t] - old_surface[t] for t in temperatures], "o-", color="#4C6A92", markersize=3)
ax.axhline(0, linewidth=0.8)
ax.set(xlabel="Temperature (K)", ylabel="F(extended) - F(original)", title=f"Surface melt fraction difference (y = {y_top:g})")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUTPUT / "surface_melt_fraction_difference.png", dpi=180)
plt.close(fig)

print(f"Original regime boundary: {t_original_switch:.2f} K")
print(f"Extended implementation regime boundary: {t_cpx_out_zero:.2f} K")
print(f"Plots saved in {OUTPUT}/")