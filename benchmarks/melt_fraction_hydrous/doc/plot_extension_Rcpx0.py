import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import os

aspect_dir = "benchmark_ball_depleted"

r0 = 1.2471514
r1 = -0.17266209

# Pressure where R_cpx = r0 + r1 * P = 0
P_rcpx0 = -r0 / r1

# Representative temperatures for different melting values.
targets_C = [1700, 1750, 1800, 1900]

rows = []

for fname in glob.glob(os.path.join(aspect_dir, "solution", "*.gnuplot")):
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) == 7:
                rows.append([float(x) for x in parts])

data = np.array(rows)

P_GPa = data[:, 4] / 1e9
T_C = data[:, 5] - 273.15
F = data[:, 6]

# Pressure range around R_cpx = 0
mask = (P_GPa >= 6.5) & (P_GPa <= 8.0)

P_GPa = P_GPa[mask]
T_C = T_C[mask]
F = F[mask]

T_unique = np.unique(np.round(T_C, 2))

plt.figure(figsize=(8, 5))

print(f"Rcpx = 0 at P = {P_rcpx0:.3f} GPa")

for T_target in targets_C:

    T_sel = T_unique[np.argmin(np.abs(T_unique - T_target))]

    m = np.isclose(np.round(T_C, 2), T_sel)

    P_line = P_GPa[m]
    F_line = F[m]

    order = np.argsort(P_line)

    P_line = P_line[order]
    F_line = F_line[order]

    idx_rcpx0 = np.argmin(np.abs(P_line - P_rcpx0))
    print(f"Target T = {T_target:.0f} °C -> "f"using {T_sel:.2f} °C, "f"F near Rcpx=0 = {F_line[idx_rcpx0]:.5f}")

    plt.plot(P_line, F_line, marker="o", markevery=20, markersize=4, label=f"{T_sel:.2f} °C")

plt.axvline(P_rcpx0, linestyle="--", linewidth=1.5, label=r"$R_{cpx}=0$")

plt.xlabel("Pressure (GPa)")
plt.ylabel("Melt fraction F")
plt.xlim(6.5, 8.0)
plt.ylim(bottom=0)
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig("extension_Rcpx0_depleted.png", dpi=300)

plt.close()