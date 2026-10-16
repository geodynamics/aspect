import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob
import os

sys.path.insert(0, "BDD21/src/Ball_Duvernay_Davies_2022")
from Melt import Katz

model = Katz()

cases = [
    ("Primitive", "primitive", 8.0, {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.17127832,
        "r0": 0.99133219,
        "r1": -0.12357699,
        "X_H2O_bulk": 0.028,
    }),

    ("50/50", "5050", 7.64, {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.1702248,
        "r0": 1.09774656,
        "r1": -0.14365651,
        "X_H2O_bulk": 0.02,
    }),

    ("Depleted", "depleted", 7.22, {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.16783327,
        "r0": 1.2471514,
        "r1": -0.17266209,
        "X_H2O_bulk": 0.01,
    }),
]

for name, case_id, P_max, mantle_params in cases:

    aspect_dir = f"benchmark_ball_{case_id}"

    # Read ASPECT output
    rows = []

    for fname in glob.glob(
        os.path.join(aspect_dir, "solution", "*.gnuplot")
    ):
        with open(fname) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.split()

                if len(parts) == 7:
                    rows.append([float(x) for x in parts])

    data = np.array(rows)

    P_pa = data[:, 4]
    T_K = data[:, 5]
    F_aspect = data[:, 6]

    P_GPa = P_pa / 1e9

    mask = (P_pa >= 0) & (P_GPa < P_max)

    P_pa = P_pa[mask]
    T_K = T_K[mask]
    F_aspect = F_aspect[mask]

    # Original Ball implementation
    F_ball = []

    for P, T in zip(P_pa, T_K):
        F = model.KatzPT(P/1e9, T, inputConst=mantle_params)
        F_ball.append(F)

    F_ball = np.array(F_ball)

    deltaF = F_aspect - F_ball

    # P-T plot
    P_GPa = P_pa/1e9
    T_C = T_K - 273.15

    plt.figure(figsize=(10, 8))
    sc = plt.scatter(P_GPa, T_C, c=deltaF, cmap="bwr", vmin=-2e-4, vmax=2e-4, marker="s", s=10)

    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Temperature (°C)")
    plt.xlim(0, 8)
    plt.ylim(1100, 2000)

    cbar = plt.colorbar(sc)
    cbar.ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    cbar.set_label("Delta F = F_ASPECT - F_Ball")

    plt.tight_layout()
    plt.savefig(f"Ball_{case_id}_PT_comparison.png", dpi=300)

    plt.close()