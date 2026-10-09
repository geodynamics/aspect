import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

sys.path.insert(0, "BDD21/src/Ball_Duvernay_Davies_2022")
from Melt import Katz

model = Katz()

cases = [
    ("Primitive", "primitive", {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.17127832,
        "r0": 0.99133219,
        "r1": -0.12357699,
        "X_H2O_bulk": 0.028,
    }),

    ("50/50", "5050", {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.1702248,
        "r0": 1.09774656,
        "r1": -0.14365651,
        "X_H2O_bulk": 0.02,
    }),

    ("Depleted", "depleted", {
        "B1": 1520 + 273.15,
        "beta2": 1.2,
        "M_cpx": 0.16783327,
        "r0": 1.2471514,
        "r1": -0.17266209,
        "X_H2O_bulk": 0.01,
    }),
]

target_P = 1e9

colors_original = [
    "lightskyblue",
    "moccasin",
    "lightgreen",
]

colors_aspect = [
    "tab:blue",
    "tab:orange",
    "tab:green",
]

plt.figure(figsize=(8, 6))

for i, (name, case_id, mantle_params) in enumerate(cases):

    aspect_dir = f"benchmark_ball_{case_id}"

    # ASPECT implementation
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

    P = data[:, 4]
    T_C = data[:, 5] - 273.15
    F_aspect = data[:, 6]

    P_unique = np.unique(P)
    P_selected = P_unique[np.argmin(np.abs(P_unique - target_P))]

    mask = P == P_selected

    T = T_C[mask]
    aspect = F_aspect[mask]

    # Original Ball implementation
    original = []

    for T_C_value in T:
        F_original = model.KatzPT(P_selected / 1e9,T_C_value + 273.15, inputConst=mantle_params)
        original.append(F_original)

    original = np.array(original)

    order = np.argsort(T)

    plt.plot(T[order], original[order], color=colors_original[i], label=f"Ball original {name}")
    plt.plot(T[order], aspect[order], "--", color=colors_aspect[i], label=f"ASPECT {name}")

plt.xlabel("Temperature (°C)")
plt.ylabel("Melt fraction F")
plt.legend()
plt.tight_layout()

plt.savefig("Ball_1GPa_comparison.png", dpi=300)
plt.close()