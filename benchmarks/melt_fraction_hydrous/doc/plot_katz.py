import numpy as np
import matplotlib.pyplot as plt
import glob
import os

KATZ_DIR = "katz_original_code/meltParam/src"

cases = [
    ("0 wt%", "0ppm", os.path.join(KATZ_DIR, "katz_original_M017_0wt.dat")),
    ("0.02 wt%", "200ppm", os.path.join(KATZ_DIR, "katz_original_M017_0p02wt.dat")),
    ("0.05 wt%", "500ppm", os.path.join(KATZ_DIR, "katz_original_M017_0p05wt.dat")),
    ("0.10 wt%", "1000ppm", os.path.join(KATZ_DIR, "katz_original_M017_0p10wt.dat")),
    ("0.30 wt%", "3000ppm", os.path.join(KATZ_DIR, "katz_original_M017_0p30wt.dat")),
]

resolutions = ["Res7", "Res9"]

target_P = 1e9

color_original = [
    "lightskyblue",
    "moccasin",
    "lightgreen",
    "lightcoral",
    "plum",
]

color_aspect = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
]

for resolution in resolutions:
    plt.figure(figsize=(8, 6))

    for i, (name, ppm, original_file) in enumerate(cases):

        aspect_dir = f"benchmark_katz_{ppm}_{resolution}"

        # Original Katz implementation
        original = np.loadtxt(original_file)
        T_original = original[:,0]
        F_original = original[:,1]

        # ASPECT implementation
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

        P = data[:,4]
        T_C = data[:,5] - 273.15
        F_aspect = data[:,6]

        P_unique = np.unique(P)
        P_selected = P_unique[np.argmin(np.abs(P_unique - target_P))]

        mask = P == P_selected

        T = T_C[mask]
        aspect = F_aspect[mask]

        order = np.argsort(T)

        plt.plot(T_original, F_original, color=color_original[i], label=f"Katz original {name}")
        plt.plot(T[order], aspect[order], "--", color=color_aspect[i], label=f"ASPECT {name}")

    plt.xlabel("Temperature (°C)")
    plt.ylabel("Melt fraction F")
    plt.xlim(900, 1400)
    plt.ylim(0, 0.4)
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"Katz_1GPa_comparison_{resolution}.png", dpi=300)

    plt.close()