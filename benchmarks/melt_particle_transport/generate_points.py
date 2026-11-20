#!/usr/bin/env python3

"""Generate common Gauss quadrature points for the melt particle transport benchmark.

Usage:
    python3 generate_points.py
    python3 generate_points.py --refinement 6
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--refinement", type=int, default=5,
                    help="Refinement level of the integration grid (default: 5)")
parser.add_argument("--gauss", type=int, default=2,
                    help="Number of Gauss points per direction (default: 2)")
args = parser.parse_args()

# Benchmark geometry and integration parameters.
x_extent = 200000
y_extent = 100000
x_repetitions = 2
y_repetitions = 1
gauss_points = args.gauss

# Number and size of cells in the integration grid.
nx = x_repetitions * 2**args.refinement
ny = y_repetitions * 2**args.refinement
hx = x_extent / nx
hy = y_extent / ny

# Gauss quadrature on the reference interval [0,1].
nodes, _ = np.polynomial.legendre.leggauss(gauss_points)
nodes = 0.5 * (nodes + 1.0)

# Generate quadrature points and corresponding integration weights.
points = []

for j in range(ny):
    for i in range(nx):
        for qy in range(gauss_points):
            for qx in range(gauss_points):
                points.append(((i + nodes[qx]) * hx,
                               (j + nodes[qy]) * hy))

points = np.array(points)

# Write the ASPECT parameter snippet.
with open("point_values.prm", "w") as f:
    f.write("subsection Postprocess\n")
    f.write("  subsection Point values\n")
    f.write("    set Evaluation points = ")

    f.write("; \\\n      ".join(
        f"{x:.17g},{y:.17g}" for x, y in points
    ))

    f.write("\n")
    f.write("  end\n")
    f.write("end\n")

print(f"Generated {len(points)} evaluation points "
      f"on a {nx} x {ny} integration grid.")
