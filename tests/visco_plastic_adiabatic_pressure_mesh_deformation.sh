#!/bin/bash

if [ "$1" = "screen-output" ]; then
  cat > /dev/null
  echo "Current-surface points (x, y, geometry depth, current depth, pressure used for plasticity in MPa, diffusion viscosity in Pa s):"
  # Invert the 2D Drucker-Prager yield-stress equation to recover the pressure
  # used by the material model: yield_stress = C cos(phi) + p sin(phi).
  awk '
    BEGIN {
      pi = atan2(0, -1)
      cohesion = 1e6
      friction = 10 * pi / 180
    }
    !/^#/ && $14 == 0 {
      pressure = ($9 - cohesion * cos(friction)) / sin(friction)
      if (pressure > -500 && pressure < 500)
        pressure = 0
      printf "%g %g %.0f %.0f %.3f %.6e\n", $1, $2, $13, $14, pressure / 1e6, $11
    }' \
    output-visco_plastic_adiabatic_pressure_mesh_deformation/solution/solution-00001.0000.gnuplot \
    | sort -n -k1,1 \
    | uniq
else
  cat
fi
