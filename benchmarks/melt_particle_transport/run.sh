#!/bin/bash
set -euo pipefail

# Number of MPI processes and ASPECT executable.
NP=${NP:-8}
ASPECT=${ASPECT:-./aspect-release}

# Create a common grid for all refinements we can use
# to evaluate the solution using ASPECT's point value postprocessor
# This is very slow for a large number of points (refinements above 6).
POINT_REFINEMENT=6
GAUSS_POINTS=2
python3 generate_points.py --refinement "$POINT_REFINEMENT" --gauss "$GAUSS_POINTS"

# Parameter file used for all simulations.
base_prm="melt_particle_transport.prm"

# Run the model at different global refinement levels.
for r in 3 4 5 6 7 8
do
    echo "Running refinement $r"

    prm_file="melt_particle_transport_r${r}.prm"
    output_dir="output/refinement_${r}"

    # Remove results from previous runs of this refinement.
    rm -rf -- "$output_dir"
    mkdir -p "$output_dir"

    # Copy the original parameter file.
    cp "$base_prm" "$prm_file"

    # Specify the refinement level and output directory.
    echo "set Output directory = $output_dir"   >> "$prm_file"
    echo ""                                     >> "$prm_file"
    echo "subsection Mesh refinement"           >> "$prm_file"
    echo "  set Initial global refinement = $r" >> "$prm_file"
    echo "end"                                  >> "$prm_file"

    cat point_values.prm >> "$prm_file"

    mkdir -p "$output_dir"

    # Run ASPECT
    mpirun -np "$NP" "$ASPECT" "$prm_file"
    rm "$prm_file"
done

rm point_values.prm

echo "All simulations completed."

python3 analyze_convergence.py \
    --refinements 3 4 5 6 7 8 \
    --point-refinement "$POINT_REFINEMENT" \
    --gauss "$GAUSS_POINTS"

echo "Analysis completed."
