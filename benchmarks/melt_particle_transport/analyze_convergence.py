#!/usr/bin/env python3

"""Analyze convergence of FEM and particle advection in melt_particle_transport.

Examples:
    python3 analyze_convergence.py
    python3 analyze_convergence.py --point-refinement 5 --gauss 2
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# Match the geometry and quadrature rule in generate_points.py.
X_EXTENT = 200000.0
Y_EXTENT = 100000.0
X_REPETITIONS = 2
Y_REPETITIONS = 1


def integration_weights(refinement, gauss):
    """Compute quadrature weights in the same order as generate_points.py."""

    nx = X_REPETITIONS * 2**refinement
    ny = Y_REPETITIONS * 2**refinement

    hx = X_EXTENT / nx
    hy = Y_EXTENT / ny

    # Gauss weights on [-1,1].
    _, weights = np.polynomial.legendre.leggauss(gauss)

    # Tensor-product weights, scaled to the physical cell area.
    # The ordering matches generate_points.py:
    # cell y, cell x, quadrature y, quadrature x.
    cell_weights = hx * hy / 4.0 * np.outer(weights, weights).flatten()

    # Repeat for every cell of the integration grid.
    return np.tile(cell_weights, nx * ny)


def read_point_values(filename):
    """Read the five compositional fields at the last output time."""

    data = np.loadtxt(filename, comments="#", ndmin=2)
    available_times = np.unique(data[:, 0])

    # Read the last time that is stored
    time = available_times[-1]

    selected = data[data[:, 0] == time]

    # Layout: time, x, y, other solution components, compositions.
    return float(time), selected[:, 1:3], selected[:, -5:]


def weighted_norm(values, weights):
    """Compute a quadrature approximation of the L2 norm."""

    return float(np.sqrt(np.sum(weights * values**2)))


def analyze(args):

    refinements = sorted(args.refinements)

    if len(refinements) < 3:
        raise ValueError("At least three refinement levels are required.")

    if any(b != a + 1 for a, b in zip(refinements[:-1], refinements[1:])):
        raise ValueError("Refinement levels must be consecutive.")

    # Odd-order Gauss rules include the cell center, which can lie on
    # DG element boundaries if the integration mesh is too coarse.
    if args.point_refinement < refinements[-1] and args.gauss % 2:
        raise ValueError(
            "Odd-order Gauss quadrature on a coarser integration mesh "
            "can place points on DG element boundaries. "
            "Use an even Gauss order or increase --point-refinement."
        )

    # Calculate the contribution of each evaluation point to an
    # integral over the model domain, using the Gauss quadrature rule.
    weights = integration_weights(args.point_refinement, args.gauss)

    # Read the solutions from all refinement levels.
    fields_by_refinement = []
    reference_coords = None

    for refinement in refinements:

        filename = f"output/refinement_{refinement}/point_values.txt"
        time, coordinates, fields = read_point_values(filename)
        fields_by_refinement.append(fields)

        # Make sure all refinements use the same coordinates
        if reference_coords is None:
            reference_coords = coordinates
        else:
            assert np.array_equal(coordinates, reference_coords)

    # Use the porosity from the finest simulation for all error norms.
    porosity = np.clip(fields_by_refinement[-1][:, 0], 0.0, 1.0)

    # Phase-weighted integration weights based on solid/melt fraction.
    # Only compute melt error in places where melt particles are present.
    phase_weights = {
        "solid": weights * (1.0 - porosity),
        "melt": weights * porosity * (porosity > args.melt_cutoff),
    }

    # Common normalization for FEM and particle solutions within each phase.
    average_solid_composition = 0.5 * (fields_by_refinement[-1][:, 1] + fields_by_refinement[-1][:, 2])
    average_melt_composition = 0.5 * (fields_by_refinement[-1][:, 3] + fields_by_refinement[-1][:, 4])
    phase_scales = {
        "solid": weighted_norm(average_solid_composition, phase_weights["solid"]),
        "melt": weighted_norm(average_melt_composition, phase_weights["melt"]),
    }

    results = {}
    csv_rows = []

    # The five compositional fields are assumed to be the last five
    # components of the solution vector, in the following order:
    #
    # 0: porosity
    # 1: solid_tracer_fem
    # 2: solid_tracer_particles
    # 3: melt_tracer_fem
    # 4: melt_tracer_particles

    solution_fields = {
        "solid_fem": (1, "solid"),
        "solid_particles": (2, "solid"),
        "melt_fem": (3, "melt"),
        "melt_particles": (4, "melt"),
    }

    # Analyze each solution field independently.
    for name, (column, phase) in solution_fields.items():

        solutions = [
            fields[:, column]
            for fields in fields_by_refinement
        ]

        scale = phase_scales[phase]

        # Differences between consecutive refinement levels.
        differences = np.array([
            weighted_norm(solutions[i] - solutions[i+1], phase_weights[phase])
            for i in range(len(solutions)-1)
        ])

        # Calculate convergence orders.
        # TODO: simplify
        orders = np.array([
            np.log2(coarse / fine)
            if coarse > 0 and fine > 0 else np.nan
            for coarse, fine in zip(differences[:-1], differences[1:])
        ])

        # Use the finest three refinements to estimate the order.
        p = orders[-1]

        extrapolated = None
        estimated_errors = np.full(len(refinements), np.nan)

        # Richardson extrapolation of the finest solution.
        if np.isfinite(p) and p > 0 and (2.0**p - 1.0) > 1e-8:

            extrapolated = (
                solutions[-1] +
                (solutions[-1] - solutions[-2]) / (2.0**p - 1.0)
            )

            # Relative errors compared to the extrapolated solution.
            estimated_errors = np.array([
                weighted_norm(solution - extrapolated, phase_weights[phase]) / scale
                for solution in solutions
            ])

        results[name] = {
            "errors": estimated_errors,
            "extrapolated": extrapolated,
        }

        print(
            f"{name:17s}: observed orders {orders}; "
            f"finest relative error {estimated_errors[-1]:.5g}"
        )

        # Save the results for each refinement.
        for i, refinement in enumerate(refinements):

            csv_rows.append({
                "field": name,
                "refinement": refinement,

                "cell_size": (
                    X_EXTENT / (X_REPETITIONS * 2**refinement)
                ),

                "estimated_relative_error": estimated_errors[i],
                "observed_order": p,

                "relative_difference_to_next": (
                    differences[i] / scale
                    if i < len(differences)
                    else np.nan
                ),
            })

    # Check whether FEM and particles converge towards the same solution.
    for phase in ("solid", "melt"):

        fem = results[f"{phase}_fem"]["extrapolated"]
        particles = results[f"{phase}_particles"]["extrapolated"]

        if fem is not None and particles is not None:
            discrepancy = weighted_norm(fem - particles, phase_weights[phase])

            print(
                f"{phase}: relative difference between extrapolated "
                f"solutions = {discrepancy / phase_scales[phase]:.5g}"
            )

    # Save the convergence data.
    results_dir = Path("output/results")
    results_dir.mkdir(parents=True, exist_ok=True)

    with (results_dir / "convergence.csv").open(
        "w", newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=csv_rows[0].keys()
        )

        writer.writeheader()
        writer.writerows(csv_rows)

    # Generate convergence plots.
    h = np.array([
        X_EXTENT / (X_REPETITIONS * 2**r)
        for r in refinements
    ])

    # Compare FEM and particle solutions directly at each refinement.
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    for phase, (fem_col, particle_col) in {
        "solid": (1, 2),
        "melt": (3, 4),
    }.items():

        deviations = np.array([
            weighted_norm(
                fields[:, fem_col] - fields[:, particle_col],
                phase_weights[phase]
            ) / phase_scales[phase]
            for fields in fields_by_refinement
        ])

        ax.loglog(h, deviations, "-o", label=phase.capitalize())

        print(f"{phase}: FEM-particle differences = {deviations}")


    ax.loglog(h, deviations[0] * (h / h[0])**2, "--", color="lightgray", label="Second order")
    ax.set_xlabel("Cell size h (m)")
    ax.set_ylabel("Relative FEM-particle L2 difference")
    ax.set_title("Difference between advection methods")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(results_dir / "fem_particle_difference.pdf")
    plt.close(fig)

    # Plot error with respect to extrapolated solution.
    for phase in ("solid", "melt"):

        fig, ax = plt.subplots(figsize=(6.5, 4.5))

        for method, marker in (("fem", "o"), ("particles", "s")):

            errors = results[f"{phase}_{method}"]["errors"]
            valid = np.isfinite(errors) & (errors > 0)

            if np.any(valid):
                ax.loglog(
                    h[valid],
                    errors[valid],
                    "-" + marker,
                    label=method.upper()
                )

        ax.set_xlabel("Cell size h (m)")
        ax.set_ylabel("Relative Richardson-estimated L2 error")
        ax.set_title(f"{phase.capitalize()} tracer")

        ax.grid(True, which="both", alpha=0.3)

        if ax.lines:
            ax.legend()

        fig.tight_layout()

        fig.savefig(
            results_dir / f"convergence_{phase}.pdf"
        )

        plt.close(fig)

    print(f"Results written to {results_dir}")


def main():

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--point-refinement",
        type=int,
        default=5,
        help="Refinement of the common integration grid (default: 5)"
    )

    parser.add_argument(
        "--gauss",
        type=int,
        default=2,
        help="Gauss points per direction per cell (default: 2)"
    )

    parser.add_argument(
        "--refinements",
        type=int,
        nargs="+",
        default=[3, 4, 5, 6],
        help="ASPECT refinement levels (default: 3 4 5 6)"
    )

    parser.add_argument(
        "--melt-cutoff",
        type=float,
        default=1e-6,
        help="Porosity threshold for melt error weighting (default: 1e-6)"
    )

    args = parser.parse_args()

    if args.point_refinement < 0 or args.gauss < 1:
        parser.error(
            "--point-refinement must be >= 0 and --gauss must be >= 1"
        )

    analyze(args)


if __name__ == "__main__":
    main()
