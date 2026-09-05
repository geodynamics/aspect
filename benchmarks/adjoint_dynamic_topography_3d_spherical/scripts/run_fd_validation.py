#!/usr/bin/env python3
"""Run a cellwise finite-difference validation for the adjoint benchmark."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent


def replace_parameter(text: str, parameter: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*set\s+{re.escape(parameter)}\s*=\s*).*$", re.MULTILINE)
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise ValueError(f"Could not find unique parameter line for '{parameter}'.")
    return updated


def insert_plugin_material_model(text: str, plugin: Path) -> str:
    if "set Additional shared libraries" not in text:
        text = f"set Additional shared libraries            = {plugin}\n" + text
    else:
        text = replace_parameter(text, "Additional shared libraries", str(plugin))

    material_header = "subsection Material model\n  set Model name = simple"
    if material_header not in text:
        raise ValueError("Could not find the simple material model block.")
    text = text.replace(material_header,
                        "subsection Material model\n  set Model name = adjoint cellwise perturbation",
                        1)

    return text


def write_case_prm(
    base_text: str,
    prm_path: Path,
    output_directory: Path,
) -> None:
    text = replace_parameter(base_text, "Output directory", str(output_directory))
    prm_path.write_text(text)


def run_aspect(aspect: Path, prm_path: Path, log_path: Path, environment: dict[str, str] | None = None) -> None:
    aspect_executable = aspect.resolve() if aspect.exists() else aspect
    env = os.environ.copy()
    if environment is not None:
        env.update(environment)
    with log_path.open("w") as log:
        result = subprocess.run(
            [str(aspect_executable), str(prm_path)],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )

    if result.returncode != 0:
        tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-40:])
        raise RuntimeError(f"ASPECT failed for {prm_path}:\n{tail}")


def read_objective(output_directory: Path, objective_name: str) -> float:
    kernel_file = output_directory / "adjoint_kernels_rank_00000.txt"
    with kernel_file.open() as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 3 and fields[0] == "# objective_value" and fields[1] == objective_name:
                return float(fields[2])
    raise ValueError(f"Could not find objective '{objective_name}' in {kernel_file}.")


def read_property_kernels(output_directory: Path) -> list[dict[str, float | str]]:
    kernel_file = output_directory / "adjoint_kernels_rank_00000.txt"
    contributions: dict[tuple[str, str, int], dict[str, float | str]] = {}
    with kernel_file.open() as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            objective, _term, property_name, cell_index, cell_volume, value = line.rstrip("\n").split("\t")
            if property_name not in {"density", "viscosity"}:
                continue

            key = (objective, property_name, int(cell_index))
            if key not in contributions:
                contributions[key] = {
                    "objective": objective,
                    "property": property_name,
                    "active_cell_index": int(cell_index),
                    "cell_volume": float(cell_volume),
                    "value": 0.0,
                }
            contributions[key]["value"] = float(contributions[key]["value"]) + float(value)

    return list(contributions.values())


def selected_kernel_rows(
    rows: list[dict[str, float | str]],
    objective_name: str,
    cells_per_property: int,
    selection: str,
    explicit_cells: list[int] | None,
) -> list[dict[str, float | str]]:
    selected: list[dict[str, float | str]] = []
    for property_name in ["density", "viscosity"]:
        subset = [
            row for row in rows
            if row["objective"] == objective_name and row["property"] == property_name
        ]
        if explicit_cells is not None:
            selected.extend(row for row in subset if int(row["active_cell_index"]) in explicit_cells)
        else:
            if selection == "largest":
                selected.extend(
                    sorted(
                        subset,
                        key=lambda row: abs(float(row["value"]) * float(row["cell_volume"])),
                        reverse=True,
                    )[:cells_per_property]
                )
            else:
                sorted_subset = sorted(
                    subset,
                    key=lambda row: float(row["value"]) * float(row["cell_volume"]),
                )
                if cells_per_property >= len(sorted_subset):
                    selected.extend(sorted_subset)
                else:
                    selected_indices = {
                        round(i * (len(sorted_subset) - 1) / (cells_per_property - 1))
                        for i in range(cells_per_property)
                    }
                    selected.extend(sorted_subset[i] for i in sorted(selected_indices))
    return selected


def l2_relative_error(fd_values: list[float], adjoint_values: list[float]) -> float:
    numerator = sum((adjoint - fd) * (adjoint - fd)
                    for fd, adjoint in zip(fd_values, adjoint_values))
    denominator = sum(fd * fd for fd in fd_values)
    return math.sqrt(numerator / denominator) if denominator > 0.0 else 0.0


def pearson_correlation(fd_values: list[float], adjoint_values: list[float]) -> float:
    if len(fd_values) < 2:
        return 1.0

    mean_fd = sum(fd_values) / len(fd_values)
    mean_adjoint = sum(adjoint_values) / len(adjoint_values)
    covariance = sum((fd - mean_fd) * (adjoint - mean_adjoint)
                     for fd, adjoint in zip(fd_values, adjoint_values))
    variance_fd = sum((fd - mean_fd) * (fd - mean_fd) for fd in fd_values)
    variance_adjoint = sum((adjoint - mean_adjoint) * (adjoint - mean_adjoint)
                           for adjoint in adjoint_values)
    denominator = math.sqrt(variance_fd * variance_adjoint)
    return covariance / denominator if denominator > 0.0 else 1.0


def plot_results(rows: list[dict[str, float | str]], output: Path) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "axes.linewidth": 0.8,
        }
    )

    colors = {
        "density": "#2a9d8f",
        "viscosity": "#e76f51",
    }

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), constrained_layout=True)
    fig.patch.set_facecolor("white")

    for axis, control in zip(axes, ["density", "viscosity"]):
        subset = [row for row in rows if row["control"] == control]
        fd_values = [float(row["finite_difference"]) for row in subset]
        adjoint_values = [float(row["adjoint"]) for row in subset]
        relative_errors = [float(row["relative_error"]) for row in subset]

        axis.set_facecolor("#fbfbfa")
        axis.grid(True, which="major", color="#e7e3dd", linewidth=0.7)
        axis.axhline(0.0, color="#c8c4bd", linewidth=0.8, zorder=0)
        axis.axvline(0.0, color="#c8c4bd", linewidth=0.8, zorder=0)
        axis.scatter(fd_values,
                     adjoint_values,
                     s=42,
                     color=colors[control],
                     alpha=0.9,
                     edgecolor="white",
                     linewidth=0.7,
                     zorder=3)
        if fd_values and adjoint_values:
            raw_lower = min(fd_values + adjoint_values)
            raw_upper = max(fd_values + adjoint_values)
            if raw_lower >= 0.0:
                lower = 0.0
                upper = raw_upper * 1.08
            elif raw_upper <= 0.0:
                lower = raw_lower * 1.08
                upper = 0.0
            else:
                bound = max(abs(raw_lower), abs(raw_upper))
                lower = -1.08 * bound
                upper = 1.08 * bound
            axis.plot([lower, upper],
                      [lower, upper],
                      color="#202020",
                      linewidth=1.15,
                      label="1:1",
                      zorder=2)
            axis.set_xlim(lower, upper)
            axis.set_ylim(lower, upper)
            axis.set_aspect("equal", adjustable="box")

        annotation = (
            f"n = {len(subset)}\n"
            f"max rel. err. = {max(relative_errors):.2e}\n"
            f"L2 rel. err. = {l2_relative_error(fd_values, adjoint_values):.2e}\n"
            f"r = {pearson_correlation(fd_values, adjoint_values):.6f}"
        )
        axis.text(0.04,
                  0.96,
                  annotation,
                  transform=axis.transAxes,
                  va="top",
                  ha="left",
                  fontsize=8,
                  color="#333333",
                  bbox={
                      "boxstyle": "round,pad=0.32",
                      "facecolor": "white",
                      "edgecolor": "#ded8cf",
                      "linewidth": 0.7,
                      "alpha": 0.92,
                  })
        axis.set_title(control.capitalize(), pad=8)
        axis.set_xlabel("Finite difference")
        axis.set_ylabel("Adjoint kernel prediction")
        axis.ticklabel_format(axis="both", style="sci", scilimits=(-2, 3), useOffset=False)
        axis.legend(frameon=False, loc="lower right", handlelength=2.4)
        for spine in axis.spines.values():
            spine.set_color("#222222")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, facecolor=fig.get_facecolor())


def parse_cell_list(value: str | None) -> list[int] | None:
    if value is None:
        return None
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aspect", type=Path, default=Path("./aspect"), help="Path to the ASPECT executable.")
    parser.add_argument(
        "--plugin",
        type=Path,
        default=BENCHMARK_DIR / "libadjoint_dynamic_topography_3d_spherical.so",
        help="Path to the benchmark material plugin.",
    )
    parser.add_argument(
        "--prm",
        type=Path,
        default=BENCHMARK_DIR / "adjoint_dynamic_topography.prm",
        help="Base benchmark prm file.",
    )
    parser.add_argument("--workdir", type=Path, default=Path("fd_validation_work"))
    parser.add_argument("--density-step", type=float, default=1.0)
    parser.add_argument("--viscosity-step", type=float, default=1.0e18)
    parser.add_argument("--objective", default="dynamic topography")
    parser.add_argument("--cells-per-property", type=int, default=16)
    parser.add_argument("--selection", choices=["quantile", "largest"], default="quantile")
    parser.add_argument("--cells", help="Comma-separated active cell indices. Defaults to the largest adjoint kernels.")
    parser.add_argument(
        "--csv",
        type=Path,
        default=BENCHMARK_DIR / "doc" / "fd_vs_adjoint_validation.csv",
        help="Validation table to write.",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=BENCHMARK_DIR / "doc" / "fd_vs_adjoint_validation.png",
        help="Validation figure to write.",
    )
    args = parser.parse_args()

    plugin = args.plugin.resolve()
    variant_plugin = plugin.with_name(plugin.name.replace(".so", ".release.so"))
    if not plugin.exists() and not variant_plugin.exists():
        raise FileNotFoundError(f"Benchmark plugin not found: {plugin} or {variant_plugin}")

    base_text = insert_plugin_material_model(args.prm.read_text(), plugin)
    args.workdir.mkdir(parents=True, exist_ok=True)

    baseline_dir = args.workdir / "baseline"
    if baseline_dir.exists():
        shutil.rmtree(baseline_dir)
    baseline_dir.mkdir(parents=True)
    write_case_prm(base_text, baseline_dir / "baseline.prm", baseline_dir / "output")
    run_aspect(args.aspect, baseline_dir / "baseline.prm", baseline_dir / "aspect.log")

    kernels = read_property_kernels(baseline_dir / "output")
    selected_rows = selected_kernel_rows(
        kernels,
        args.objective,
        args.cells_per_property,
        args.selection,
        parse_cell_list(args.cells),
    )

    rows: list[dict[str, float | str]] = []
    for kernel in selected_rows:
        property_name = str(kernel["property"])
        cell_index = int(kernel["active_cell_index"])
        step = args.density_step if property_name == "density" else args.viscosity_step
        objectives: dict[str, float] = {}

        for sign_name, sign in [("minus", -1.0), ("plus", 1.0)]:
            case_name = f"{property_name}_{cell_index}_{sign_name}"
            case_dir = args.workdir / case_name
            if case_dir.exists():
                shutil.rmtree(case_dir)
            case_dir.mkdir(parents=True)
            write_case_prm(
                base_text,
                case_dir / f"{case_name}.prm",
                case_dir / "output",
            )
            run_aspect(
                args.aspect,
                case_dir / f"{case_name}.prm",
                case_dir / "aspect.log",
                {
                    "ASPECT_ADJOINT_FD_CELL": str(cell_index),
                    "ASPECT_ADJOINT_FD_PROPERTY": property_name,
                    "ASPECT_ADJOINT_FD_AMPLITUDE": f"{sign * step:.16g}",
                },
            )
            objectives[sign_name] = read_objective(case_dir / "output", args.objective)

        finite_difference = (objectives["plus"] - objectives["minus"]) / (2.0 * step)
        adjoint = float(kernel["value"]) * float(kernel["cell_volume"])
        scale = max(abs(finite_difference), abs(adjoint), 1.0e-300)
        rows.append(
            {
                "control": property_name,
                "active_cell_index": cell_index,
                "cell_volume": float(kernel["cell_volume"]),
                "step": step,
                "objective_minus": objectives["minus"],
                "objective_plus": objectives["plus"],
                "finite_difference": finite_difference,
                "adjoint": adjoint,
                "absolute_error": abs(finite_difference - adjoint),
                "relative_error": abs(finite_difference - adjoint) / scale,
            }
        )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    plot_results(rows, args.figure)

    for row in rows:
        print(
            f"{row['control']} cell {row['active_cell_index']}: "
            f"FD={row['finite_difference']:.8e}, adjoint={row['adjoint']:.8e}, "
            f"relative_error={row['relative_error']:.3e}"
        )


if __name__ == "__main__":
    main()
