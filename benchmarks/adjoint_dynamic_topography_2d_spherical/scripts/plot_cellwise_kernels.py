#!/usr/bin/env python3
"""Plot cellwise dynamic-topography adjoint density and viscosity kernels."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib.pyplot as plt
import pandas as pd


def read_adjoint_table(path: Path, columns: list[str]) -> pd.DataFrame:
    rows: list[list[str]] = []
    with path.open() as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))

    data = pd.DataFrame(rows, columns=columns)
    for column in columns:
        if column not in {"objective", "term", "property", "control"}:
            data[column] = pd.to_numeric(data[column])
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output_directory",
        type=Path,
        help="ASPECT output directory containing adjoint_kernels_rank_00000.txt",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("doc/cellwise_density_viscosity_kernels.png"),
        help="Figure path to write.",
    )
    args = parser.parse_args()

    kernel_file = args.output_directory / "adjoint_kernels_rank_00000.txt"
    kernels = read_adjoint_table(
        kernel_file,
        ["objective", "term", "property", "active_cell_index", "cell_volume", "value"],
    )

    terms = kernels[kernels["property"].isin(["density", "viscosity"])].copy()
    grouped = (
        terms.groupby(["property", "active_cell_index"], as_index=False)["value"]
        .sum()
        .sort_values(["property", "active_cell_index"])
    )

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)
    for axis, property_name in zip(axes, ["density", "viscosity"]):
        subset = grouped[grouped["property"] == property_name]
        axis.axhline(0.0, color="0.3", linewidth=0.8)
        axis.plot(subset["active_cell_index"], subset["value"], marker="o", linewidth=1.2)
        axis.set_title(property_name)
        axis.set_xlabel("active cell index")
        axis.set_ylabel("cell-average kernel")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
