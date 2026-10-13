#!/usr/bin/env python3
"""Plot the MFEM Cahn-Hilliard solution and compare it with the phase field module.

Produces the figures embedded in doc/content/problems/cahn_hilliard.md:

    cahn_hilliard_snapshots.png    c at three times, from cahn_hilliard.i
    cahn_hilliard_comparison.png   the free energy, and c along the centre line,
                                   from cahn_hilliard.i and the phase field module

and prints the difference between the two solutions at the end, node by node.

The snapshots need a run long enough to show coarsening, with ASCII ParaView
output every ten steps:

    moose_test-opt -i cahn_hilliard.i end_time=200 \\
        Outputs/ParaViewDataCollection/vtk_format=ASCII \\
        Outputs/ParaViewDataCollection/time_step_interval=10 \\
        Outputs/ParaViewDataCollection/file_base=OutputData/CahnHilliardLong \\
        Outputs/CSV/file_base=OutputData/CahnHilliardLong

The comparison needs the default run of cahn_hilliard.i, also with ASCII output,
and cahn_hilliard_phase_field.i, which needs an application with the phase field
module:

    moose_test-opt -i cahn_hilliard.i Outputs/ParaViewDataCollection/vtk_format=ASCII
    phase_field-opt -i cahn_hilliard_phase_field.i

Other resolutions are compared by giving both runs a file_base suffix, e.g. for
300 elements per side, with elements_per_side=300 and

    Outputs/ParaViewDataCollection/file_base=OutputData/CahnHilliard300
    Outputs/CSV/file_base=OutputData/CahnHilliard300

for cahn_hilliard.i and Outputs/CSV/file_base=OutputData/CahnHilliardPhaseField300
for cahn_hilliard_phase_field.i, then ./plot_cahn_hilliard.py --case 300.

Usage:
    ./plot_cahn_hilliard.py                  # after running the cases above
    ./plot_cahn_hilliard.py --outdir /tmp
"""

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "OutputData")
MEDIA = os.path.join(HERE, "..", "..", "..", "..", "doc", "content", "media")

DOMAIN_SIZE = 60.0
SNAPSHOT_TIMES = (20, 60, 200)


def read_csv(path):
    """Return the columns of a MOOSE CSV file as a dict of numpy arrays."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def read_vtu(path, field):
    """Return the points and the named point field of an ASCII MFEM ParaView file."""
    with open(path) as f:
        text = f.read()
    arrays = re.findall(r"<DataArray([^>]*)>(.*?)</DataArray>", text, re.S)
    if 'format="ascii"' not in arrays[0][0]:
        sys.exit(f"{path} is not ASCII; rerun with vtk_format=ASCII")
    points = np.array(arrays[0][1].split(), float).reshape(-1, 3)
    for header, body in arrays:
        if f'Name="{field}"' in header:
            return points[:, :2], np.array(body.split(), float)
    sys.exit(f"no field {field} in {path}")


def nodal_image(points, values, spacing):
    """Gather values on the nodes of a uniform grid, which MFEM repeats per element."""
    index = np.rint(points / spacing).astype(int)
    n = index.max() + 1
    image = np.full((n, n), np.nan)
    image[index[:, 1], index[:, 0]] = values
    return image


def vtu_at(base, time):
    """Return the serial ParaView file written at the given time, found from Run0.pvd."""
    pvd = os.path.join(base, "Run0", "Run0.pvd")
    if not os.path.exists(pvd):
        sys.exit(f"missing {pvd}; see the docstring for how to produce it")
    with open(pvd) as f:
        for step, cycle in re.findall(r'timestep="([^"]+)"[^>]*file="([^/"]+)/', f.read()):
            if np.isclose(float(step), time):
                return os.path.join(base, "Run0", cycle, "proc000000.vtu")
    sys.exit(f"no output at t = {time} in {pvd}")


def element_spacing(points):
    return np.min(np.diff(np.unique(np.round(points[:, 0], 9))))


def plot_snapshots(outdir):
    base = os.path.join(DATA, "CahnHilliardLong")
    fig, axes = plt.subplots(1, len(SNAPSHOT_TIMES), figsize=(12, 4.2), constrained_layout=True)
    for ax, time in zip(axes, SNAPSHOT_TIMES):
        points, c = read_vtu(vtu_at(base, time), "c")
        image = nodal_image(points, c, element_spacing(points))
        mesh = ax.imshow(
            image,
            origin="lower",
            extent=(0, DOMAIN_SIZE, 0, DOMAIN_SIZE),
            cmap="RdBu_r",
            vmin=-1.05,
            vmax=1.05,
        )
        ax.set_title(f"t = {time}")
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(mesh, ax=axes, shrink=0.85, label="c")
    path = os.path.join(outdir, "cahn_hilliard_snapshots.png")
    fig.savefig(path, dpi=110)
    print(f"wrote {path}")


def compare(outdir, case=""):
    mfem_base = os.path.join(DATA, "CahnHilliard" + case)
    pf_base = os.path.join(DATA, "CahnHilliardPhaseField" + case)
    mfem = read_csv(mfem_base + ".csv")
    pf = read_csv(pf_base + ".csv")
    pf_rows = pf["time"] > 0  # the phase field module also reports t = 0

    mfem_line = read_csv(sorted(glob.glob(mfem_base + "_c_line_*.csv"))[-1])
    pf_line = read_csv(
        sorted(glob.glob(pf_base + "_c_line_*.csv"))[-1]
    )

    fig, (ax_energy, ax_line) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    ax_energy.plot(pf["time"][pf_rows], pf["total_free_energy"][pf_rows], "o-", label="Phase field module")
    ax_energy.plot(mfem["time"], mfem["total_free_energy"], "x--", label="MFEM")
    ax_energy.set_xlabel("t")
    ax_energy.set_ylabel("Total free energy F")
    ax_energy.legend()
    ax_energy.grid(alpha=0.3)

    ax_line.plot(pf_line["x"], pf_line["c"], "-", label="Phase field module")
    ax_line.plot(mfem_line["x_0"], mfem_line["c_0"], "--", label="MFEM")
    ax_line.set_xlabel("x (y = 30)")
    ax_line.set_ylabel(f"c at t = {mfem['time'][-1]:g}")
    ax_line.legend()
    ax_line.grid(alpha=0.3)

    if not case:
        path = os.path.join(outdir, "cahn_hilliard_comparison.png")
        fig.savefig(path, dpi=110)
        print(f"wrote {path}")

    # Node by node at the end. MFEM writes each node once per element sharing it.
    points, c_mfem = read_vtu(vtu_at(mfem_base, mfem["time"][-1]), "c")
    spacing = element_spacing(points)
    image_mfem = nodal_image(points, c_mfem, spacing)
    nodal = read_csv(
        sorted(glob.glob(pf_base + "_c_nodal_*.csv"))[-1]
    )
    image_pf = nodal_image(np.column_stack((nodal["x"], nodal["y"])), nodal["c"], spacing)
    difference = image_mfem - image_pf
    print(f"t = {mfem['time'][-1]:g}, {image_pf.shape[0] - 1} elements per side")
    print(f"  max |c_MFEM - c_PF|      = {np.nanmax(np.abs(difference)):.3e}")
    print(f"  RMS (c_MFEM - c_PF)      = {np.sqrt(np.nanmean(difference**2)):.3e}")
    for time, f_mfem in zip(mfem["time"], mfem["total_free_energy"]):
        f_pf = pf["total_free_energy"][np.isclose(pf["time"], time)][0]
        print(f"  t = {time:4g}  F_MFEM = {f_mfem:.4f}  F_PF = {f_pf:.4f}  "
              f"relative difference = {(f_mfem - f_pf) / f_pf:+.2e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default=MEDIA, help="where to write the figures")
    parser.add_argument("--no-snapshots", action="store_true", help="skip the long run")
    parser.add_argument("--no-comparison", action="store_true", help="skip the comparison")
    parser.add_argument(
        "--case",
        default="",
        help="suffix of the file_base of both runs, to compare other resolutions; no figure is drawn",
    )
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if not args.no_snapshots:
        plot_snapshots(args.outdir)
    if not args.no_comparison:
        compare(args.outdir, args.case)


if __name__ == "__main__":
    main()
