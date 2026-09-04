#!/usr/bin/env python3
"""Repeat the TEAM Problem 4 convergence study and Richardson-extrapolate.

Runs team4_induced_field.i over a sequence of refinement levels, extracts the two
figures of merit at each level, estimates the observed order of convergence, and
extrapolates to the zero-discretisation-error limit. Prints the results both as a
readable summary and as a markdown table for doc/content/problems/team_problem_4.md.

Two studies are available:

  timestep  halve dt from --dt0 for --levels levels, on the supplied mesh. This is
            the one that matters: MFEMTransient offers only implicit Euler, so the
            time discretisation is first-order and dominates the error.

  mesh      uniformly refine the mesh, at a fixed --dt. Each level multiplies the
            element count by eight, so level 2 is already a large job; the default
            of two levels bounds the spatial error without trying to resolve an
            order from it.

Figures of merit, both quoted for the whole brick rather than the modelled octant:

  peak current   the largest total current circulating around the hole. The octant
                 is cut at z = 0, so the postprocessor value is doubled.
  peak power     the largest ohmic dissipation, summed over all eight octants.

Usage:
    ./run_convergence_study.py --exe /path/to/moose_test-opt
    ./run_convergence_study.py --study mesh --levels 2 --exe ...
    ./run_convergence_study.py --reuse            # re-analyse without re-running
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = HERE / "team4_induced_field.i"

# The octant is one eighth of the brick and is cut at z = 0 through the current path.
OCTANT_TO_FULL_SECTION = 2.0
OCTANT_TO_WHOLE_BRICK = 8.0


def find_executable(explicit):
    """Resolve the MOOSE application to run."""
    for candidate in (explicit, os.environ.get("MOOSE_APP")):
        if candidate:
            path = Path(candidate).expanduser()
            if path.is_file():
                return path
            found = shutil.which(str(candidate))
            if found:
                return Path(found)
            sys.exit(f"error: executable not found: {candidate}")

    moose_dir = os.environ.get("MOOSE_DIR")
    if moose_dir:
        path = Path(moose_dir) / "test" / "moose_test-opt"
        if path.is_file():
            return path

    sys.exit(
        "error: no MOOSE application given.\n"
        "Pass --exe /path/to/app-opt, or set MOOSE_APP, or set MOOSE_DIR so that\n"
        "$MOOSE_DIR/test/moose_test-opt exists. It must be built with MFEM support."
    )


def levels_for(args):
    """Return [(tag, [extra cli args], description)] for the chosen study."""
    if args.study == "timestep":
        return [
            (
                f"dt_{args.dt0 / 2 ** n:g}",
                [f"Executioner/dt={args.dt0 / 2 ** n:.12g}"],
                f"dt = {args.dt0 / 2 ** n:g} s",
            )
            for n in range(args.levels)
        ]

    # Both apps must be refined together: MultiAppMFEMCopyTransfer copies degrees of
    # freedom directly and needs the two meshes to agree exactly.
    return [
        (
            f"refine_{n}",
            [
                f"Executioner/dt={args.dt:.12g}",
                f"Mesh/uniform_refine={n}",
                f"MultiApps/external_source_field/cli_args=Mesh/uniform_refine={n}",
            ],
            f"{n} uniform refinement(s)",
        )
        for n in range(args.levels)
    ]


def run_level(exe, tag, extra, workdir, mpi):
    """Run one level and return (csv path, wall time in seconds)."""
    base = workdir / tag
    command = [str(exe), "-i", str(INPUT)]
    if mpi > 1:
        command = ["mpiexec", "-n", str(mpi)] + command
    command += extra + [
        "Outputs/ParaViewDataCollection/enable=false",
        f"Outputs/ReportedPostprocessors/file_base={base}",
    ]

    start = time.monotonic()
    result = subprocess.run(
        command, cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    elapsed = time.monotonic() - start

    if result.returncode != 0:
        log = workdir / f"{tag}.log"
        log.write_text(result.stdout)
        sys.exit(f"error: {tag} failed (exit {result.returncode}); output in {log}")
    return Path(f"{base}.csv"), elapsed


def figures_of_merit(csv_path):
    """Return (peak total current in A, peak whole-brick power in W, steps)."""
    with open(csv_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    peak_current = 0.0
    peak_power = 0.0
    for row in rows:
        inner = abs(float(row["InnerLimbCurrent"]))
        outer = abs(float(row["OuterLimbCurrent"]))
        peak_current = max(peak_current, 0.5 * (inner + outer) * OCTANT_TO_FULL_SECTION)
        peak_power = max(peak_power, float(row["OhmicHeating"]) * OCTANT_TO_WHOLE_BRICK)
    return peak_current, peak_power, len(rows) - 1


def observed_orders(values):
    """Order of convergence from each triple of successive halvings."""
    orders = []
    for i in range(len(values) - 2):
        d1 = values[i + 1] - values[i]
        d2 = values[i + 2] - values[i + 1]
        orders.append(math.log2(d1 / d2) if d2 and d1 / d2 > 0 else float("nan"))
    return orders


def richardson(values, order, ratio=2.0):
    """Extrapolate to zero discretisation error from the two finest levels."""
    if len(values) < 2:
        return float("nan")
    return values[-1] + (values[-1] - values[-2]) / (ratio**order - 1.0)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--study", choices=("timestep", "mesh"), default="timestep")
    parser.add_argument("--levels", type=int, help="number of refinement levels")
    parser.add_argument(
        "--dt0",
        type=float,
        default=4e-3,
        help="coarsest dt for --study timestep; halving it --levels times lands on the shipped dt",
    )
    parser.add_argument("--dt", type=float, default=1e-3, help="fixed dt for --study mesh")
    parser.add_argument(
        "--order",
        type=float,
        default=1.0,
        help="order assumed for the extrapolation (default 1, the order of implicit Euler)",
    )
    parser.add_argument("--exe", help="MOOSE application built with MFEM support")
    parser.add_argument("--mpi", type=int, default=1, help="MPI ranks per run")
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "ConvergenceStudy", help="where CSVs are written"
    )
    parser.add_argument(
        "--reuse", action="store_true", help="reuse CSVs already present instead of re-running"
    )
    args = parser.parse_args()

    if args.levels is None:
        args.levels = 6 if args.study == "timestep" else 2
    if args.levels < 2:
        sys.exit("error: --levels must be at least 2 for an extrapolation")

    exe = None if args.reuse else find_executable(args.exe)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for tag, extra, description in levels_for(args):
        csv_path = args.output_dir / f"{tag}.csv"
        if args.reuse:
            if not csv_path.is_file():
                sys.exit(f"error: --reuse given but {csv_path} does not exist")
            elapsed = float("nan")
        else:
            print(f"running {description} ...", flush=True)
            csv_path, elapsed = run_level(exe, tag, extra, args.output_dir, args.mpi)
        current, power, steps = figures_of_merit(csv_path)
        results.append((description, steps, current, power, elapsed))

    currents = [r[2] for r in results]
    powers = [r[3] for r in results]
    current_limit = richardson(currents, args.order)
    power_limit = richardson(powers, args.order)

    print(f"\n{'level':<28} {'steps':>6} {'peak I (kA)':>12} {'peak P (W)':>11} {'err(I)':>8} "
          f"{'runtime':>8}")
    for description, steps, current, power, elapsed in results:
        err = 100.0 * abs(current - current_limit) / current_limit
        runtime = "--" if math.isnan(elapsed) else f"{elapsed:.0f} s"
        print(f"{description:<28} {steps:>6} {current / 1e3:>12.3f} {power:>11.2f} "
              f"{err:>7.1f}% {runtime:>8}")
    print(f"{'Richardson extrapolation':<28} {'--':>6} {current_limit / 1e3:>12.3f} "
          f"{power_limit:>11.2f} {'--':>8} {'--':>8}")

    orders = observed_orders(currents)
    if orders:
        print("\nobserved order of convergence in peak current: "
              + ", ".join(f"{p:.2f}" for p in orders))
    print(f"extrapolation assumes order {args.order:g} (--order to change)")

    label = "$\\Delta t$ (s)" if args.study == "timestep" else "Refinements"
    print(f"\nmarkdown table:\n\n| {label} | Steps | Peak current (kA) | Peak power (W) "
          "| Error in peak current | Runtime |")
    print("| :- | -: | -: | -: | -: | -: |")
    for description, steps, current, power, elapsed in results:
        err = 100.0 * abs(current - current_limit) / current_limit
        runtime = "--" if math.isnan(elapsed) else f"{elapsed:.0f} s"
        value = description.replace("dt = ", "").replace(" s", "")
        print(f"| {value} | {steps} | {current / 1e3:.3f} | {power:.2f} | {err:.1f}% | {runtime} |")
    print(f"| Richardson extrapolation | -- | {current_limit / 1e3:.3f} | {power_limit:.2f} "
          "| -- | -- |")


if __name__ == "__main__":
    main()
