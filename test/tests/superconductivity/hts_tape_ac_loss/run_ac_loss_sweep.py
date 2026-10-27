#!/usr/bin/env python3
"""Sweep the HTS tape AC loss benchmark over the transport current and compare.

Runs hts_tape.i at each ratio I_m/I_c the benchmark asks for, 0.1 to 0.9, and
compares the loss per cycle with the Norris formula for a thin strip in the
critical state,

    Q = (mu0 I_c^2 / pi) [(1 - F) ln(1 - F) + (1 + F) ln(1 + F) - F^2],  F = I_m/I_c,

and the current density across the tape at the peak of the current with the
critical state profile of the same strip (Norris; Brandt and Indenbom),

    J/J_c = (2/pi) arctan(sqrt((a^2 - b^2) / (b^2 - x^2)))   for |x| < b,
    J/J_c = 1                                                for b < |x| < a,

with a the half-width and b = a sqrt(1 - F^2) the edge of the region the current
has not yet reached.

Backward Euler makes the loss first order in the time step, so each ratio is run
at --steps-per-cycle and at twice that, and Richardson-extrapolated to the limit.
The results are printed as the markdown table in
doc/content/problems/hts_tape_ac_loss.md, and the figures there are written with
--plot.

Usage:
    ./run_ac_loss_sweep.py --exe /path/to/moose_test-opt --plot
    ./run_ac_loss_sweep.py --reuse --plot     # re-analyse without re-running
"""

import argparse
import csv
import glob
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = HERE / "hts_tape.i"
OUT = HERE / "OutputData" / "sweep"
MEDIA = HERE.parents[3] / "doc" / "content" / "media"

MU0 = 4e-7 * math.pi
I_C = 300.0  # A
HALF_WIDTH = 6e-3  # m
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
PROFILE_RATIOS = [0.5, 0.9]


def norris_loss(f):
    """Loss per cycle per unit length of a thin strip in the critical state, J/m."""
    return MU0 * I_C**2 / math.pi * ((1 - f) * math.log(1 - f) + (1 + f) * math.log(1 + f) - f**2)


def bean_profile(x, f):
    """J/J_c across a thin strip carrying a transport current F I_c, at its peak."""
    b = HALF_WIDTH * math.sqrt(1 - f**2)
    if abs(x) >= b:
        return 1.0
    return 2 / math.pi * math.atan(math.sqrt((HALF_WIDTH**2 - b**2) / (b**2 - x**2)))


def find_executable(explicit):
    for candidate in (explicit, os.environ.get("MOOSE_APP")):
        if candidate:
            found = shutil.which(str(candidate)) or (candidate if Path(candidate).is_file() else None)
            if found:
                return found
            sys.exit(f"error: executable not found: {candidate}")
    moose_dir = os.environ.get("MOOSE_DIR")
    if moose_dir and (Path(moose_dir) / "test" / "moose_test-opt").is_file():
        return str(Path(moose_dir) / "test" / "moose_test-opt")
    sys.exit("error: no MOOSE application given; pass --exe, or set MOOSE_APP or MOOSE_DIR")


def tag(n, f, steps):
    return f"n{n:g}_F{f:g}_s{steps}"


def run(exe, n, f, steps, extra):
    """Run one case, one at a time; the solve is serial and small."""
    base = OUT / tag(n, f, steps)
    cmd = [exe, "-i", INPUT.name, f"n_value={n}", f"current_ratio={f}", f"steps_per_cycle={steps}",
           "Outputs/ParaViewDataCollection/enable=false", "Solvers/newton/print_level=0",
           f"Outputs/ReportedPostprocessors/file_base={base}_loss",
           f"Outputs/PeakProfile/file_base={base}_profile"] + extra
    print(f"running n = {n:g}, I_m/I_c = {f:g} at {steps} steps per cycle", flush=True)
    subprocess.run(cmd, cwd=HERE, check=True, stdout=subprocess.DEVNULL)


def loss(n, f, steps):
    with open(OUT / f"{tag(n, f, steps)}_loss.csv") as fh:
        return float(list(csv.DictReader(fh))[-1]["LossPerCycle"])


def profile(n, f, steps):
    (path,) = glob.glob(str(OUT / f"{tag(n, f, steps)}_profile_current_density_profile_*.csv"))
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return [float(r["x_0"]) for r in rows], [float(r["j_field_0"]) for r in rows]


def plot(results, n_values, steps):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    j_c = I_C / (2 * HALF_WIDTH * 1e-6)
    MEDIA.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    fs = [0.05 + 0.005 * i for i in range(181)]
    ax.loglog(fs, [norris_loss(f) for f in fs], "k-", label="Norris, thin strip")
    for n, marker in zip(n_values, "os^"):
        rs = [r for r in results if r[0] == n]
        ax.loglog([r[1] for r in rs], [r[4] for r in rs], marker, mfc="none",
                  label=f"MOOSE, n = {n:g}")
    ax.set_xlabel("$I_m/I_c$")
    ax.set_ylabel("loss per cycle (J/m)")
    ax.set_title("AC transport loss of a 12 mm HTS tape at 50 Hz")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(MEDIA / "hts_tape_ac_loss.png", dpi=150)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for f, colour in zip(PROFILE_RATIOS, ("C0", "C3")):
        for n, style in zip(n_values, ("-", ":", "-.")):
            x, j = profile(n, f, 2 * steps)
            ax.plot([1e3 * xi for xi in x], [ji / j_c for ji in j], style, color=colour,
                    label=f"MOOSE, n = {n:g}, $I_m/I_c$ = {f:g}")
        xs = [HALF_WIDTH * i / 600 for i in range(601)]
        ax.plot([1e3 * xi for xi in xs], [bean_profile(xi, f) for xi in xs], "--",
                color="k", lw=0.8, label="critical state" if f == PROFILE_RATIOS[0] else None)
    ax.set_xlabel("distance from the centre of the tape (mm)")
    ax.set_ylabel("$J/J_c$")
    ax.set_title("Current density across the tape at the peak of the current")
    ax.set_xlim(0, 1e3 * HALF_WIDTH)
    ax.set_ylim(0, 1.2)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    fig.savefig(MEDIA / "hts_tape_current_density.png", dpi=150)
    print(f"wrote figures to {MEDIA}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", help="MOOSE application built with MFEM")
    p.add_argument("--steps-per-cycle", type=int, default=400,
                   help="coarser of the two time steps run at each ratio")
    p.add_argument("--ratios", type=float, nargs="+", default=RATIOS)
    p.add_argument("--n-values", type=float, nargs="+", default=[25, 101],
                   help="power law exponents; the critical state is the limit of large n")
    p.add_argument("--reuse", action="store_true", help="analyse existing output only")
    p.add_argument("--plot", action="store_true", help="write the documentation figures")
    p.add_argument("extra", nargs="*", help="further command line arguments for the input")
    args = p.parse_args()

    s = args.steps_per_cycle
    if not args.reuse:
        exe = find_executable(args.exe)
        OUT.mkdir(parents=True, exist_ok=True)
        for n in args.n_values:
            for f in args.ratios:
                for steps in (s, 2 * s):
                    run(exe, n, f, steps, args.extra)

    results = []
    for n in args.n_values:
        for f in args.ratios:
            q1, q2 = loss(n, f, s), loss(n, f, 2 * s)
            results.append((n, f, q1, q2, 2 * q2 - q1, norris_loss(f)))

    print(f"\n| $n$ | $I_m/I_c$ | {s} steps (J/m) | {2 * s} steps (J/m) | extrapolated (J/m) "
          f"| Norris (J/m) | difference |")
    print("| - | - | - | - | - | - | - |")
    for n, f, q1, q2, q, qn in results:
        print(f"| {n:g} | {f:g} | {q1:.4e} | {q2:.4e} | {q:.4e} | {qn:.4e} "
              f"| {100 * (q / qn - 1):+.1f}% |")

    if args.plot:
        plot(results, args.n_values, s)


if __name__ == "__main__":
    main()
