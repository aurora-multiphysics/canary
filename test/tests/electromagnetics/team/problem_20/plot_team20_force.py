#!/usr/bin/env python3
"""Plot the computed force on the pole against the measurements.

Reads the postprocessor values written by team20.i at the benchmark's four
excitations and produces the figure embedded in
doc/content/problems/team_problem_20.md:

    team20_force.png   the force on the pole against ampere-turns

The test suite writes all four cases; to produce them by hand,

    <app> -i team20.i
    <app> -i team20.i ampere_turns=3000 \\
          Outputs/ReportedPostprocessors/file_base=OutputData/TEAM20CSV_3000AT

and likewise for 4500 and 5000 AT.

The problem statement publishes no measurements. They are reported in the summary
of results (Takahashi, Nakata and Morishige, COMPEL 14, 1995), which is not openly
available, and the values below are read off Fig. 6 of Hafla et al., ICS
Newsletter 13(2), 2006, which plots them. They are good to about half a newton.

Usage:
    ./plot_team20_force.py                  # after running the cases above
    ./plot_team20_force.py --outdir /tmp
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the other TEAM figures.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
MEASURED = "#e8710a"
COMPUTED = "#2a78d6"

# Ampere-turns, the CSV each is written to, and the measured force (N), read off
# the published figure.
CASES = [
    (1000, "TEAM20CSV.csv", 8.2),
    (3000, "TEAM20CSV_3000AT.csv", 54.5),
    (4500, "TEAM20CSV_4500AT.csv", 75.1),
    (5000, "TEAM20CSV_5000AT.csv", 80.2),
]


def read_force(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return float(rows[-1]["PoleForceZ"])


def _style(ax):
    """The axis styling used by the other figures in the documentation."""
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, solid_capstyle="butt")
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(which="both", colors=TEXT_SECONDARY, labelsize=9, length=0)


def plot(csv_dir, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    at, measured, computed = [], [], []
    for amps, name, meas in CASES:
        path = csv_dir / name
        if not path.exists():
            sys.exit(f"{path} not found; run the {amps} AT case first (see --help)")
        at.append(amps)
        measured.append(meas)
        # The pole is pulled down, towards the yoke under it; the benchmark
        # tabulates the size of the force.
        computed.append(-read_force(path))

    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    fig.suptitle("Force on the centre pole, TEAM Problem 20", color=TEXT_PRIMARY,
                 fontsize=12, x=0.012, ha="left", y=0.975)

    ax.plot(at, measured, color=MEASURED, linewidth=2.0, marker="o", markersize=6.0,
            markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
    ax.plot(at, computed, color=COMPUTED, linewidth=2.0, marker="o", markersize=6.0,
            markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
    # Direct labels at the right-hand end rather than a legend.
    for series, colour, label, dy in ((measured, MEASURED, "Measured", -9),
                                      (computed, COMPUTED, "Computed", 9)):
        ax.annotate(label, (at[-1], series[-1]), xytext=(8, dy), textcoords="offset points",
                    color=colour, fontsize=9, va="center")

    ax.set_xlim(0, 5600)
    ax.set_ylim(0, 100)
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    ax.set_xlabel("Excitation (AT)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Downward force on the pole (N)", color=TEXT_SECONDARY, fontsize=10)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)

    for amps, meas, comp in zip(at, measured, computed):
        print("%5d AT  measured %5.1f N  computed %6.2f N  (%+.1f%%)"
              % (amps, meas, comp, 100 * (comp / meas - 1)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv-dir", type=Path, default=HERE / "OutputData",
                   help="directory holding the postprocessor CSVs (default: %(default)s)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="directory the figure is written to (default: %(default)s)")
    args = p.parse_args()
    out = args.outdir / "team20_force.png"
    plot(args.csv_dir, out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
