#!/usr/bin/env python3
"""Plot the computed flux density up the centre plate against the measurements.

Reads the postprocessor values written by team13.i and produces the figure embedded
in doc/content/problems/team_problem_13.md.

The benchmark was measured on two specimens whose plate-to-channel gap came out at
0.52 mm and 0.47 mm against a nominal 0.5 mm, and that gap carries most of the
reluctance of the magnetic circuit. The two sets of readings are therefore drawn as
a band rather than as a single curve: the width of the band is the experimental
spread the model is being asked to land inside.

Usage:
    ./plot_team13_flux.py                       # after running team13.i
    ./plot_team13_flux.py --csv other.csv --output /tmp/team13.png
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_CSV = HERE / "OutputData" / "TEAM13CSV.csv"
DEFAULT_PNG = REPO_ROOT / "doc" / "content" / "media" / "team13_flux.png"

# Chart surface and ink, shared with the TEAM Problem 4 figure.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
COMPUTED = "#2a78d6"
MEASURED = "#e8710a"

# Measured average flux density (T) at 1000 AT, from the problem statement, for the
# two specimens, whose gaps G were 0.52 mm and 0.47 mm.
PLATE_Z = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
PLATE_G52 = [1.333, 1.329, 1.286, 1.225, 1.129, 0.985, 0.655]
PLATE_G47 = [1.354, 1.339, 1.304, 1.245, 1.138, 0.982, 0.674]
PLATE_KEYS = ["PlateB%03d" % int(z) for z in PLATE_Z]

# The two sections in the channel: position 14, halfway along the top flange, and
# position 25, at the foot of the web.
CHANNEL = [
    ("FlangeB060", "flange at x = 60 mm (position 14)", 0.809, 0.819),
    ("WebB000", "web at z = 0 (position 25)", 0.985, 0.995),
]


def read_values(path):
    """Return the last row of the postprocessor CSV as a dict of floats."""
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError:
        sys.exit(
            f"error: {path} not found.\nRun team13.i first; it writes OutputData/TEAM13CSV.csv."
        )
    if not rows:
        sys.exit(f"error: {path} holds no rows")

    missing = [k for k in PLATE_KEYS + [c[0] for c in CHANNEL] if k not in rows[-1]]
    if missing:
        sys.exit(f"error: {path} is missing column(s): {', '.join(missing)}")
    return {k: float(v) for k, v in rows[-1].items() if v != ""}


def plot(computed, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, solid_capstyle="butt")
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)

    lower = [min(a, b) for a, b in zip(PLATE_G52, PLATE_G47)]
    upper = [max(a, b) for a, b in zip(PLATE_G52, PLATE_G47)]
    ax.fill_between(PLATE_Z, lower, upper, color=MEASURED, alpha=0.18, linewidth=0)
    ax.plot(
        PLATE_Z,
        [0.5 * (a + b) for a, b in zip(PLATE_G52, PLATE_G47)],
        color=MEASURED,
        linewidth=2.0,
        marker="o",
        markersize=8,
        markeredgecolor=SURFACE,
        markeredgewidth=2.0,
        label="Measured (band spans the two specimens)",
    )
    ax.plot(
        PLATE_Z,
        [computed[k] for k in PLATE_KEYS],
        color=COMPUTED,
        linewidth=2.0,
        marker="s",
        markersize=8,
        markeredgecolor=SURFACE,
        markeredgewidth=2.0,
        label="Computed",
    )

    ax.set_xlabel("Height up the plate, z (mm)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Average flux density (T)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_title(
        "Flux density in the centre plate at 1000 AT",
        color=TEXT_PRIMARY,
        fontsize=12,
        loc="left",
        pad=12,
    )
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=0)
    ax.set_xlim(-2, 62)
    ax.set_ylim(0, max(upper) * 1.2)

    legend = ax.legend(loc="lower left", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)


def report(computed):
    print("  centre plate (positions 1 to 7)")
    print("    %-8s %10s %10s %10s" % ("z (mm)", "computed", "G=0.52", "G=0.47"))
    worst = 0.0
    for z, key, g52, g47 in zip(PLATE_Z, PLATE_KEYS, PLATE_G52, PLATE_G47):
        c = computed[key]
        band = (min(g52, g47), max(g52, g47))
        gap = 0.0 if band[0] <= c <= band[1] else min(abs(c - band[0]), abs(c - band[1])) / c
        worst = max(worst, gap)
        print("    %-8.0f %10.4f %10.3f %10.3f" % (z, c, g52, g47))
    print("  channel")
    for key, label, g52, g47 in CHANNEL:
        print("    %-34s %10.4f %10.3f %10.3f" % (label, computed[key], g52, g47))
    print("  worst distance from the measured band, over the plate: %.1f%%" % (100 * worst))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--csv", type=Path, default=DEFAULT_CSV, help=f"input CSV (default: {DEFAULT_CSV})"
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_PNG, help=f"output PNG (default: {DEFAULT_PNG})"
    )
    parser.add_argument("--no-plot", action="store_true", help="only print the comparison")
    args = parser.parse_args()

    computed = read_values(args.csv)
    if not args.no_plot:
        plot(computed, args.output)
        print(f"wrote {args.output}")
    report(computed)


if __name__ == "__main__":
    main()
