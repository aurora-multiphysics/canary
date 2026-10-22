#!/usr/bin/env python3
"""Plot the total circulating current in the FELIX brick against time.

Reads the postprocessor history written by team4_induced_field.i and produces the
figure embedded in doc/content/problems/team_problem_4.md.

The modelled octant is cut at z = 0, so each `semi_xsection` sideset spans only half
the brick thickness and the current a full cross-section of the brick carries is twice
what the postprocessors report. The two sidesets sit on opposite limbs of the loop
around the hole, so their currents are equal and opposite; this takes the mean of
their magnitudes and reports the spread as a consistency check.

Usage:
    ./plot_team4_current.py                       # after running team4_induced_field.i
    ./plot_team4_current.py --csv other.csv --output /tmp/team4.png
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_CSV = HERE / "OutputData" / "TEAM4CSV.csv"
DEFAULT_PNG = REPO_ROOT / "doc" / "content" / "media" / "team4_current.png"

# Chart surface and ink. Kept as named roles so the whole figure restyles from here.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
SERIES = "#2a78d6"

# The octant cuts the brick in half at z = 0.
OCTANT_TO_FULL_SECTION = 2.0


def read_history(path):
    """Return (times in ms, total current in A, worst relative limb mismatch)."""
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError:
        sys.exit(
            f"error: {path} not found.\n"
            "Run team4_induced_field.i first; it writes OutputData/TEAM4CSV.csv."
        )

    required = {"time", "InnerLimbCurrent", "OuterLimbCurrent"}
    missing = required - set(rows[0] if rows else {})
    if missing:
        sys.exit(f"error: {path} is missing column(s): {', '.join(sorted(missing))}")

    times, totals, mismatch = [], [], 0.0
    for row in rows:
        inner = abs(float(row["InnerLimbCurrent"]))
        outer = abs(float(row["OuterLimbCurrent"]))
        mean = 0.5 * (inner + outer)
        if mean > 0.0:
            mismatch = max(mismatch, abs(inner - outer) / mean)
        times.append(float(row["time"]) * 1e3)
        totals.append(mean * OCTANT_TO_FULL_SECTION)
    return times, totals, mismatch


def plot(times, totals, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator

    fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # Recessive grid, drawn behind the data, hairline and solid.
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, solid_capstyle="butt")
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)

    ax.plot(
        times,
        totals,
        color=SERIES,
        linewidth=2.0,
        solid_capstyle="round",
        solid_joinstyle="round",
    )

    # Label the extreme only; the axes carry every other value.
    peak_i = max(range(len(totals)), key=totals.__getitem__)
    peak_t, peak_a = times[peak_i], totals[peak_i]
    ax.plot(
        peak_t,
        peak_a,
        marker="o",
        markersize=8,
        color=SERIES,
        markeredgecolor=SURFACE,
        markeredgewidth=2.0,
        zorder=3,
    )
    ax.annotate(
        f"{peak_a / 1e3:.2f} kA at {peak_t:.1f} ms",
        xy=(peak_t, peak_a),
        xytext=(10, 8),
        textcoords="offset points",
        color=TEXT_PRIMARY,
        fontsize=10,
    )

    ax.set_xlabel("Time (ms)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Current (kA)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_title(
        "Total current circulating around the hole in the FELIX brick",
        color=TEXT_PRIMARY,
        fontsize=12,
        loc="left",
        pad=12,
    )

    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1e3:,.1f}"))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=0)
    ax.set_xlim(0, max(times))
    ax.set_ylim(0, max(totals) * 1.15)

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    return peak_t, peak_a


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--csv", type=Path, default=DEFAULT_CSV, help=f"input CSV (default: {DEFAULT_CSV})"
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_PNG, help=f"output PNG (default: {DEFAULT_PNG})"
    )
    args = parser.parse_args()

    times, totals, mismatch = read_history(args.csv)
    peak_t, peak_a = plot(times, totals, args.output)

    print(f"wrote {args.output}")
    print(f"  peak total current   {peak_a:,.1f} A at {peak_t:.3f} ms")
    print(f"  limb current mismatch {mismatch:.2e} (relative, worst over the run)")


if __name__ == "__main__":
    main()
