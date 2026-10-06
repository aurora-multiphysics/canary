#!/usr/bin/env python3
"""Compare the TEAM Problem 12 results with the FELIX measurements.

Reads the postprocessor histories written by team12.i and the measurements in
team12_experiment.csv, and produces the figures embedded in
doc/content/problems/team_problem_12.md:

- team12_coupled.png: deflection at the free end and 167 mm from it, and the circulating
  current, for each steady field B_x;
- team12_uncoupled.png: the circulating current in copper, aluminium and brass beams
  with no steady field, when the beam does not move.

The measurements were digitised from Figs. 5, 8 and 9 of the LSTC test case document
for TEAM Problem 12 (EM-VAL-3.1), which reproduces them from Turner, Hua and Lee (1986).

The problem statement does not fix the sense of the steady field relative to the
measured deflection, and reversing B_x reverses the deflection and leaves everything
else unchanged. Both the deflections and the circulating current are therefore reported
with the sign that makes their first swing positive, as the measurements are.

Usage:
    ./plot_team12.py          # after running the tests that write OutputData/TEAM12CSV_*.csv
    ./plot_team12.py --output-dir /tmp
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_OUTPUT_DIR = REPO_ROOT / "doc" / "content" / "media"
EXPERIMENT_CSV = HERE / "team12_experiment.csv"

STEADY_FIELDS = ("0.2", "0.5", "0.7", "0.9")
MATERIALS = ("copper", "aluminium", "brass")

# Chart surface and ink. Kept as named roles so the whole figure restyles from here.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
# Categorical slots, in fixed order, validated together for colour vision deficiency.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")

# The experiment's own marker shapes: crosses at the free end, boxed pluses 167 mm from
# it, and triangle, cross and boxed plus for copper, aluminium and brass.
STATION_MARKERS = {"tip": "x", "station": "P"}
MATERIAL_MARKERS = {"copper": "^", "aluminium": "x", "brass": "P"}


def read_simulation(path):
    """Return {column: [values]} from a MOOSE postprocessor CSV, times in ms, or None if
    that case has not been run, in which case only the measurements are drawn."""
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError:
        print(f"warning: {path} not found; plotting the measurements only", file=sys.stderr)
        return None
    columns = defaultdict(list)
    for row in rows:
        for key, value in row.items():
            columns[key].append(float(value))
    columns["time"] = [t * 1e3 for t in columns["time"]]
    return columns


def read_experiment():
    """Return {(case, quantity): ([times in ms], [values])}."""
    data = defaultdict(lambda: ([], []))
    with open(EXPERIMENT_CSV, newline="") as handle:
        for row in csv.DictReader(handle):
            times, values = data[(row["case"], row["quantity"])]
            times.append(float(row["time"]) * 1e3)
            values.append(float(row["value"]))
    return data


def first_swing_sign(values):
    """The sign of the first excursion beyond a tenth of the largest one."""
    threshold = 0.1 * max(abs(v) for v in values)
    first = next(v for v in values if abs(v) > threshold)
    return 1.0 if first > 0 else -1.0


def oriented(values, sign=None):
    """values, signed so that their first swing is positive, or by sign if given."""
    sign = first_swing_sign(values) if sign is None else sign
    return [sign * v for v in values]


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0)
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=0)
    ax.axhline(0.0, color=GRID, linewidth=1.0, zorder=0)


def plot_measured(ax, times, values, colour, marker):
    ax.plot(
        times,
        values,
        linestyle="none",
        marker=marker,
        markersize=7,
        markerfacecolor="none" if marker != "x" else colour,
        markeredgecolor=colour,
        markeredgewidth=1.5,
        zorder=3,
    )


def plot_coupled(experiment, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(
        2, len(STEADY_FIELDS), figsize=(13.0, 6.4), dpi=200, sharex=True, sharey="row"
    )
    fig.patch.set_facecolor(SURFACE)

    for col, bx in enumerate(STEADY_FIELDS):
        sim = read_simulation(HERE / "OutputData" / f"TEAM12CSV_bx_{bx}.csv")
        case = f"bx_{bx}"

        ax = axes[0, col]
        style(ax)
        # Both stations take the sign of the free end, so they stay in phase.
        sign = first_swing_sign(sim["TipDeflection"]) if sim else 1.0
        for (name, column, quantity), colour in zip(
            (
                ("tip", "TipDeflection", "tip_deflection_mm"),
                ("station", "StationDeflection", "station_deflection_mm"),
            ),
            SERIES,
        ):
            if sim:
                ax.plot(
                    sim["time"],
                    oriented([d * 1e3 for d in sim[column]], sign),
                    color=colour,
                    linewidth=2.0,
                )
            plot_measured(ax, *experiment[(case, quantity)], colour, STATION_MARKERS[name])
        ax.set_title(f"$B_x$ = {bx} T", color=TEXT_PRIMARY, fontsize=11, loc="left")

        ax = axes[1, col]
        style(ax)
        if sim:
            ax.plot(sim["time"], oriented(sim["CirculatingCurrent"]), color=SERIES[0], linewidth=2.0)
        plot_measured(ax, *experiment[(case, "current_A")], SERIES[0], "x")
        ax.set_xlabel("Time (ms)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_xlim(0, 150)

    axes[0, 0].set_ylabel("Deflection (mm)", color=TEXT_SECONDARY, fontsize=10)
    axes[1, 0].set_ylabel("Circulating current (A)", color=TEXT_SECONDARY, fontsize=10)

    legend = [
        Line2D([], [], color=SERIES[0], linewidth=2.0, label="Free end, MOOSE"),
        Line2D([], [], color=SERIES[0], marker="x", linestyle="none", label="Free end, measured"),
        Line2D([], [], color=SERIES[1], linewidth=2.0, label="167 mm from free end, MOOSE"),
        Line2D(
            [], [], color=SERIES[1], marker="P", markerfacecolor="none", linestyle="none",
            label="167 mm from free end, measured",
        ),
    ]
    fig.legend(
        handles=legend, loc="upper center", ncol=4, frameon=False, fontsize=9,
        labelcolor=TEXT_PRIMARY, bbox_to_anchor=(0.5, 1.0),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)


def plot_uncoupled(experiment, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    for material, colour in zip(MATERIALS, SERIES):
        sim = read_simulation(HERE / "OutputData" / f"TEAM12CSV_{material}.csv")
        plot_measured(ax, *experiment[(material, "current_A")], colour, MATERIAL_MARKERS[material])
        if not sim:
            continue
        current = oriented(sim["CirculatingCurrent"])
        ax.plot(sim["time"], current, color=colour, linewidth=2.0)
        # Direct label at the peak of each curve.
        peak = max(range(len(current)), key=current.__getitem__)
        ax.annotate(
            material.capitalize(),
            xy=(sim["time"][peak], current[peak]),
            xytext=(10, 4),
            textcoords="offset points",
            color=TEXT_PRIMARY,
            fontsize=10,
        )

    ax.set_xlabel("Time (ms)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Circulating current (A)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_title(
        "Circulating current with no steady field: lines MOOSE, markers measured",
        color=TEXT_PRIMARY,
        fontsize=11,
        loc="left",
        pad=12,
    )
    ax.set_xlim(0, 40)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"directory for the PNGs (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()

    experiment = read_experiment()
    plot_coupled(experiment, args.output_dir / "team12_coupled.png")
    plot_uncoupled(experiment, args.output_dir / "team12_uncoupled.png")
    print(f"wrote {args.output_dir / 'team12_coupled.png'}")
    print(f"wrote {args.output_dir / 'team12_uncoupled.png'}")


if __name__ == "__main__":
    main()
