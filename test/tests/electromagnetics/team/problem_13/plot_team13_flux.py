#!/usr/bin/env python3
"""Plot the computed flux density in the steel against the measurements.

Reads the postprocessor values written by team13.i and produces the two figures
embedded in doc/content/problems/team_problem_13.md:

    team13_flux.png       the profile up the centre plate, positions 1 to 7
    team13_positions.png  every one of the 25 sections of Table 1

Both show the 1000 AT case, which is the one the problem statement publishes
measurements for, and the second order Nedelec solution of it alongside. The
benchmark's other excitation, 3000 AT, has no measured column to be read against,
so it is left off unless --csv-3000 asks for it.

The measurements were taken on two specimens whose plate-to-channel gap came out at
0.52 mm and 0.47 mm against a nominal 0.5 mm, and that gap carries most of the
reluctance of the magnetic circuit. The two sets of readings are therefore drawn as
a band rather than as a single curve: the width of the band is the experimental
spread the model is being asked to land inside.

The 3000 AT results are not written by the test, which runs the 1000 AT case only.
Produce them with

    <app> -i team13.i ampere_turns=3000 \
          Outputs/ReportedPostprocessors/file_base=OutputData/TEAM13CSV_3000AT

Usage:
    ./plot_team13_flux.py                       # after running the cases above
    ./plot_team13_flux.py --csv other.csv --outdir /tmp
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_CSV = HERE / "OutputData" / "TEAM13CSV.csv"
DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the TEAM Problem 4 figure.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
MEASURED = "#e8710a"
COMPUTED = "#2a78d6"        # the base case: 1000 AT, first order
COMPUTED_SECOND = "#7b4fd1"  # second order Nedelec
COMPUTED_3000 = "#12866e"    # 3000 AT

# Table 1 of the problem statement: the 25 sections, in the order it numbers them,
# with the postprocessor reporting each and the flux density measured on the two
# specimens (gaps G = 0.52 mm and G = 0.47 mm) at 1000 AT.
#
#   number, postprocessor, part, position along the part, G=0.52, G=0.47
POSITIONS = [
    (1, "PlateB000", "plate", 0.0, 1.333, 1.354),
    (2, "PlateB010", "plate", 10.0, 1.329, 1.339),
    (3, "PlateB020", "plate", 20.0, 1.286, 1.304),
    (4, "PlateB030", "plate", 30.0, 1.225, 1.245),
    (5, "PlateB040", "plate", 40.0, 1.129, 1.138),
    (6, "PlateB050", "plate", 50.0, 0.985, 0.982),
    (7, "PlateB060", "plate", 60.0, 0.655, 0.674),
    (8, "FlangeB002", "flange", 2.1, 0.259, 0.263),
    (9, "FlangeB010", "flange", 10.0, 0.453, 0.451),
    (10, "FlangeB020", "flange", 20.0, 0.554, 0.563),
    (11, "FlangeB030", "flange", 30.0, 0.637, 0.641),
    (12, "FlangeB040", "flange", 40.0, 0.698, 0.706),
    (13, "FlangeB050", "flange", 50.0, 0.755, 0.763),
    (14, "FlangeB060", "flange", 60.0, 0.809, 0.819),
    (15, "FlangeB080", "flange", 80.0, 0.901, 0.907),
    (16, "FlangeB100", "flange", 100.0, 0.945, 0.958),
    (17, "FlangeB110", "flange", 110.0, 0.954, 0.968),
    (18, "FlangeB122", "flange", 122.1, 0.956, 0.968),
    (19, "WebB060", "web", 60.0, 0.960, 0.971),
    (20, "WebB050", "web", 50.0, 0.965, 0.973),
    (21, "WebB040", "web", 40.0, 0.970, 0.982),
    (22, "WebB030", "web", 30.0, 0.974, 0.985),
    (23, "WebB020", "web", 20.0, 0.981, 0.991),
    (24, "WebB010", "web", 10.0, 0.984, 0.995),
    (25, "WebB000", "web", 0.0, 0.985, 0.995),
]

# The runs of consecutive positions that follow one part of the magnetic circuit.
# Positions are only joined up within a run: 7 to 8 steps from the top of the plate
# to the tip of the flange, which is a jump across the gap, not a profile.
PARTS = [
    ("plate", "centre plate", "height up the plate, z (mm)"),
    ("flange", "top flange", "along the flange, x (mm)"),
    ("web", "web", "height up the web, z (mm)"),
]


def read_values(path, required=True):
    """Return the last row of the postprocessor CSV as a dict of floats.

    An optional case that has not been run yet returns None, and is left off the
    figures rather than stopping them being drawn.
    """
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError:
        if not required:
            return None
        sys.exit(
            f"error: {path} not found.\nRun team13.i first; it writes OutputData/TEAM13CSV.csv."
        )
    if not rows:
        sys.exit(f"error: {path} holds no rows")

    missing = [p[1] for p in POSITIONS if p[1] not in rows[-1]]
    if missing:
        sys.exit(f"error: {path} is missing column(s): {', '.join(missing)}")
    return {k: float(v) for k, v in rows[-1].items() if v != ""}


def _suffix(cases):
    """Name the excitation in the title when there is only one to name."""
    amps = {c[0] for c in cases}
    return f", at {amps.pop()} AT" if len(amps) == 1 else ""


def _top(cases):
    """A common y limit for both figures, so they can be read against each other."""
    peak = max(max(v for k, v in values.items() if k != "time") for _, _, values, _, _ in cases)
    peak = max(peak, max(max(p[4], p[5]) for p in POSITIONS))
    return 0.1 * (int(peak / 0.1) + 3)


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, solid_capstyle="butt")
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=0)


def _series(ax, xs, rows, cases, label=False):
    """Draw the measured band and one computed curve per excitation.

    cases is a list of (ampere turns, legend label, values, colour, marker); labels
    are attached only on the first call so that a figure split into several runs of
    positions still gets one legend entry per series.
    """
    lower = [min(r[4], r[5]) for r in rows]
    upper = [max(r[4], r[5]) for r in rows]
    ax.fill_between(xs, lower, upper, color=MEASURED, alpha=0.18, linewidth=0)
    ax.plot(xs, [0.5 * (r[4] + r[5]) for r in rows], color=MEASURED, linewidth=2.0,
            marker="o", markersize=7, markeredgecolor=SURFACE, markeredgewidth=2.0,
            label="Measured, 1000 AT (band spans the two specimens)" if label else None)
    for _, series_label, values, colour, marker in cases:
        ax.plot(xs, [values[r[1]] for r in rows], color=colour, linewidth=2.0,
                marker=marker, markersize=7, markeredgecolor=SURFACE, markeredgewidth=2.0,
                label=series_label if label else None)


def plot_plate(cases, output):
    """The profile up the centre plate, positions 1 to 7."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    rows = [p for p in POSITIONS if p[2] == "plate"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    _series(ax, [r[3] for r in rows], rows, cases, label=True)

    ax.set_xlabel("Height up the plate, z (mm)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Flux density through the section (T)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_title("Flux density in the centre plate, positions 1 to 7" + _suffix(cases),
                 color=TEXT_PRIMARY, fontsize=12, loc="left", pad=12)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlim(-2, 62)
    ax.set_ylim(0, _top(cases))
    legend = ax.legend(loc="lower left", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def plot_positions(cases, output):
    """All 25 sections of Table 1, against the position number the benchmark gives."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)

    top = _top(cases)
    for i, (key, title, _) in enumerate(PARTS):
        rows = [p for p in POSITIONS if p[2] == key]
        xs = [r[0] for r in rows]
        _series(ax, xs, rows, cases, label=(i == 0))
        # A rule between the parts: consecutive numbers either side of it are on
        # opposite sides of an air gap, so the curves are not joined across it.
        if i:
            ax.axvline(xs[0] - 0.5, color=GRID, linewidth=1.0, zorder=0)
        ax.text(0.5 * (xs[0] + xs[-1]), top * 0.97, title, ha="center", va="top",
                color=TEXT_SECONDARY, fontsize=9)

    # Label only the two that need explaining; the axes carry every other value.
    # Both sit in the empty space under the flange curve, with a leader back to the
    # point, so neither runs into the other or into the data.
    notes = [
        (7, "averages over the turn\nout into the flange", (4.6, 0.26)),
        (8, "band average, not the\nsection at the tip", (10.4, 0.14)),
    ]
    for number, note, xytext in notes:
        row = next(p for p in POSITIONS if p[0] == number)
        ax.annotate(note, xy=(number, cases[0][2][row[1]]), xytext=xytext,
                    textcoords="data", color=TEXT_SECONDARY, fontsize=8,
                    ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", color=GRID, linewidth=1.2,
                                    shrinkA=2, shrinkB=6))

    ax.set_xlabel("Position number in Table 1 of the problem statement",
                  color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Flux density through the section (T)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_title("Flux density in the steel, all 25 measured sections" + _suffix(cases),
                 color=TEXT_PRIMARY, fontsize=12, loc="left", pad=12)
    ax.set_xticks([p[0] for p in POSITIONS])
    ax.set_xticklabels([str(p[0]) for p in POSITIONS], fontsize=8)
    ax.set_xlim(0.3, 25.7)
    ax.set_ylim(0, top)
    legend = ax.legend(loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def report(cases):
    header = "  %-4s %-8s %-9s" % ("No.", "part", "position")
    for _, series_label, _, _, _ in cases:
        header += " %13s" % series_label.replace("Computed, ", "").replace("Computed", "computed")
    header += " %9s %9s %8s" % ("G=0.52", "G=0.47", "error")
    print(header)
    worst = (0.0, None)
    for number, key, part, where, g52, g47 in POSITIONS:
        mean = 0.5 * (g52 + g47)
        err = 100.0 * (cases[0][2][key] / mean - 1.0)
        if abs(err) > abs(worst[0]):
            worst = (err, number)
        line = "  %-4d %-8s %-9s" % (number, part, f"{where:g} mm")
        for _, _, values, _, _ in cases:
            line += " %13.4f" % values[key]
        print(line + " %9.3f %9.3f %7.1f%%" % (g52, g47, err))
    print("  %s: worst departure from the measured mean %+.1f%% at position %d"
          % (cases[0][1], worst[0], worst[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV,
                        help=f"1000 AT postprocessor CSV (default: {DEFAULT_CSV})")
    parser.add_argument("--csv-3000", type=Path, default=None,
                        help="3000 AT postprocessor CSV; drawn alongside 1000 AT when given")
    parser.add_argument("--csv-second-order", type=Path, default=None,
                        help="second order Nedelec postprocessor CSV; drawn alongside when given")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                        help=f"directory for the figures (default: {DEFAULT_OUTDIR})")
    parser.add_argument("--no-plot", action="store_true", help="only print the comparison")
    args = parser.parse_args()

    second = read_values(args.csv_second_order) if args.csv_second_order else None
    at3000 = read_values(args.csv_3000) if args.csv_3000 else None
    # Name the base series for whatever it is being contrasted with.
    base = "Computed"
    if second is not None:
        base = "Computed, first order"
    elif at3000 is not None:
        base = "Computed, 1000 AT"
    cases = [(1000, base, read_values(args.csv), COMPUTED, "s")]
    if second is not None:
        cases.append((1000, "Computed, second order", second, COMPUTED_SECOND, "^"))
    if at3000 is not None:
        cases.append((3000, "Computed, 3000 AT", at3000, COMPUTED_3000, "D"))

    if not args.no_plot:
        plot_plate(cases, args.outdir / "team13_flux.png")
        plot_positions(cases, args.outdir / "team13_positions.png")
        print(f"wrote {args.outdir / 'team13_flux.png'}")
        print(f"wrote {args.outdir / 'team13_positions.png'}")
    report(cases)


if __name__ == "__main__":
    main()
