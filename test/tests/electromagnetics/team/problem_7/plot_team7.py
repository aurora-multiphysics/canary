#!/usr/bin/env python3
"""Compare the computed fields with the TEAM Problem 7 measurements.

Reads the sampled values written by team7.i at 50 Hz and 200 Hz, prints them
against Tables 2 and 3 of the problem statement, and draws the figures embedded
in doc/content/problems/team_problem_7.md:

    team7_bz.png   B_z along A1-B1 and A2-B2, 15 mm above the plate
    team7_jy.png   J_y along A3-B3 and A4-B4, on the top and bottom of the plate

The solve is in the frequency domain with fields varying as exp(i omega t), and
the coil current is real, so it peaks at omega t = 0. The field at an instant is
the real part of the phasor times exp(i omega t), which makes the benchmark's two
columns

    omega t = 0      Re(phasor)
    omega t = 90     Re(i * phasor) = -Im(phasor)

The test runs the 50 Hz case, and a second, heavy, test the 200 Hz one. Produce
both by hand with

    <app> -i team7.i
    <app> -i team7.i frequency=200 Outputs/ReportedPostprocessors/file_base=OutputData/TEAM7CSV_200Hz

Usage:
    ./plot_team7.py                  # print the comparison and draw the figures
    ./plot_team7.py --no-plot        # print the comparison only
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_BASES = {50: HERE / "OutputData" / "TEAM7CSV",
                 200: HERE / "OutputData" / "TEAM7CSV_200Hz"}
DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the TEAM Problem 4 and 13 figures.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
MEASURED = "#e8710a"
COMPUTED = "#2a78d6"

GAUSS = 1.0e-4  # T
MA_PER_M2 = 1.0e6  # A/m2, the unit Table 3 is given in

# Tables 2 and 3 of the problem statement, the measured values: for each line, the
# positions along it in mm and, at each of 50 Hz and 200 Hz, the omega t = 0 and
# omega t = 90 degree columns. B_z is in gauss and J_y in 1e6 A/m2.
BZ_X = [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0, 180.0, 198.0,
        216.0, 234.0, 252.0, 270.0, 288.0]
JY_X = [0.0, 18.0, 126.0, 144.0, 162.0, 180.0, 198.0, 216.0, 234.0, 252.0, 270.0,
        288.0]

MEASURED_VALUES = {
    "bz_a1_b1": {
        50: ([-4.90, -17.88, -22.13, -20.19, -15.67, 0.36, 43.64, 78.11, 71.55, 60.44,
              53.91, 52.62, 53.81, 56.91, 59.24, 52.78, 27.61],
             [-1.16, 2.84, 4.15, 4.00, 3.07, 2.31, 1.89, 4.97, 12.61, 14.15, 13.04,
              12.40, 12.05, 12.27, 12.66, 9.96, 2.36]),
        200: ([-3.63, -18.46, -23.62, -21.59, -16.09, 0.23, 44.35, 75.53, 63.42, 53.20,
               48.66, 47.31, 48.31, 51.26, 53.61, 46.11, 24.96],
              [-1.38, 1.20, 2.15, 1.63, 1.10, 0.27, -2.28, -1.40, 4.17, 3.94, 4.86,
               4.09, 3.69, 4.60, 3.48, 4.10, 0.98]),
    },
    "bz_a2_b2": {
        50: ([-1.83, -8.50, -13.60, -15.21, -14.48, -5.62, 28.77, 60.34, 61.84, 56.64,
              53.40, 52.36, 53.93, 56.82, 59.48, 52.08, 26.56],
             [-1.63, -0.60, -0.43, 0.11, 1.26, 3.40, 6.53, 10.25, 11.83, 11.83, 11.01,
              10.58, 10.80, 10.54, 10.62, 9.03, 1.79]),
        200: ([-0.86, -7.00, -11.58, -13.36, -13.77, -6.74, 24.63, 53.19, 54.89, 50.72,
               48.03, 47.13, 48.25, 51.35, 53.35, 45.37, 24.01],
              [-1.35, -0.71, -0.81, -0.67, 0.15, 1.39, 2.67, 3.00, 4.01, 3.80, 4.00,
               3.02, 2.20, 2.78, 1.58, 1.37, 0.93]),
    },
    # The two J_y tables are assigned the other way round from the headings they
    # carry in the compumag.org copy of the problem statement, which puts the first
    # on the top face (A3-B3, z = 19 mm) and the second on the bottom (A4-B4, z = 0).
    # Those headings cannot be right: the second table has much the larger currents,
    # at 200 Hz three times larger at x = 288 mm and twenty times larger mid-plate,
    # and a plate three skin depths thick carries its larger currents on the face
    # nearer the coil, which is the top. The published comparisons agree: the
    # z = 19 mm measurements plotted in, for example, Kapidani et al.,
    # arXiv:2310.14313, Fig. 8, are the second table's.
    "jy_a3_b3": {
        50: ([0.461, 0.621, 1.573, 0.556, 0.237, 0.097, -0.034, -0.157, -0.305, -0.478,
              -0.660, -1.217],
             [-0.662, -0.644, -1.027, -0.757, -0.364, -0.149, 0.015, 0.154, 0.311, 0.508,
              0.747, 1.034]),
        200: ([1.057, 1.597, 4.163, 1.143, 0.672, 0.307, -0.050, -0.370, -0.749, -1.205,
               -1.575, -2.583],
              [-0.915, -1.036, -2.328, -1.193, -0.613, -0.259, 0.061, 0.334, 0.674, 1.064,
               1.404, 2.331]),
    },
    "jy_a4_b4": {
        50: ([0.249, 0.685, -0.015, -0.103, -0.061, -0.004, 0.051, 0.095, 0.135, 0.104,
              -0.321, -0.687],
             [-0.629, -0.873, -0.593, -0.249, -0.101, -0.001, 0.087, 0.182, 0.322, 0.555,
              0.822, 0.855]),
        200: ([0.427, 0.794, 1.401, -0.035, 0.005, -0.011, 0.007, 0.027, 0.042, 0.043,
               0.050, -0.321],
              [-0.623, -0.755, -1.304, -0.229, -0.041, -0.014, -0.002, -0.000, 0.008,
               0.033, 0.116, 0.893]),
    },
}

# Each line: the sampler reporting it, its title, the column of the sampled vector
# that holds the component the benchmark measures, the positions along it, and the
# factor that turns the sampled value into the table's unit.
LINES = [
    ("bz_a1_b1", "B$_z$ along A1-B1 (y = 72 mm, z = 34 mm)", "b_field", 2, BZ_X, GAUSS),
    ("bz_a2_b2", "B$_z$ along A2-B2 (y = 144 mm, z = 34 mm)", "b_field", 2, BZ_X, GAUSS),
    ("jy_a3_b3", "J$_y$ along A3-B3 (y = 72 mm, z = 19 mm)", "j_field", 1, JY_X, MA_PER_M2),
    ("jy_a4_b4", "J$_y$ along A4-B4 (y = 72 mm, z = 0 mm)", "j_field", 1, JY_X, MA_PER_M2),
]


def read_line(base, sampler, variable, component, unit):
    """Return the omega t = 0 and 90 degree values along one line, in table units."""
    path = Path(f"{base}_{sampler}_0001.csv")
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError:
        return None
    re_col, im_col = f"{variable}_real_{component}", f"{variable}_imag_{component}"
    return ([float(r[re_col]) / unit for r in rows],
            [-float(r[im_col]) / unit for r in rows])


def compare(results):
    """Print the computed values against the measured ones.

    Most of these profiles pass through zero, so each difference is scaled by the
    largest measured magnitude along that line, over both instants, rather than by
    the value at the point. Scaling by the peak of each column separately would
    make the omega t = 90 degree columns, which at 200 Hz are an order of magnitude
    smaller than the omega t = 0 ones, look far worse than they are.
    """
    for freq, lines in results.items():
        for sampler, title, _, _, xs, _ in LINES:
            computed = lines.get(sampler)
            if computed is None:
                continue
            unit = "G" if sampler.startswith("bz") else "1e6 A/m2"
            print(f"\n{title.replace('$', '')} at {freq} Hz  [{unit}]")
            print(f"{'x (mm)':>7} | {'wt=0 meas':>10} {'computed':>10} | "
                  f"{'wt=90 meas':>10} {'computed':>10}")
            meas0, meas90 = MEASURED_VALUES[sampler][freq]
            for i, x in enumerate(xs):
                print(f"{x:7.1f} | {meas0[i]:10.3f} {computed[0][i]:10.3f} | "
                      f"{meas90[i]:10.3f} {computed[1][i]:10.3f}")
            peak = max(abs(v) for v in meas0 + meas90)
            for label, meas, comp in (("wt=0", meas0, computed[0]),
                                      ("wt=90", meas90, computed[1])):
                worst = max(abs(c - m) for c, m in zip(comp, meas)) / peak
                rms = (sum((c - m) ** 2 for c, m in zip(comp, meas)) / len(meas)) ** 0.5 / peak
                print(f"  {label:>5}: rms difference {100 * rms:5.1f}% and largest "
                      f"{100 * worst:5.1f}% of the peak measured value on the line")


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
    ax.axhline(0.0, color=TEXT_SECONDARY, linewidth=0.6)


def plot(results, samplers, filename, ylabel, outdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    freqs = sorted(results)
    fig, axes = plt.subplots(2 * len(samplers), len(freqs), figsize=(5.2 * len(freqs), 10.5),
                             squeeze=False, facecolor=SURFACE)
    for col, freq in enumerate(freqs):
        for i, sampler in enumerate(samplers):
            _, title, _, _, xs, _ = next(line for line in LINES if line[0] == sampler)
            meas = MEASURED_VALUES[sampler][freq]
            comp = results[freq].get(sampler)
            for j, phase in enumerate(("0", "90")):
                ax = axes[2 * i + j][col]
                _style(ax)
                ax.plot(xs, meas[j], linestyle="none", marker="o", markersize=5,
                        color=MEASURED, label="measured")
                if comp is not None:
                    # The J_y lines cross the hole between x = 18 and 126 mm, where
                    # there is no current, so the curve is broken there rather than
                    # drawn straight across it.
                    px, py = [], []
                    for k, (x, v) in enumerate(zip(xs, comp[j])):
                        if k and x - xs[k - 1] > 50.0:
                            px.append(float("nan"))
                            py.append(float("nan"))
                        px.append(x)
                        py.append(v)
                    ax.plot(px, py, color=COMPUTED, linewidth=1.6, marker=".",
                            label="MFEM-MOOSE")
                ax.set_title(f"{title}, {freq} Hz, $\\omega t$ = {phase}°",
                             fontsize=9.5, color=TEXT_PRIMARY, loc="left")
                ax.set_ylabel(ylabel, fontsize=9, color=TEXT_SECONDARY)
                if 2 * i + j == 2 * len(samplers) - 1:
                    ax.set_xlabel("x (mm)", fontsize=9, color=TEXT_SECONDARY)
    axes[0][0].legend(frameon=False, fontsize=9, labelcolor=TEXT_SECONDARY)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / filename
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-50", type=Path, default=DEFAULT_BASES[50],
                   help="file_base of the 50 Hz run")
    p.add_argument("--base-200", type=Path, default=DEFAULT_BASES[200],
                   help="file_base of the 200 Hz run")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--no-plot", action="store_true", help="print the comparison only")
    args = p.parse_args()

    results = {}
    for freq, base in ((50, args.base_50), (200, args.base_200)):
        lines = {s: read_line(base, s, var, comp, unit) for s, _, var, comp, _, unit in LINES}
        if any(v is not None for v in lines.values()):
            results[freq] = lines
    if not results:
        sys.exit("error: no results found. Run team7.i first; it writes OutputData/TEAM7CSV_*.csv.")

    compare(results)
    if not args.no_plot:
        plot(results, ["bz_a1_b1", "bz_a2_b2"], "team7_bz.png", "B$_z$ (G)", args.outdir)
        plot(results, ["jy_a3_b3", "jy_a4_b4"], "team7_jy.png", "J$_y$ (10$^6$ A/m$^2$)",
             args.outdir)


if __name__ == "__main__":
    main()
