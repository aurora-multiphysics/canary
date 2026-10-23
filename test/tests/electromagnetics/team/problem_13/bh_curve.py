#!/usr/bin/env python3
"""Turn the TEAM Problem 13 B-H curve into expressions for the MOOSE input.

The benchmark specifies the steel by 27 tabulated (B, H) points up to B = 1.8 T
and, above that, by the analytic form of its Eq. (1):

    B = mu0 H + (a H^2 + b H + c)   for 1.8 <= B <= 2.22 T
    B = mu0 H + Ms                  above that

with a = -2.381e-10, b = 2.327e-5, c = 1.590 and Ms = 2.16 T.  Table 4 of the
problem statement asks which interpolation was used between the tabulated
points and lists splines as the first option, so that is what is done here: a
monotone cubic (Fritsch-Carlson) through all 27 points, joined at B = 1.8 T to
Eq. (1) inverted for H(B).

The solver needs the reluctivity nu = H/B and, for the Jacobian, B dnu/dB.  Both
follow from H(B) and its derivative,

    nu = H / B,        B dnu/dB = H'(B) - H/B,

so this script emits H(B) and H'(B) and the input file forms the rest.  Run it
with no arguments to print the two expressions, or --check to report how well
they reproduce the specification.

Two small liberties are taken with the specification, both noted by --check:

* Eq. (1) gives H = 9424.35 A/m at B = 1.8 T where the table says 9423, so the
  last tabulated point is moved onto Eq. (1) to make the join continuous.
* The parabola a H^2 + b H + c peaks at 2.158556 T rather than at the quoted
  Ms = 2.16 T.  Taking that peak as the saturation magnetisation keeps B(H)
  continuously differentiable where the two branches meet; using 2.16 exactly
  would put a 1.4 mT step in the curve and make H(B) non-monotone.
"""

import argparse
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the other TEAM Problem 13 figures; see
# plot_team13_flux.py.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
# The tabulated points and the curve through them are one series, so they share a
# colour; the points are drawn without a surface-coloured ring, which would break
# the line wherever they sit close together.
CURVE = "#2a78d6"

MU0 = 4.0e-7 * math.pi
A, B_, C = -2.381e-10, 2.327e-5, 1.590

# B (T), H (A/m) from Fig. 3 of the problem statement.
TABLE = [
    (0.0, 0.0), (0.0025, 16.0), (0.0050, 30.0), (0.0125, 54.0), (0.025, 93.0),
    (0.05, 143.0), (0.10, 191.0), (0.20, 210.0), (0.30, 222.0), (0.40, 233.0),
    (0.50, 247.0), (0.60, 258.0), (0.70, 272.0), (0.80, 289.0), (0.90, 313.0),
    (1.00, 342.0), (1.10, 377.0), (1.20, 433.0), (1.30, 509.0), (1.40, 648.0),
    (1.50, 933.0), (1.55, 1228.0), (1.60, 1934.0), (1.65, 2913.0),
    (1.70, 4993.0), (1.75, 7189.0), (1.80, 9423.0),
]

H_PEAK = -B_ / (2 * A)                 # where the magnetisation parabola peaks
M_SAT = C - B_ * B_ / (4 * A)          # its peak value, 2.158556 T
B_SAT = MU0 * H_PEAK + M_SAT           # 2.219963 T, quoted as 2.22 in the statement
S = B_ + MU0


def h_upper(b):
    """Eq. (1) inverted for H, in the form that avoids cancellation."""
    disc = max(S * S - 4 * A * (C - b), 0.0)
    return 2.0 * (C - b) / (-S - math.sqrt(disc))


def knots():
    """Table points, with the last moved onto Eq. (1) so the join is continuous."""
    bs = [b for b, _ in TABLE]
    hs = [h for _, h in TABLE]
    hs[-1] = h_upper(bs[-1])
    return bs, hs


def slopes(bs, hs):
    """Fritsch-Carlson slopes, clamped at the right end onto Eq. (1)'s slope."""
    n = len(bs)
    h = [bs[i + 1] - bs[i] for i in range(n - 1)]
    d = [(hs[i + 1] - hs[i]) / h[i] for i in range(n - 1)]
    m = [0.0] * n
    for i in range(1, n - 1):
        if d[i - 1] * d[i] > 0.0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    m[0] = d[0]
    m[-1] = min(1.0 / (MU0 + 2 * A * hs[-1] + B_), 3.0 * d[-1])
    return m


def _n(v):
    return "%.7g" % v


def expressions():
    """Return the H(B) and dH/dB expressions, in terms of a symbol named B."""
    bs, hs = knots()
    m = slopes(bs, hs)
    h_segs, d_segs = [], []
    for i in range(len(bs) - 1):
        w = bs[i + 1] - bs[i]
        c2 = (3 * (hs[i + 1] - hs[i]) / w - 2 * m[i] - m[i + 1]) / w
        c3 = (m[i] + m[i + 1] - 2 * (hs[i + 1] - hs[i]) / w) / (w * w)
        t = "B" if bs[i] == 0.0 else "(B-%s)" % _n(bs[i])
        h_segs.append((bs[i + 1], "%s+%s*(%s+%s*(%s+%s*%s))"
                       % (_n(hs[i]), t, _n(m[i]), t, _n(c2), t, _n(c3))))
        d_segs.append((bs[i + 1], "%s+%s*(%s+%s*%s)"
                       % (_n(m[i]), t, _n(2 * c2), t, _n(3 * c3))))

    inv = "2*(%s-B)/(-%s-sqrt(max(%s-4*(%s)*(%s-B),0)))" % (_n(C), _n(S), _n(S * S), _n(A), _n(C))
    h_tail = "if(B>=%s,(B-%s)/%s,%s)" % (_n(B_SAT), _n(M_SAT), _n(MU0), inv)
    d_tail = "if(B>=%s,%s,1/(%s+2*(%s)*(%s)))" % (_n(B_SAT), _n(1.0 / MU0), _n(S), _n(A), inv)

    def wrap(segs, tail):
        e = tail
        for bmax, seg in reversed(segs):
            e = "if(B<%s,%s,%s)" % (_n(bmax), seg, e)
        return e

    return wrap(h_segs, h_tail), wrap(d_segs, d_tail)


# --- evaluation, used only by --check ---------------------------------------
def h_of_b(b):
    bs, hs = knots()
    m = slopes(bs, hs)
    for i in range(len(bs) - 1):
        if b < bs[i + 1]:
            w, t = bs[i + 1] - bs[i], b - bs[i]
            c2 = (3 * (hs[i + 1] - hs[i]) / w - 2 * m[i] - m[i + 1]) / w
            c3 = (m[i] + m[i + 1] - 2 * (hs[i + 1] - hs[i]) / w) / (w * w)
            return hs[i] + t * (m[i] + t * (c2 + t * c3))
    return (b - M_SAT) / MU0 if b >= B_SAT else h_upper(b)


def check():
    bs, hs = knots()
    print("Eq.(1) at B=1.8 gives H=%.2f A/m; the table says 9423 (%.3f%% apart)"
          % (hs[-1], 100 * (hs[-1] / 9423.0 - 1)))
    print("magnetisation parabola peaks at %.6f T (Ms is quoted as 2.16 T),"
          % M_SAT)
    print("which puts the fully saturated branch at B >= %.6f T (quoted as 2.22 T)\n" % B_SAT)

    worst = max(abs(h_of_b(b) - h) for b, h in TABLE[:-1])
    print("tabulated points reproduced to %.3g A/m" % worst)

    n, prev, ok = 400000, -1.0, True
    for i in range(n + 1):
        b = 2.6 * i / n
        v = h_of_b(b)
        if v < prev:
            ok = False
        prev = v
    print("H(B) monotone increasing on [0, 2.6]: %s" % ok)

    print("\n     B (T)      H spec      H here")
    for b, h in TABLE[1:]:
        print("   %7.4f  %10.1f  %10.1f" % (b, h, h_of_b(b)))
    for b in (1.85, 1.9, 2.0, 2.1, 2.2, B_SAT, 2.3, 2.4):
        ref = (b - M_SAT) / MU0 if b >= B_SAT else h_upper(b)
        print("   %7.4f  %10.1f  %10.1f   (Eq. 1)" % (b, ref, h_of_b(b)))

    nu = [h_of_b(b) / b for b in (0.001, 0.1, 0.5, 1.0, 1.05, 1.5, 1.8, 2.2)]
    print("\nreluctivity nu = H/B ranges from %.0f (at B ~ 1.05 T) to %.0f at B = 2.2 T"
          % (min(nu), nu[-1]))


def _sci(v):
    """A coefficient in mathtext scientific notation, sign included."""
    e = math.floor(math.log10(abs(v)))
    return r"%.4g{\times}10^{%d}" % (v / 10.0 ** e, e)


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


def _curve(b_max=2.6, n=6000):
    """The curve as (H, B) pairs, split where Eq. (1) takes over from the spline."""
    spline, analytic = [], []
    for i in range(n + 1):
        b = b_max * i / n
        point = (h_of_b(b), b)
        # 1.8 T is the last tabulated point; above it the curve is Eq. (1). The
        # join is repeated in both lists so the two segments meet on the plot.
        if b <= TABLE[-1][0]:
            spline.append(point)
        if b >= TABLE[-1][0]:
            analytic.append(point)
    return spline, analytic


def _draw(ax, spline, analytic=None):
    """The spline, the analytic branch above it if asked for, and the table points."""
    ax.plot([h for h, _ in spline], [b for _, b in spline], color=CURVE,
            linewidth=2.0, zorder=2)
    if analytic:
        ax.plot([h for h, _ in analytic], [b for _, b in analytic], color=CURVE,
                linewidth=2.0, linestyle=(0, (4, 2)), zorder=2)
    # The first tabulated point is H = 0, which has no place on a log axis.
    first = 1 if ax.get_xscale() == "log" else 0
    ax.plot([h for _, h in TABLE[first:]], [b for b, _ in TABLE[first:]],
            linestyle="none", marker="o", markersize=5.0, color=CURVE,
            markeredgecolor="none", zorder=3)
    ax.set_xlabel("Magnetic field strength, H (A/m)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Flux density, B (T)", color=TEXT_SECONDARY, fontsize=10)


def _legend(ax):
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=5.0, color=CURVE,
               markeredgecolor="none", label="Reference data"),
        Line2D([], [], color=CURVE, linewidth=2.0, label="Cubic spline fit"),
    ]
    legend = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)


def _new_figure(title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=12, x=0.012, ha="left", y=0.975)
    return fig, ax


def _save(fig, output, top):
    """Lay the figure out once everything is on it, then write it."""
    import matplotlib.pyplot as plt

    fig.tight_layout(rect=(0, 0, 1, top))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def plot_knee(output):
    """B against H on linear axes, over the range where the permeability varies.

    Carrying the axis to the last tabulated point instead would put twenty-two of
    the twenty-seven points in its leftmost tenth; plot_saturation covers the rest.
    """
    spline, _ = _curve()
    fig, ax = _new_figure("B-H curve of the TEAM Problem 13 steel, to the knee")
    _draw(ax, spline)
    ax.set_xlim(-25, 1000)
    ax.set_ylim(0, 1.55)
    _legend(ax)
    _save(fig, output, 0.93)


def plot_saturation(output):
    """B against H on a log field axis, over the whole range through to saturation."""
    spline, analytic = _curve()
    fig, ax = _new_figure("B-H curve of the TEAM Problem 13 steel, to saturation")
    ax.set_xscale("log")  # before _draw, which drops the H = 0 point for a log axis
    _draw(ax, spline, analytic)
    ax.set_xlim(10, 1.0e6)
    ax.set_ylim(0, 2.6)
    # Above the last tabulated point the curve is B = mu0 H + M(H), with the
    # magnetisation the benchmark's quadratic, held at its peak once it gets there.
    fig.text(0.988, 0.900,
             r"$M(H) = %sH^2 + %sH + %.3f$" % (_sci(A), _sci(B_), C) + "\n"
             r"above B = %.1f T (dashed), held at $M_s$ = %.4f T" % (TABLE[-1][0], M_SAT),
             color=TEXT_SECONDARY, fontsize=8.5, ha="right", va="top", linespacing=1.5)
    _legend(ax)
    _save(fig, output, 0.82)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true",
                   help="report how faithfully the expressions follow the specification")
    p.add_argument("--plot", action="store_true",
                   help="write a figure of the B-H curve instead of the expressions")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="directory the figure is written to (default: %(default)s)")
    args = p.parse_args()
    if args.check:
        check()
        return
    if args.plot:
        for name, draw in (("team13_bh_curve.png", plot_knee),
                           ("team13_bh_saturation.png", plot_saturation)):
            draw(args.outdir / name)
            print(f"wrote {args.outdir / name}")
        return
    h, d = expressions()
    print("# H(B), A/m -- paste as the expression of the h_of_b function")
    print(h)
    print()
    print("# dH/dB, A/m/T -- paste as the expression of the dh_db function")
    print(d)


if __name__ == "__main__":
    main()
