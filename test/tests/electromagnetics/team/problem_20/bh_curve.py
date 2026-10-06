#!/usr/bin/env python3
"""Turn the TEAM Problem 20 B-H curve into expressions for the MOOSE input.

The benchmark specifies the steel by the 38 (B, H) points of its Table 1, up to
B = 2.3 T and H = 135000 A/m, and above that by its Eq. (1),

    B = mu0 H + Ms,

with Ms the saturation magnetisation, quoted as 2.16 T.  Table 6 of the problem
statement asks which interpolation was used between the tabulated points and
lists splines first, so that is what is done here: a monotone cubic
(Fritsch-Carlson) through all 38 points, joined at B = 2.3 T to Eq. (1).

The solver needs the reluctivity nu = H/B and, for the Jacobian, B dnu/dB.  Both
follow from H(B) and its derivative,

    nu = H / B,        B dnu/dB = H'(B) - H/B,

so this script emits H(B) and H'(B) and the input file forms the rest.  Run it
with no arguments to print the two expressions, or --check to report how well
they reproduce the specification.

One liberty is taken with the specification, and --check reports it: Eq. (1)
with Ms = 2.16 T does not pass through the last tabulated point.  At H = 135000
A/m it gives B = 2.3296 T where the table says 2.3, a step of 30 mT.  Ms is
instead taken as 2.3 - mu0 * 135000 = 2.1304 T, so that the line continues from
the last point, and the spline's slope there is clamped to 1/mu0 so that the
curve is continuously differentiable across the join.
"""

import argparse
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the other TEAM figures.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
CURVE = "#2a78d6"

MU0 = 4.0e-7 * math.pi
MS_QUOTED = 2.16

# B (T), H (A/m) from Table 1 of the problem statement.
TABLE = [
    (0.0, 0.0), (0.01, 27.0), (0.025, 58.0), (0.05, 100.0), (0.10, 153.0),
    (0.15, 185.0), (0.20, 205.0), (0.30, 233.0), (0.40, 255.0), (0.50, 285.0),
    (0.60, 320.0), (0.70, 355.0), (0.80, 405.0), (0.90, 470.0), (1.00, 555.0),
    (1.10, 673.0), (1.20, 836.0), (1.30, 1065.0), (1.35, 1220.0), (1.40, 1420.0),
    (1.45, 1720.0), (1.50, 2130.0), (1.55, 2670.0), (1.60, 3480.0), (1.65, 4500.0),
    (1.70, 5950.0), (1.75, 7650.0), (1.80, 10100.0), (1.85, 13000.0),
    (1.90, 15900.0), (1.95, 21100.0), (2.00, 26300.0), (2.05, 32900.0),
    (2.10, 42700.0), (2.15, 61700.0), (2.20, 84300.0), (2.25, 110000.0),
    (2.30, 135000.0),
]

B_LAST, H_LAST = TABLE[-1]
M_SAT = B_LAST - MU0 * H_LAST  # 2.1304 T: Eq. (1) through the last point


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
    m[-1] = min(1.0 / MU0, 3.0 * d[-1])
    return m


def _n(v):
    return "%.7g" % v


def _segments():
    bs = [b for b, _ in TABLE]
    hs = [h for _, h in TABLE]
    m = slopes(bs, hs)
    for i in range(len(bs) - 1):
        w = bs[i + 1] - bs[i]
        c2 = (3 * (hs[i + 1] - hs[i]) / w - 2 * m[i] - m[i + 1]) / w
        c3 = (m[i] + m[i + 1] - 2 * (hs[i + 1] - hs[i]) / w) / (w * w)
        yield bs[i], bs[i + 1], hs[i], m[i], c2, c3


def expressions():
    """Return the H(B) and dH/dB expressions, in terms of a symbol named B."""
    h_segs, d_segs = [], []
    for b0, b1, h0, m0, c2, c3 in _segments():
        t = "B" if b0 == 0.0 else "(B-%s)" % _n(b0)
        h_segs.append((b1, "%s+%s*(%s+%s*(%s+%s*%s))"
                       % (_n(h0), t, _n(m0), t, _n(c2), t, _n(c3))))
        d_segs.append((b1, "%s+%s*(%s+%s*%s)" % (_n(m0), t, _n(2 * c2), t, _n(3 * c3))))

    h_tail = "(B-%s)/%s" % (_n(M_SAT), _n(MU0))
    d_tail = _n(1.0 / MU0)

    def wrap(segs, tail):
        e = tail
        for bmax, seg in reversed(segs):
            e = "if(B<%s,%s,%s)" % (_n(bmax), seg, e)
        return e

    return wrap(h_segs, h_tail), wrap(d_segs, d_tail)


# --- evaluation, used only by --check and --plot ----------------------------
def h_of_b(b):
    for b0, b1, h0, m0, c2, c3 in _segments():
        if b < b1:
            t = b - b0
            return h0 + t * (m0 + t * (c2 + t * c3))
    return (b - M_SAT) / MU0


def check():
    print("Eq.(1) with Ms = %.2f T gives B = %.4f T at H = %.0f A/m; the table says %.2f"
          % (MS_QUOTED, MU0 * H_LAST + MS_QUOTED, H_LAST, B_LAST))
    print("so Ms is taken as %.4f T, which puts Eq.(1) through the last point\n" % M_SAT)

    worst = max(abs(h_of_b(b) - h) for b, h in TABLE[:-1])
    print("tabulated points reproduced to %.3g A/m" % worst)

    n, prev, ok = 400000, -1.0, True
    for i in range(n + 1):
        v = h_of_b(2.6 * i / n)
        if v < prev:
            ok = False
        prev = v
    print("H(B) monotone increasing on [0, 2.6]: %s" % ok)

    eps = 1e-7
    left = (h_of_b(B_LAST - eps) - h_of_b(B_LAST - 2 * eps)) / eps
    print("slope either side of the join at %.1f T: %.0f and %.0f A/m/T (1/mu0 = %.0f)"
          % (B_LAST, left, 1.0 / MU0, 1.0 / MU0))

    nu = [(h_of_b(b) / b, b) for b in [0.001 * k for k in range(1, 2400)]]
    lo = min(nu)
    print("\nreluctivity nu = H/B runs from %.0f as B -> 0, down to %.0f at B = %.2f T,"
          % (h_of_b(1e-6) / 1e-6, lo[0], lo[1]))
    print("and up to %.0f at B = 2.3 T" % (H_LAST / B_LAST))


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


def plot(output):
    """B against H on a log field axis, over the whole range through to saturation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    fig.suptitle("B-H curve of the TEAM Problem 20 steel", color=TEXT_PRIMARY,
                 fontsize=12, x=0.012, ha="left", y=0.975)
    ax.set_xscale("log")

    n, b_max = 6000, 2.6
    spline = [(h_of_b(b_max * i / n), b_max * i / n) for i in range(1, n + 1)
              if b_max * i / n <= B_LAST]
    tail = [(h_of_b(b), b) for b in (B_LAST, b_max)]
    ax.plot([h for h, _ in spline], [b for _, b in spline], color=CURVE, linewidth=2.0)
    ax.plot([h for h, _ in tail], [b for _, b in tail], color=CURVE, linewidth=2.0,
            linestyle=(0, (4, 2)))
    ax.plot([h for _, h in TABLE[1:]], [b for b, _ in TABLE[1:]], linestyle="none",
            marker="o", markersize=5.0, color=CURVE, markeredgecolor="none")
    ax.set_xlim(10, 1.0e6)
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("Magnetic field strength, H (A/m)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Flux density, B (T)", color=TEXT_SECONDARY, fontsize=10)
    fig.text(0.988, 0.915,
             r"$B = \mu_0 H + M_s$ above B = %.1f T (dashed), $M_s$ = %.4f T"
             % (B_LAST, M_SAT),
             color=TEXT_SECONDARY, fontsize=8.5, ha="right", va="top")
    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=5.0, color=CURVE,
               markeredgecolor="none", label="Table 1 of the problem statement"),
        Line2D([], [], color=CURVE, linewidth=2.0, label="Monotone cubic spline"),
    ]
    legend = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)

    fig.tight_layout(rect=(0, 0, 1, 0.9))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


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
        out = args.outdir / "team20_bh_curve.png"
        plot(out)
        print(f"wrote {out}")
        return
    h, d = expressions()
    print("# H(B), A/m -- paste as the expression of the h_of_b function")
    print(h)
    print()
    print("# dH/dB, A/m/T -- paste as the expression of the dh_db function")
    print(d)


if __name__ == "__main__":
    main()
