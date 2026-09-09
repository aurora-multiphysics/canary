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


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true",
                   help="report how faithfully the expressions follow the specification")
    args = p.parse_args()
    if args.check:
        check()
        return
    h, d = expressions()
    print("# H(B), A/m -- paste as the expression of the h_of_b function")
    print(h)
    print()
    print("# dH/dB, A/m/T -- paste as the expression of the dh_db function")
    print(d)


if __name__ == "__main__":
    main()
