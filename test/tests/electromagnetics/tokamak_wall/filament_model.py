#!/usr/bin/env python3
"""Axisymmetric filament model of the tokamak wall induction example.

This is the independent reference that tokamak_wall.i is verified against; the
comparison is documented in doc/content/problems/tokamak_wall_induction.md.

When the plasma current and the wall are both axisymmetric, the induced current
is purely toroidal and the problem reduces to circuits: cut the wall's poloidal
cross-section into small cells, and each cell becomes a closed toroidal loop
carrying a current I_k, coupled to every other loop by mutual inductance and
linked by the plasma's poloidal flux. With no electrostatic potential (nothing
varies toroidally) Faraday's law around loop k is

    sum_j L_kj dI_j/dt + R_k I_k = -d(Psi_k)/dt,

with L the inductance matrix of the wall loops, R_k the toroidal resistance of
cell k, and Psi_k the flux of the plasma current through loop k. The mutual
inductance of two coaxial circular loops is Maxwell's formula in complete
elliptic integrals; the self inductance of a cell is that of a thin ring with
the cell's cross-section, through its geometric mean distance. As the cells are
refined this converges to the exact axisymmetric eddy current solution, and it
shares nothing with the finite element model: no mesh, no vector potential, no
far-field boundary.

The plasma is cut into filaments too, carrying the parabolic current profile of
plasma_current_density.i, and moves bodily with its magnetic axis, so
Psi_k(t) = I_p(t) M_k(z_axis(t)). Psi is sampled finely enough in time to be
taken as piecewise linear, which makes its derivative piecewise constant, and
the system is then integrated exactly through the modes of the generalised
eigenproblem R v = lambda L v. It can also be integrated with implicit Euler at
a given time step, which reproduces the time discretisation of the finite
element model and so separates its spatial error from its temporal one.

Run with --help for the options; the defaults match tokamak_wall.i and
plasma_current_density.i.
"""

import argparse
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.special import ellipe, ellipkm1

MU0 = 4e-7 * np.pi


# ---------------------------------------------------------------------------
# Circular loop formulae
# ---------------------------------------------------------------------------
def _elliptic(R1, Z1, R2, Z2):
    """k^2 for two coaxial loops, and K(k^2) and E(k^2). K is evaluated from the
    complementary parameter, which keeps its accuracy for neighbouring loops."""
    denom = (R1 + R2) ** 2 + (Z1 - Z2) ** 2
    m = 4 * R1 * R2 / denom
    m1 = ((R1 - R2) ** 2 + (Z1 - Z2) ** 2) / denom  # 1 - m
    return m, ellipkm1(m1), ellipe(m)


def mutual_inductance(R1, Z1, R2, Z2):
    """Mutual inductance of two coaxial circular filaments (Maxwell)."""
    m, K, E = _elliptic(R1, Z1, R2, Z2)
    k = np.sqrt(m)
    return MU0 * np.sqrt(R1 * R2) * ((2 / k - k) * K - 2 / k * E)


def loop_field_z(Rl, Zl, R, Z):
    """B_z per unit current at (R, Z) from a circular filament of radius Rl at
    height Zl."""
    dz = Z - Zl
    m, K, E = _elliptic(Rl, Zl, R, Z)
    near = (Rl - R) ** 2 + dz ** 2
    return MU0 / (2 * np.pi * np.sqrt((R + Rl) ** 2 + dz ** 2)) * (
        K + (Rl ** 2 - R ** 2 - dz ** 2) / near * E)


def rectangle_gmd(w, h):
    """Geometric mean distance of a w by h rectangle from itself."""
    return np.exp(0.5 * np.log(w ** 2 + h ** 2)
                  - w ** 2 / (12 * h ** 2) * np.log(1 + h ** 2 / w ** 2)
                  - h ** 2 / (12 * w ** 2) * np.log(1 + w ** 2 / h ** 2)
                  + 2 * w / (3 * h) * np.arctan(h / w)
                  + 2 * h / (3 * w) * np.arctan(w / h)
                  - 25.0 / 12.0)


def annulus_cells(R0, r_in, r_out, n_r, n_theta):
    """Cells of the annulus r_in < r < r_out about (R0, 0). Returns the centre of
    each, its minor radius, its area, its width across the annulus and its length
    along it, and the Gauss points and weights of a 3 x 3 rule over it."""
    edges = np.linspace(r_in, r_out, n_r + 1)
    r = 0.5 * (edges[1:] + edges[:-1])
    dtheta = 2 * np.pi / n_theta
    theta = (np.arange(n_theta) + 0.5) * dtheta
    r, theta = (x.ravel() for x in np.meshgrid(r, theta, indexing="ij"))
    w = (r_out - r_in) / n_r
    # Gauss-Legendre points in (r, theta) about each centre, weighted by area.
    xg, wg = np.polynomial.legendre.leggauss(3)
    gr = r[:, None, None] + 0.5 * w * xg[None, :, None]
    gt = theta[:, None, None] + 0.5 * dtheta * xg[None, None, :]
    gw = (wg[:, None] * wg[None, :])[None] * gr
    gw = (gw / gw.sum(axis=(1, 2), keepdims=True)).reshape(r.size, -1)
    gR = (R0 + gr * np.cos(gt)).reshape(r.size, -1)
    gZ = (gr * np.sin(gt)).reshape(r.size, -1)
    return (R0 + r * np.cos(theta), r * np.sin(theta), r, r * dtheta * w,
            np.full(r.size, w), r * dtheta, (gR, gZ, gw))


# ---------------------------------------------------------------------------
# The plasma
# ---------------------------------------------------------------------------
class Plasma:
    """The plasma current of plasma_current_density.i: a parabolic profile over a
    circle of minor radius a, whose magnetic axis is at R0 and rises linearly by
    displacement between displacement_start and displacement_end, carrying a
    current that is piecewise linear in time."""

    def __init__(self, R0, a, waveform, displacement=0.0, displacement_start=0.0,
                 displacement_end=1.0, cells=(12, 48)):
        self.R, self.Z, rho, area, *_ = annulus_cells(R0, 0.0, a, *cells)
        weight = (1 - rho ** 2 / a ** 2) * area
        self.weight = weight / weight.sum()
        self.waveform = tuple(map(np.asarray, waveform))
        self.displacement = (np.array([0.0, displacement_start, displacement_end]),
                             np.array([0.0, 0.0, displacement]))

    def current(self, t):
        return np.interp(t, *self.waveform)

    def axis_height(self, t):
        return np.interp(t, *self.displacement)

    @property
    def moves(self):
        return self.displacement[1][-1] != 0.0

    def knots(self):
        """Times between which the current and the axis height are both linear."""
        return np.union1d(self.waveform[0], self.displacement[0])

    def flux(self, R, Z, t):
        """Poloidal flux of the plasma through the loops (R, Z), at each time t."""
        heights = np.unique(self.axis_height(np.atleast_1d(t)))
        per_amp = {h: mutual_inductance(R[:, None], Z[:, None], self.R[None, :], self.Z[None, :] + h)
                   @ self.weight for h in heights}
        return np.array([self.current(s) * per_amp[self.axis_height(s)] for s in np.atleast_1d(t)])

    def field_z(self, R, Z, t):
        """B_z from the plasma at (R, Z), at each time t."""
        return np.array([self.current(s) * (loop_field_z(self.R, self.Z + self.axis_height(s), R, Z)
                                            @ self.weight) for s in np.atleast_1d(t)])


# ---------------------------------------------------------------------------
# The wall
# ---------------------------------------------------------------------------
class FilamentModel:
    def __init__(self, R0, b, d, conductivity, cells=(4, None)):
        n_r, n_theta = cells
        if n_theta is None:
            # Square cells: as long around the wall as they are across it.
            n_theta = int(round(2 * np.pi * (b + 0.5 * d) / (d / n_r)))
        self.R, self.Z, _, area, w, h, (gR, gZ, gw) = annulus_cells(R0, b, b + d, n_r, n_theta)
        self.resistance = 2 * np.pi * self.R / (conductivity * area)

        # Mutual inductances between cell centres, except for cells within a few cell
        # sizes of each other, where treating each as a filament at its centre is not
        # accurate enough to keep the matrix positive definite: those are averaged
        # over Gauss points in both cells, which accounts for their finite size.
        n = self.R.size
        L = np.empty((n, n))
        for i in range(n):
            with np.errstate(divide="ignore", invalid="ignore"):  # the diagonal, set below
                L[i] = mutual_inductance(self.R[i], self.Z[i], self.R, self.Z)
            near = np.hypot(self.R - self.R[i], self.Z - self.Z[i]) < 4 * max(w[i], h[i])
            near[i] = False
            idx = np.nonzero(near)[0]
            if idx.size:
                m = mutual_inductance(gR[i][None, :, None], gZ[i][None, :, None],
                                      gR[idx][:, None, :], gZ[idx][:, None, :])
                L[i, idx] = np.einsum("a,jab,jb->j", gw[i], m, gw[idx])
        L[np.diag_indices(n)] = MU0 * self.R * (np.log(8 * self.R / rectangle_gmd(w, h)) - 2)
        self.inductance = 0.5 * (L + L.T)

        # Modes: R v = lambda L v, normalised so that V^T L V = 1 and V^T R V = diag(lambda).
        self.rates, self.modes = eigh(np.diag(self.resistance), self.inductance)

    @property
    def time_constants(self):
        return 1.0 / self.rates

    def solve_exact(self, plasma, times, max_step=1e-4):
        """Wall filament currents at the given times. The plasma flux is taken as
        linear between the plasma's own knots, the output times and, while the
        plasma moves, points no more than max_step apart."""
        knots = np.union1d(plasma.knots(), times)
        if plasma.moves:
            t0, t1 = plasma.displacement[0][1:]
            knots = np.union1d(knots, np.linspace(t0, t1, int(np.ceil((t1 - t0) / max_step)) + 1))
        knots = knots[(knots >= 0) & (knots <= times[-1])]
        psi = plasma.flux(self.R, self.Z, knots)
        V, lam = self.modes, self.rates
        # Modal amplitudes c, with I = V c, obey dc/dt + lambda c = -V^T dPsi/dt.
        c = np.zeros(lam.size)
        out = np.zeros((len(times), lam.size))
        out_index = {t: i for i, t in enumerate(times)}
        for k in range(1, len(knots)):
            dt = knots[k] - knots[k - 1]
            drive = -V.T @ (psi[k] - psi[k - 1]) / dt
            decay = np.exp(-lam * dt)
            c = c * decay + drive / lam * (1 - decay)
            if knots[k] in out_index:
                out[out_index[knots[k]]] = c
        return out @ V.T

    def solve_implicit_euler(self, plasma, times):
        """As solve_exact, but integrated by implicit Euler with steps between the
        given times, as the finite element model is."""
        L, R = self.inductance, np.diag(self.resistance)
        psi = plasma.flux(self.R, self.Z, times)
        out = np.zeros((len(times), L.shape[0]))
        inverses = {}
        for n in range(1, len(times)):
            dt = times[n] - times[n - 1]
            key = round(dt, 12)
            if key not in inverses:
                inverses[key] = np.linalg.inv(L / dt + R)
            out[n] = inverses[key] @ (L @ out[n - 1] - (psi[n] - psi[n - 1])) / dt
        return out

    def field_z(self, wall_currents, R, Z):
        """B_z at (R, Z) from the wall currents."""
        return wall_currents @ loop_field_z(self.R, self.Z, R, Z)

    def power(self, wall_currents):
        return (wall_currents ** 2) @ self.resistance


# Points at which B_z is reported, matching the point samplers in tokamak_wall.i:
# the centre of the machine, and the outboard midplane between plasma and wall.
PROBES = {"bz_centre": (0.0, 0.0), "bz_outboard": (1.35, 0.0)}


def build(args):
    plasma = Plasma(args.major_radius, args.plasma_radius,
                    ([0.0, args.ramp_end, args.flat_top_end, args.quench_end],
                     [0.0, args.peak_current, args.peak_current, 0.0]),
                    args.vertical_displacement, args.displacement_start, args.displacement_end)
    wall = FilamentModel(args.major_radius, args.wall_radius, args.wall_thickness,
                         args.conductivity, (args.wall_radial_cells, args.wall_poloidal_cells))
    return plasma, wall


def reference(args, times, implicit_euler=False, models=None):
    """Time histories of the quantities that tokamak_wall.i reports, at the given
    times, which must start from zero."""
    plasma, wall = models or build(args)
    I = (wall.solve_implicit_euler if implicit_euler else wall.solve_exact)(plasma, times)
    cols = {"time": times, "plasma_current": plasma.current(times), "wall_current": I.sum(axis=1),
            "wall_current_moment": I @ wall.Z, "wall_power": wall.power(I)}
    for name, (R, Z) in PROBES.items():
        cols[name] = wall.field_z(I, R, Z) + plasma.field_z(R, Z, times)
    return cols


def make_parser():
    """The command line options, whose defaults match tokamak_wall.i and
    plasma_current_density.i."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--major-radius", type=float, default=1.0)
    p.add_argument("--plasma-radius", type=float, default=0.3)
    p.add_argument("--wall-radius", type=float, default=0.4)
    p.add_argument("--wall-thickness", type=float, default=0.03)
    p.add_argument("--conductivity", type=float, default=1.0 / 7.4e-7,
                   help="of the wall (S/m); 316L stainless steel at room temperature")
    p.add_argument("--peak-current", type=float, default=1e6, help="A")
    p.add_argument("--ramp-end", type=float, default=0.010, help="s")
    p.add_argument("--flat-top-end", type=float, default=0.030, help="s")
    p.add_argument("--quench-end", type=float, default=0.032, help="s")
    p.add_argument("--vertical-displacement", type=float, default=0.0, help="m")
    p.add_argument("--displacement-start", type=float, default=0.020, help="s")
    p.add_argument("--displacement-end", type=float, default=0.030, help="s")
    p.add_argument("--end-time", type=float, default=0.060, help="s")
    p.add_argument("--dt", type=float, default=2.5e-4,
                   help="output interval and implicit Euler step (s)")
    p.add_argument("--wall-radial-cells", type=int, default=4,
                   help="cells through the wall thickness")
    p.add_argument("--wall-poloidal-cells", type=int,
                   help="cells around the wall; by default, enough to make them square")
    p.add_argument("--implicit-euler", action="store_true",
                   help="integrate with implicit Euler at --dt instead of exactly")
    p.add_argument("--csv", help="write the time histories here")
    p.add_argument("--convergence", action="store_true",
                   help="print the wall time constant and peak wall current as the cells are refined")
    return p


def main():
    args = make_parser().parse_args()
    times = np.linspace(0, args.end_time, int(round(args.end_time / args.dt)) + 1)

    if args.convergence:
        print(f"{'cells':>12} {'tau_wall (ms)':>14} {'peak I_wall (kA)':>17}")
        for nr in (1, 2, 4, 8):
            args.wall_radial_cells = nr
            plasma, wall = build(args)
            cols = reference(args, times, models=(plasma, wall))
            print(f"{nr:>4} x {wall.R.size // nr:<5} {1e3 * wall.time_constants.max():14.5f} "
                  f"{1e-3 * np.abs(cols['wall_current']).max():17.4f}")
        return

    plasma, wall = build(args)
    cols = reference(args, times, args.implicit_euler, models=(plasma, wall))
    tau = np.sort(wall.time_constants)[::-1]
    print(f"wall time constants (ms): {', '.join(f'{1e3 * t:.4f}' for t in tau[:4])}, ...")
    k = np.argmax(np.abs(cols["wall_current"]))
    print(f"peak wall current {cols['wall_current'][k]:.6g} A at t = {times[k]:.6g} s")
    if args.csv:
        np.savetxt(args.csv, np.column_stack(list(cols.values())), delimiter=",",
                   header=",".join(cols), comments="", fmt="%.10g")
        print(f"wrote {args.csv}")


if __name__ == "__main__":
    sys.exit(main())
