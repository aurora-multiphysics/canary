#!/usr/bin/env python3
"""Semi-analytic resistance of a circular contact spot between two coaxial rods.

The reference solution for the electrical contact example, documented in
doc/content/problems/electrical_contact.md. It is independent of the finite element
model: the potential in each rod is written as a Fourier-Bessel series, and the
current distribution through the contact spot is found by minimising the power
dissipated, which is what the true current distribution does (Thomson's theorem).

Two rods of radius b, with resistivities rho_1 and rho_2, are pressed together end
to end. They touch over a spot of radius a on the axis, and are insulated from
each other everywhere else on the contact plane. The spot is covered by a film of
contact conductance C (S/m^2), across which the potential jumps by q / C, with q
the current density through it. Each rod runs a length l from the contact plane to
an electrode far enough away that the current there is uniform. The rods are
otherwise alike: the analysis needs only the same radius and length for both.

A current density q(r) through the spot sets up, in a rod of resistivity rho,

  phi(r, z) = rho (C_0 - I z / (pi b^2) + sum_n c_n J0(k_n r) exp(-k_n z)),

with J1(k_n b) = 0, so that no current crosses the side of the rod, and

  c_n k_n = Q_n / N_n,   Q_n = int_0^a q J0(k_n r) r dr,   N_n = b^2 J0(k_n b)^2 / 2.

The power dissipated, over and above that of a uniform current, is

  rho * 2 pi sum_n Q_n^2 / (k_n N_n),

and the same current density flows into both rods, so with the film it is

  R I^2 = (rho_1 + rho_2) (I^2 l / (pi b^2) + 2 pi sum_n Q_n^2 / (k_n N_n))
          + (1 / C) int_0^a q^2 2 pi r dr.

q is found by minimising R over a finite family of distributions: any one of them
gives an upper bound on R, which converges from above as the family grows. A clean
contact (C infinite) is an equipotential spot, with an inverse square root
singularity in q at its edge. q is expanded there in (1 - r^2/a^2)^(m - 1/2), whose
Hankel transforms are known in closed form (Sonine's integral), and a dozen terms
are converged to rounding error. Those distributions would dissipate infinite power
in a film, so with one q is taken constant on each of a set of annuli instead,
graded towards the edge of the spot, where the current crowds into a layer of
width about 1 / (C (rho_1 + rho_2)). Their Hankel transforms are closed form too,
and the error falls as the square of the width of the annuli, so the result is
extrapolated from two sets of them. The Bessel series is summed to a few hundred
thousand terms; what is left of it falls as one over the number of terms, and is
extrapolated in the same way.

Two simple distributions bound R. The equipotential spot is the one of least
constriction resistance and uniform current density the one of least film
resistance, so R is at least their sum, which is Holm's formula,

  R >= R_bulk + (rho_1 + rho_2) psi_eq(a/b) / (4 a) + 1 / (C pi a^2),

and at most the resistance of a uniform current density, which is an upper bound
like any other trial distribution,

  R <= R_bulk + (rho_1 + rho_2) psi_uniform(a/b) * 8 / (3 pi^2 a) + 1 / (C pi a^2).

Run with --help for the options. Lengths are in metres.
"""

import argparse
import math

import numpy as np
from scipy import special


def bessel_j1_zeros(count):
    """The first `count` positive zeros of J1: exact for the first thousand, then
    McMahon's expansion, which is good to rounding error from there on."""
    exact = min(count, 1000)
    zeros = special.jn_zeros(1, exact)
    if count > exact:
        beta = (np.arange(exact + 1, count + 1) + 0.25) * math.pi
        m = 4.0  # 4 nu^2 with nu = 1
        zeros = np.concatenate([
            zeros,
            beta - (m - 1) / (8 * beta) - 4 * (m - 1) * (7 * m - 31) / (3 * (8 * beta) ** 3)
            - 32 * (m - 1) * (83 * m**2 - 982 * m + 3779) / (15 * (8 * beta) ** 5)])
    return zeros


def constriction_matrix(transforms, a, b, terms, chunk=20000):
    """2 pi sum_n Q_n^i Q_n^j / (k_n N_n) over a family of current densities, each
    given by a function of k returning its Hankel transforms Q_n: the power that
    pairs of them dissipate in a semi-infinite rod of unit resistivity, over and
    above a uniform current. Summed in chunks to keep the memory down."""
    lam = bessel_j1_zeros(terms)
    matrix = 0
    for start in range(0, terms, chunk):
        k = lam[start:start + chunk] / b
        norm = 0.5 * b**2 * special.j0(k * b) ** 2
        q = transforms(k)
        matrix = matrix + 2 * math.pi * (q / (k * norm)) @ q.T
    return matrix


def least_power(matrix, current):
    """The least of alpha^T M alpha subject to current^T alpha = 1: the resistance
    of the best combination of the family carrying unit current, and the
    combination, alpha."""
    x = np.linalg.solve(matrix, current)
    resistance = 1.0 / (current @ x)
    return resistance, resistance * x


def clean_spot(a, b, rho, terms, basis=12):
    """Constriction resistance of an equipotential spot, in the expansion
    q = sum_m alpha_m (1 - r^2/a^2)^(m - 1/2), and a function giving q(r) for unit
    current."""
    mus = np.arange(basis) - 0.5

    def transforms(k):  # Sonine's first integral
        ka = k[None, :] * a
        mu = mus[:, None]
        return a**2 * 2**mu * special.gamma(mu + 1) * special.jv(mu + 1, ka) / ka ** (mu + 1)

    resistance, alpha = least_power(rho * constriction_matrix(transforms, a, b, terms),
                                    math.pi * a**2 / (mus + 1))
    return resistance, lambda r: alpha @ (1 - (np.asarray(r)[None, :] / a) ** 2) ** mus[:, None]


def filmed_spot(a, b, rho, conductance, terms, annuli, grading=3):
    """Resistance of the spot and its film, with q constant on each of `annuli`
    annuli, graded towards the edge of the spot, and the radii bounding the annuli
    and q on each for unit current."""
    r = a * (1 - (1 - np.linspace(0, 1, annuli + 1)) ** grading)
    area = math.pi * np.diff(r**2)

    def transforms(k):  # int r J0(k r) dr = r J1(k r) / k
        f = r[:, None] * special.j1(k[None, :] * r[:, None]) / k[None, :]
        return np.diff(f, axis=0)

    matrix = rho * constriction_matrix(transforms, a, b, terms) + np.diag(area) / conductance
    resistance, q = least_power(matrix, area)
    return resistance, r, q


def uniform_spot(a, b, rho, conductance, terms):
    """Resistance of the spot and its film for a uniform current density."""
    area = math.pi * a**2
    q = constriction_matrix(lambda k: (a * special.j1(k * a) / k)[None, :], a, b, terms)
    return (rho * q[0, 0] + area / conductance) / area**2


def extrapolate(coarse, fine, ratio):
    """Richardson extrapolation of an error falling as 1/ratio per refinement."""
    return fine + (fine - coarse) / (ratio - 1)


def spot_resistance(a, b, rho, conductance=math.inf, terms=200000, annuli=200):
    """R - R_bulk: the resistance of the spot and its film over and above that of
    the rods carrying a uniform current, extrapolated in the number of Bessel terms
    and, with a film, in the width of the annuli."""
    if math.isinf(conductance):
        return extrapolate(clean_spot(a, b, rho, terms)[0], clean_spot(a, b, rho, 4 * terms)[0], 4)
    coarse = filmed_spot(a, b, rho, conductance, terms, annuli)[0]
    fine = filmed_spot(a, b, rho, conductance, terms, 2 * annuli)[0]
    return extrapolate(coarse, fine, 4)


def spot_current_density(a, b, rho, conductance=math.inf, terms=200000, annuli=400):
    """The current density through the spot for unit current, at radii r, as a
    pair of arrays (r, q). For a filmed spot these are the midpoints of the annuli."""
    if math.isinf(conductance):
        r = a * np.linspace(0, 1, 401)[:-1]
        return r, clean_spot(a, b, rho, terms)[1](r)
    _, edges, q = filmed_spot(a, b, rho, conductance, terms, annuli)
    return 0.5 * (edges[1:] + edges[:-1]), q


def contact_resistance(spot_radius, rod_radius, rod_length, sigma_1, sigma_2,
                       conductance=math.inf, terms=200000):
    """The resistance between the electrodes, with the bounds of Holm's formula
    (below) and of a uniform current density over the spot (above). rod_length is
    the length of each rod, from its electrode to its face against the film."""
    a, b = spot_radius, rod_radius
    rho = 1 / sigma_1 + 1 / sigma_2
    bulk = rho * rod_length / (math.pi * b**2)
    film = 0.0 if math.isinf(conductance) else 1 / (conductance * math.pi * min(a, b) ** 2)
    if a >= b:  # one dimensional
        return dict(resistance=bulk + film, lower=bulk + film, upper=bulk + film,
                    bulk=bulk, film=film, constriction=0.0)
    clean = spot_resistance(a, b, rho, math.inf, terms)
    exact = spot_resistance(a, b, rho, conductance, terms)
    upper = math.inf if math.isinf(conductance) else uniform_spot(a, b, rho, conductance, terms)
    # The bounds are the sum of the resistances of their distributions: the
    # equipotential spot and the film, below, and the uniform current, above.
    return dict(resistance=bulk + exact, lower=bulk + clean + film, upper=bulk + upper,
                bulk=bulk, film=film, constriction=clean, psi=clean / (rho / (4 * a)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spot-radius", type=float, default=1e-3, help="a (m)")
    p.add_argument("--rod-radius", type=float, default=5e-3, help="b (m)")
    p.add_argument("--rod-length", type=float, default=20e-3 - 0.5e-6,
                   help="length of each rod, from its electrode to its face against "
                        "the film (m)")
    p.add_argument("--sigma-1", type=float, default=5.8e7, help="conductivity of rod 1 (S/m)")
    p.add_argument("--sigma-2", type=float, default=3.77e7, help="conductivity of rod 2 (S/m)")
    p.add_argument("--conductance", type=float, default=math.inf,
                   help="contact conductance of the film over the spot (S/m^2); "
                        "inf for a clean contact")
    p.add_argument("--terms", type=int, default=200000,
                   help="terms in the Bessel series before extrapolation")
    args = p.parse_args()

    res = contact_resistance(args.spot_radius, args.rod_radius, args.rod_length,
                             args.sigma_1, args.sigma_2, args.conductance, args.terms)
    for k, v in res.items():
        print(f"{k:13s} {v:.8e}")


if __name__ == "__main__":
    main()
