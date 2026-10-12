"""
A FREEGSnke machine built around the torus of MOOSE's embedded_concentric_torus.e,
so that FREEGSnke and the MFEM model in torus_vessel.i describe the same vessel.

- The vessel is the TorusSheath block: a shell of 316L stainless steel centred on
  (R, Z) = (5, 0) m between minor radii 0.5 m and 1 m. It is split into
  N_POLOIDAL poloidal sectors and N_RADIAL radial layers, each a passive
  structure, so that FREEGSnke can redistribute current around and through it.
- The plasma lives in the non-conducting TorusCore block, of minor radius 0.5 m,
  inside a circular limiter of radius 0.45 m.
- Up-down symmetric poloidal field coil pairs supply the vertical field and the
  shaping field that elongates the plasma; outside the shell they do not overlap
  the MFEM conductors. The PR pair is wired antisymmetrically and makes a radial
  field, used to displace the plasma.
"""

import numpy as np

R0 = 5.0
VESSEL_INNER = 0.5
VESSEL_OUTER = 1.0
LIMITER_RADIUS = 0.45

# 316L stainless steel at room temperature [Ohm m], the inverse of the
# conductivity in torus_vessel.i.
VESSEL_RESISTIVITY = 1 / 1.35e6
# Copper [Ohm m].
COPPER_RESISTIVITY = 1.55e-8

N_POLOIDAL = 24
N_RADIAL = 2


def vessel_label(i, j):
    """Label of the passive structure in poloidal sector i and radial layer j."""
    return f"vessel_{i:02d}_{j}"


def vessel_labels():
    return [vessel_label(i, j) for i in range(N_POLOIDAL) for j in range(N_RADIAL)]


def passive_coils():
    """The vessel shell as quadrilateral sectors of the annulus."""
    passives = []
    angles = np.linspace(0, 2 * np.pi, N_POLOIDAL + 1)
    radii = np.linspace(VESSEL_INNER, VESSEL_OUTER, N_RADIAL + 1)
    for i in range(N_POLOIDAL):
        for j in range(N_RADIAL):
            corners = [
                (radii[j], angles[i]),
                (radii[j + 1], angles[i]),
                (radii[j + 1], angles[i + 1]),
                (radii[j], angles[i + 1]),
            ]
            passives.append(
                {
                    "R": [R0 + r * np.cos(a) for r, a in corners],
                    "Z": [r * np.sin(a) for r, a in corners],
                    "name": vessel_label(i, j),
                    "resistivity": VESSEL_RESISTIVITY,
                }
            )
    return passives


def _coil_pair(R, Z, polarity_lower):
    """A pair of single-winding 0.2 m square coils at (R, +Z) and (R, -Z)."""
    return {
        side: {
            "R": [R],
            "Z": [sign * Z],
            "dR": 0.2,
            "dZ": 0.2,
            "resistivity": COPPER_RESISTIVITY,
            "polarity": polarity,
            "multiplier": 1,
        }
        for side, sign, polarity in (("upper", 1, 1), ("lower", -1, polarity_lower))
    }


def active_coils():
    return {
        "PF1": _coil_pair(4.0, 1.6, 1),
        "PF2": _coil_pair(6.0, 1.6, 1),
        "PF3": _coil_pair(6.9, 0.7, 1),
        "PF4": _coil_pair(3.2, 0.7, 1),
        "PR": _coil_pair(6.5, 1.3, -1),
    }


def limiter():
    angles = np.linspace(0, 2 * np.pi, 73)[:-1]
    return [
        {"R": R0 + LIMITER_RADIUS * np.cos(a), "Z": LIMITER_RADIUS * np.sin(a)} for a in angles
    ]


def wall():
    angles = np.linspace(0, 2 * np.pi, 73)
    return [{"R": R0 + VESSEL_INNER * np.cos(a), "Z": VESSEL_INNER * np.sin(a)} for a in angles]


def build():
    from freegsnke import build_machine

    return build_machine.tokamak(
        active_coils_data=active_coils(),
        passive_coils_data=passive_coils(),
        limiter_data=limiter(),
        wall_data=wall(),
    )
