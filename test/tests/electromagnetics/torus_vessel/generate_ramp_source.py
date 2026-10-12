"""
Writes the source fields for torus_vessel.i: a uniform vertical field ramped at
a constant rate, B_z = Bdot t, as produced by distant poloidal field coils.

Its flux per radian is psi_s = B_z R^2 / 2, which satisfies the vacuum
Grad-Shafranov equation, so

    dpsi_dt = Bdot R^2 / 2,   b_r = 0,   b_z = Bdot t,   b_phi = 0.

The files are written with gridded_fields.py, the same path used for
FREEGSnke exports, so the test exercises the production data chain. Run from
this directory.
"""

import numpy as np

from gridded_fields import FIELDS, write_fields

# Field ramp rate [T/s].
BDOT = 0.1

# The grid covers the poloidal footprint of the TorusSheath block, R in [4, 6]
# and Z in [-1, 1], with a margin of one grid spacing.
SPACING = 0.05
R = np.linspace(3.95, 6.05, 43)
Z = np.linspace(-1.05, 1.05, 43)

# The fields are linear in time, so two times interpolate them exactly over the
# run in torus_vessel.i, which ends at 10 s.
TIMES = np.array([0.0, 10.0])

RR = np.broadcast_to(R[:, None], (len(R), len(Z)))
fields = {name: np.zeros((len(TIMES), len(R), len(Z))) for name in FIELDS}
for n, t in enumerate(TIMES):
    fields["dpsi_dt"][n] = BDOT * RR**2 / 2
    fields["b_z"][n] = BDOT * t

assert np.allclose(np.diff(R), SPACING) and np.allclose(np.diff(Z), SPACING)
write_fields(".", R, Z, TIMES, fields, file_base="ramp_source")
