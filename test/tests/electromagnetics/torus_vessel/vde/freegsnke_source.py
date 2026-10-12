"""
Records the magnetic source field of a FREEGSnke evolution for MOOSE-MFEM
eddy-current and load calculations on 3D in-vessel components.

The source is everything FREEGSnke knows about except the passive structures
that the MFEM model resolves itself: the active coils, any passive structures
left out of the MFEM mesh, and the plasma. Leaving the resolved structures out
keeps their field from being counted twice, once by FREEGSnke's circuit model
and once by the MFEM eddy currents.

Conventions follow FreeGS: (R, phi, Z) is right handed, psi is the poloidal flux
per radian, psi = R A_phi, B_R = -(1/R) dpsi/dZ and B_Z = (1/R) dpsi/dR. In a 3D
Cartesian model with the symmetry axis along z, phi_hat = (-y, x, 0) / R, so the
source electric field driving the eddy currents is

    E_s = -(1/R) dpsi_s/dt phi_hat = dpsi_s/dt (y, -x, 0) / R^2.

Fields are written by gridded_fields.py and mapped onto the 3D mesh by
SolutionUserObject and Axisymmetric2D3DSolutionFunction.
"""

import numpy as np
from freegs4e.gradshafranov import Greens, GreensBr, GreensBz
from scipy.interpolate import RectBivariateSpline

from gridded_fields import write_fields


class SourceFieldRecorder:
    """
    Accumulates the source field on a fixed (R, Z) grid at each recorded time.

    Parameters
    ----------
    R, Z : array-like
        Monotonically increasing grid coordinates [m]. The grid must cover the
        poloidal footprint of every component in the MFEM model; it need not
        coincide with, or lie inside, the FREEGSnke computational grid.
    modelled_passives : iterable of str
        Labels of the FREEGSnke passive structures that the MFEM model resolves.
        Their currents are excluded from the source.
    """

    def __init__(self, R, Z, modelled_passives):
        self.R = np.asarray(R, dtype=float)
        self.Z = np.asarray(Z, dtype=float)
        if np.any(np.diff(self.R) <= 0) or np.any(np.diff(self.Z) <= 0):
            raise ValueError("R and Z must be strictly increasing.")
        if self.R[0] <= 0:
            raise ValueError("R must be positive: the source field is singular on the axis.")
        self._RR, self._ZZ = np.meshgrid(self.R, self.Z, indexing="ij")
        self.modelled_passives = set(modelled_passives)
        self.times = []
        self._psi = []
        self._b_r = []
        self._b_z = []
        self._b_phi = []

    def record(self, t, eq, profiles):
        """
        Evaluate the source field of an equilibrium at time t.

        Call once for the initial equilibrium and then after every nlstepper
        step, passing the stepper's eq1 and profiles1.
        """
        if self.times and t <= self.times[-1]:
            raise ValueError("Times must be recorded in increasing order.")

        labels = {label for label, _ in eq.tokamak.coils}
        unknown = self.modelled_passives - labels
        if unknown:
            raise ValueError(f"Unknown passive structure labels: {sorted(unknown)}")

        RR, ZZ = self._RR, self._ZZ
        psi = np.zeros_like(RR)
        b_r = np.zeros_like(RR)
        b_z = np.zeros_like(RR)

        for label, coil in eq.tokamak.coils:
            if label in self.modelled_passives:
                continue
            psi += coil.psi(RR, ZZ)
            b_r += coil.Br(RR, ZZ)
            b_z += coil.Bz(RR, ZZ)

        # Inside the FREEGSnke grid the plasma flux is FREEGSnke's own solution.
        # Outside it the plasma is summed as one filament per grid cell carrying
        # current; those cells lie inside the limiter, so the Green's functions
        # are never evaluated at a filament.
        inside = (
            (RR >= eq.R_1D[0])
            & (RR <= eq.R_1D[-1])
            & (ZZ >= eq.Z_1D[0])
            & (ZZ <= eq.Z_1D[-1])
        )
        spline = RectBivariateSpline(eq.R_1D, eq.Z_1D, eq.plasma_psi)
        R_in, Z_in = RR[inside], ZZ[inside]
        psi[inside] += spline.ev(R_in, Z_in)
        b_r[inside] -= spline.ev(R_in, Z_in, dy=1) / R_in
        b_z[inside] += spline.ev(R_in, Z_in, dx=1) / R_in

        outside = ~inside
        if np.any(outside):
            R_out, Z_out = RR[outside], ZZ[outside]
            jtor = profiles.jtor
            mask = jtor != 0.0
            currents = jtor[mask] * eq.dR * eq.dZ
            for Rc, Zc, current in zip(eq.R[mask], eq.Z[mask], currents):
                psi[outside] += current * Greens(Rc, Zc, R_out, Z_out)
                b_r[outside] += current * GreensBr(Rc, Zc, R_out, Z_out)
                b_z[outside] += current * GreensBz(Rc, Zc, R_out, Z_out)

        self.times.append(float(t))
        self._psi.append(psi)
        self._b_r.append(b_r)
        self._b_z.append(b_z)
        # Outside the plasma the poloidal current function takes its vacuum value.
        self._b_phi.append(profiles.fvac() / RR)

    def fields(self):
        """
        Return the recorded fields as arrays of shape (n_times, n_R, n_Z).

        dpsi_dt at each time is the backward difference over the step that
        ends there, which matches an implicit Euler step over the same
        interval. The evolution is taken to start from rest, so it is zero at
        the first time.
        """
        if len(self.times) < 2:
            raise RuntimeError("At least two times are needed to form dpsi/dt.")
        psi = np.array(self._psi)
        dt = np.diff(self.times)[:, None, None]
        dpsi_dt = np.zeros_like(psi)
        dpsi_dt[1:] = np.diff(psi, axis=0) / dt
        return {
            "dpsi_dt": dpsi_dt,
            "b_r": np.array(self._b_r),
            "b_z": np.array(self._b_z),
            "b_phi": np.array(self._b_phi),
        }

    def write(self, directory, file_base="source"):
        """
        Write the recorded fields with gridded_fields.write_fields: one
        PiecewiseMultilinear data file per field, a CSV of the times, and the
        MOOSE input that converts them to <file_base>_out.e.
        """
        write_fields(directory, self.R, self.Z, self.times, self.fields(), file_base)
