"""
Runs a vertical displacement event (VDE) in FREEGSnke on the machine of
torus_machine.py and exports its source field for torus_vessel.i.

The equilibrium is an elongated plasma limited at the top and bottom of the
limiter, found by an inverse solve. Elongation makes it vertically unstable; the
vessel slows the instability to the timescale of its eddy currents. The poloidal
field coils are held at constant voltage V = R I, so they would keep their
currents but for induction, and a constant voltage on the antisymmetric PR coil
pair builds up a radial field that starts the plasma moving.

Writes to the output directory:

- the source field excluding the vessel, via freegsnke_source.py;
- freegsnke_vde.csv, FREEGSnke's history of the plasma and of the vessel's total
  toroidal current and vertical force, for comparison with the MFEM solution.

Usage: python run_vde.py OUTPUT_DIR [--dt DT] [--steps N] [--kick-voltage V]
       [--max-axis-shift Z] [--n-poloidal N] [--n-radial N]
"""

import argparse
import csv
import os
import sys

import numpy as np

# The FREEGSnke scripts live in vde/, beside gridded_fields.py in the example directory.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]
import torus_machine as tm
from freegs4e.gradshafranov import GreensBr
from freegsnke import GSstaticsolver, equilibrium_update, nonlinear_solve
from freegsnke.inverse import Inverse_optimizer
from freegsnke.jtor_update import ConstrainBetapIp

from freegsnke_source import SourceFieldRecorder

# Plasma current [A], poloidal beta, and R B_phi of the vacuum toroidal field
# [T m], giving 2 T on the magnetic axis.
IP = 1.5e5
BETAP = 0.3
FVAC = 10.0
# Plasma resistivity [Ohm m], low enough that the plasma current decays slowly
# compared with the vertical motion.
PLASMA_RESISTIVITY = 1e-8
# Half-width of the target plasma boundary [m]; its half-height is the limiter
# radius, so the plasma touches the limiter at top and bottom with elongation 1.38.
HALF_WIDTH = 0.32


def equilibrium():
    tokamak = tm.build()
    # The grid extends vertically to include the X-points of the shaping field,
    # at Z = +-1.03 m. FREEGSnke's time-dependent solver needs the core mask they
    # define even for a limited plasma.
    eq = equilibrium_update.Equilibrium(
        tokamak=tokamak, Rmin=4.2, Rmax=5.8, Zmin=-1.3, Zmax=1.3, nx=65, ny=129
    )
    profiles = ConstrainBetapIp(eq=eq, betap=BETAP, Ip=IP, fvac=FVAC, alpha_m=1.8, alpha_n=1.2)
    solver = GSstaticsolver.NKGSsolver(eq)

    eq.tokamak.set_coil_current("PR", 0.0)
    eq.tokamak["PR"].control = False
    theta = np.linspace(0, 2 * np.pi, 13)[:-1]
    isoflux_set = np.array(
        [[list(tm.R0 + HALF_WIDTH * np.cos(theta)), list(tm.LIMITER_RADIUS * np.sin(theta))]]
    )
    solver.solve(
        eq=eq,
        profiles=profiles,
        constrain=Inverse_optimizer(isoflux_set=isoflux_set),
        target_relative_tolerance=1e-6,
        target_relative_psit_update=1e-3,
    )
    solver.solve(eq=eq, profiles=profiles, constrain=None, target_relative_tolerance=1e-9)
    return eq, profiles, solver


class VesselDiagnostics:
    """FREEGSnke's total toroidal current in, and vertical force on, the vessel."""

    def __init__(self, tokamak):
        self.labels = set(tm.vessel_labels())
        R, Z, weight, index = [], [], [], []
        for k, label in enumerate(tokamak.coils_list):
            if label not in self.labels:
                continue
            coil = tokamak.coils_dict[label]
            n_filaments = len(coil["coords"][0])
            R.extend(coil["coords"][0])
            Z.extend(coil["coords"][1])
            # Filament current per unit structure current. A refined passive
            # structure stores one polarity and multiplier for all its filaments.
            weight.extend(
                np.broadcast_to(
                    np.asarray(coil["polarity"]) * np.asarray(coil["multiplier"]), (n_filaments,)
                )
            )
            index.extend([k] * n_filaments)
        self.R, self.Z = np.array(R), np.array(Z)
        self.weight, self.index = np.array(weight), np.array(index)

    def __call__(self, eq, profiles):
        currents = eq.tokamak.getCurrentsVec()
        filament_currents = currents[self.index] * self.weight

        # Radial field at the filaments from everything but the vessel, whose
        # forces on itself sum to zero.
        b_r = np.zeros_like(self.R)
        for label, coil in eq.tokamak.coils:
            if label not in self.labels:
                b_r += coil.Br(self.R, self.Z)
        mask = profiles.jtor != 0
        for Rc, Zc, current in zip(eq.R[mask], eq.Z[mask], profiles.jtor[mask] * eq.dR * eq.dZ):
            b_r += current * GreensBr(Rc, Zc, self.R, self.Z)

        # f = J x B with J along phi_hat, so f_z = -J_phi B_R.
        vertical_force = -np.sum(2 * np.pi * self.R * filament_currents * b_r)
        return filament_currents.sum(), vertical_force


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir")
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--kick-voltage", type=float, default=0.5)
    # Vertical shift of the magnetic axis [m] at which to stop. The plasma touches
    # the limiter at top and bottom, so as it moves it is scraped off; by a shift of
    # about 0.14 m the next step's plasma region vanishes and FREEGSnke stalls.
    parser.add_argument("--max-axis-shift", type=float, default=0.12)
    parser.add_argument("--n-poloidal", type=int, default=tm.N_POLOIDAL,
                        help="Number of poloidal sectors the vessel is split into.")
    parser.add_argument("--n-radial", type=int, default=tm.N_RADIAL,
                        help="Number of radial layers the vessel is split into.")
    args = parser.parse_args()
    tm.N_POLOIDAL = args.n_poloidal
    tm.N_RADIAL = args.n_radial

    eq, profiles, solver = equilibrium()
    stepping = nonlinear_solve.nl_solver(
        eq=eq,
        profiles=profiles,
        GSStaticSolver=solver,
        full_timestep=args.dt,
        max_internal_timestep=args.dt,
        plasma_resistivity=PLASMA_RESISTIVITY,
        verbose=True,
    )
    print("Growth rates [1/s]:", stepping.linearised_sol.growth_rates)

    n_active = stepping.evol_metal_curr.n_active_coils
    voltages = (stepping.vessel_currents_vec * stepping.evol_metal_curr.coil_resist)[:n_active]
    voltages[eq.tokamak.coils_list.index("PR")] = args.kick_voltage

    # The export grid covers the vessel's poloidal footprint with a margin.
    recorder = SourceFieldRecorder(
        np.linspace(3.95, 6.05, 43), np.linspace(-1.05, 1.05, 43), tm.vessel_labels()
    )
    diagnostics = VesselDiagnostics(eq.tokamak)

    stepping.initialize_from_ICs(eq, profiles)
    os.makedirs(args.output_dir, exist_ok=True)
    with open(os.path.join(args.output_dir, "freegsnke_vde.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["time", "axis_r", "axis_z", "plasma_current", "pr_current", "vessel_current",
             "vessel_vertical_force"]
        )
        t = 0.0
        for step in range(args.steps + 1):
            if step > 0:
                # As the plasma is scraped off against the limiter it shrinks until
                # the solver fails; keep the history recorded up to that point.
                try:
                    stepping.nlstepper(active_voltage_vec=voltages, linear_only=False)
                except Exception as error:
                    print(f"Stopping at t = {t:.3f} s: the step failed with: {error!r}")
                    break
                t += stepping.dt_step
            eq1, profiles1 = stepping.eq1, stepping.profiles1
            recorder.record(t, eq1, profiles1)
            vessel_current, vessel_force = diagnostics(eq1, profiles1)
            axis = eq1.opt[0]
            writer.writerow(
                [t, axis[0], axis[1], profiles1.Ip, eq1.tokamak["PR"].current, vessel_current,
                 vessel_force]
            )
            f.flush()
            print(f"t = {t:.3f} s, axis Z = {axis[1]:.4f} m, Ip = {profiles1.Ip:.4g} A, "
                  f"vessel I = {vessel_current:.4g} A, F_z = {vessel_force:.4g} N", flush=True)
            if abs(axis[1]) > args.max_axis_shift:
                print(f"Stopping at t = {t:.3f} s: the magnetic axis has moved more than "
                      f"{args.max_axis_shift} m.")
                break

    recorder.write(args.output_dir, file_base="vde_source")


if __name__ == "__main__":
    main()
