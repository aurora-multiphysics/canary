#!/usr/bin/env python3
"""Compare tokamak_wall.i with the filament model, and draw the figure for the docs.

Reads the postprocessor history written by tokamak_wall.i, evaluates the
filament model of filament_model.py at the same times, and prints how far apart
they are. The filament model is evaluated twice: exactly in time, which is the
true axisymmetric solution, and with implicit Euler at the finite element time
step, which differs from the finite element model only in space. The comparison
is documented in doc/content/problems/tokamak_wall_induction.md.

    python3 plot_tokamak_wall.py [OutputData/TokamakWall.csv ...] [--figure out.png]

Several CSV files can be given, for example from runs at different mesh
refinements or time steps; each is compared in turn and the first is plotted.
For a run with a vertical displacement event, pass the same
--vertical-displacement as the run was given. The filament model needs scipy;
the figure needs matplotlib.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import filament_model  # noqa: E402

QUANTITIES = {
    "wall_current": ("Wall current", 1e-3, "kA"),
    "wall_current_moment": ("Vertical moment of the wall current", 1e-3, "kA m"),
    "wall_power": ("Ohmic power in the wall", 1e-6, "MW"),
    "bz_centre": ("$B_z$ at the centre, $R = 0$", 1e3, "mT"),
    "bz_outboard": ("$B_z$ at the outboard probe, $R = 1.35$ m", 1e3, "mT"),
}


def read_csv(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


def energy(t, p):
    """Energy dissipated, by the trapezium rule."""
    return float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(t)))


def compare(path, args):
    fe = read_csv(path)
    t = fe["time"]
    models = filament_model.build(args)
    exact = filament_model.reference(args, t, models=models)
    euler = filament_model.reference(args, t, implicit_euler=True, models=models)
    print(f"\n{path}: {len(t) - 1} steps to t = {t[-1] * 1e3:g} ms")
    print(f"  plasma current received / specified, at peak: "
          f"{fe['plasma_current'].max() / args.peak_current:.5f}")
    print(f"  {'quantity':<14}{'peak (FE)':>14}{'vs exact':>11}{'vs implicit Euler':>19}")
    print("  (the largest difference over the run, relative to the peak of the reference)")
    for name, (_, scale, unit) in QUANTITIES.items():
        k = np.argmax(np.abs(fe[name]))
        peak = np.abs(exact[name]).max()
        if peak < 1e-9 * np.abs(exact["wall_current"]).max():
            # Zero by symmetry in the reference; the finite element value is the error.
            print(f"  {name:<14}{fe[name][k] * scale:>10.4g} {unit:<4} (zero by symmetry)")
            continue
        err_exact = np.abs(fe[name] - exact[name]).max() / peak
        err_euler = np.abs(fe[name] - euler[name]).max() / np.abs(euler[name]).max()
        print(f"  {name:<14}{fe[name][k] * scale:>10.4g} {unit:<4}{100 * err_exact:>9.2f}%"
              f"{100 * err_euler:>18.2f}%")
    w_fe, w_ex, w_ie = (energy(t, d["wall_power"]) for d in (fe, exact, euler))
    print(f"  energy dissipated in the wall: {w_fe * 1e-3:.2f} kJ; exact {w_ex * 1e-3:.2f} kJ "
          f"({100 * (w_fe / w_ex - 1):+.2f}%), implicit Euler {w_ie * 1e-3:.2f} kJ "
          f"({100 * (w_fe / w_ie - 1):+.2f}%)")
    return fe


def plot(fe, args, panels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    fe_colour, plasma_colour = "#2a78d6", "#eb6834"
    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
        "xtick.color": muted, "ytick.color": muted, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": grid,
        "grid.linewidth": 0.6, "figure.facecolor": surface, "axes.facecolor": surface,
    })
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True, constrained_layout=True)
    t_ms = fe["time"] * 1e3
    # Dense reference curve, so that the quench is drawn smoothly.
    t_dense = np.linspace(0, fe["time"][-1], 2401)
    dense = filament_model.reference(args, t_dense)
    marker_every = max(1, len(t_ms) // 40)
    for ax, name in zip(axes.flat, panels):
        title, scale, unit = QUANTITIES[name]
        if name == "wall_current":
            ax.plot(t_ms, fe["plasma_current"] * scale, color=plasma_colour, lw=1.5,
                    label="Plasma current (applied)")
        ax.plot(t_dense * 1e3, dense[name] * scale, color=ink, lw=1.2,
                label="Filament model (reference)")
        ax.plot(t_ms, fe[name] * scale, color=fe_colour, lw=1.2, marker="o", ms=4,
                markevery=marker_every, markeredgecolor=surface, markeredgewidth=0.8,
                label="MFEM, tokamak_wall.i")
        ax.set_title(title, loc="left", color=ink, fontsize=10)
        ax.set_ylabel(unit)
        ax.axhline(0, color=muted, lw=0.6)
    for ax in axes[1]:
        ax.set_xlabel("time (ms)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=3, frameon=False)
    fig.savefig(path, dpi=150)
    print(f"\nwrote {path}")


def plot_geometry(args, path):
    """The poloidal cross-section, with the probes and a displaced plasma."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Wedge

    ink, muted, surface = "#0b0b0b", "#52514e", "#fcfcfb"
    wall_colour, plasma_colour = "#52514e", "#eb6834"
    R0, a, b, d = args.major_radius, args.plasma_radius, args.wall_radius, args.wall_thickness
    fig, ax = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    fig.set_facecolor(surface)
    ax.set_facecolor(surface)
    ax.add_patch(Wedge((R0, 0), b + d, 0, 360, width=d, color=wall_colour))
    # The plasma's parabolic current profile, shaded by current density.
    r = np.linspace(0, a, 40)
    for r_out, r_in in zip(r[:0:-1], r[-2::-1]):
        shade = 1 - (0.5 * (r_out + r_in) / a) ** 2
        ax.add_patch(Wedge((R0, 0), r_out, 0, 360, width=r_out - r_in, lw=0,
                           color=plasma_colour, alpha=0.15 + 0.75 * shade))
    ax.add_patch(Circle((R0, 0.08), a, fill=False, ls="--", lw=1.0, color=plasma_colour))
    for name, (R, Z) in filament_model.PROBES.items():
        ax.plot(R, Z, marker="x", color=ink, ms=7, mew=1.6)
    ax.axvline(0, color=muted, lw=0.8, ls="-.")
    ax.annotate("plasma, $r < a$", (R0, -0.05), ha="center", color=ink)
    ax.annotate("wall", (R0 + (b + d) * np.cos(2.3), (b + d) * np.sin(2.3)), xytext=(0.2, 0.62),
                arrowprops=dict(arrowstyle="-", color=muted, lw=0.8), color=ink)
    ax.annotate("vacuum", (R0 + 0.36 * np.cos(-0.9), 0.36 * np.sin(-0.9)), xytext=(1.35, -0.62),
                arrowprops=dict(arrowstyle="-", color=muted, lw=0.8), color=ink)
    ax.annotate("displaced plasma\n(vertical displacement event)", (R0 + a * 0.7, 0.08 + a * 0.71),
                xytext=(1.3, 0.55), arrowprops=dict(arrowstyle="-", color=muted, lw=0.8), color=ink)
    ax.annotate("probe", filament_model.PROBES["bz_outboard"], xytext=(1.55, 0.12), color=ink,
                arrowprops=dict(arrowstyle="-", color=muted, lw=0.8))
    ax.annotate("probe", filament_model.PROBES["bz_centre"], xytext=(0.06, 0.12), color=ink)
    ax.annotate("$z$ axis", (0, -0.62), xytext=(0.04, -0.62), color=muted)
    ax.set_xlim(-0.1, 2.0)
    ax.set_ylim(-0.7, 0.75)
    ax.set_aspect("equal")
    ax.set_xlabel("$R$ (m)")
    ax.set_ylabel("$z$ (m)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", nargs="*", default=["OutputData/TokamakWall.csv"])
    p.add_argument("--vertical-displacement", type=float, default=0.0,
                   help="as given to plasma_current_density.i (m)")
    p.add_argument("--panels", nargs=4, choices=QUANTITIES,
                   default=["wall_current", "wall_power", "bz_centre", "bz_outboard"],
                   help="the quantities to plot")
    p.add_argument("--figure", help="write the comparison figure here")
    p.add_argument("--geometry", help="write a figure of the cross-section here, and stop")
    opts = p.parse_args()
    if opts.geometry:
        plot_geometry(filament_model.make_parser().parse_args([]), opts.geometry)
        return

    args = filament_model.make_parser().parse_args([])
    args.vertical_displacement = opts.vertical_displacement
    results = [compare(path, args) for path in opts.csv]
    if opts.figure:
        plot(results[0], args, opts.panels, opts.figure)


if __name__ == "__main__":
    main()
