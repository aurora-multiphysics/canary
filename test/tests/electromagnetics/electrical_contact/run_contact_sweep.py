#!/usr/bin/env python3
"""Sweep the electrical contact example over the contact conductance and compare.

Runs electrical_contact.i at a range of contact conductances C, from a film that
dominates the resistance of the contact to one that barely registers, and a clean
contact, and compares the resistance the contact adds to that of the rods with:

  - the semi-analytic solution of contact_reference.py;
  - Holm's formula, the constriction resistance of an equipotential spot plus the
    resistance 1 / (C pi a^2) of a film crossed uniformly, which is a lower bound;
  - the resistance of a uniform current density through the spot and film, an
    upper bound.

It also compares the current density through the spot with that of the reference,
and with --convergence runs the default and clean contacts on a sequence of meshes.
The results are printed as the markdown tables in
doc/content/problems/electrical_contact.md, and the figures there are written with
--plot.

The reference takes some thirty seconds a conductance, and is cached in
OutputData/sweep/reference.json; the finite element runs take a second each.

Usage:
    ./run_contact_sweep.py --exe /path/to/moose_test-opt --plot
    ./run_contact_sweep.py --exe /path/to/moose_test-opt --convergence
    ./run_contact_sweep.py --reuse --plot     # re-analyse without re-running
"""

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import contact_reference as ref

HERE = Path(__file__).resolve().parent
MEDIA = HERE.parents[3] / "doc" / "content" / "media"
INPUT = HERE / "electrical_contact.i"
OUT = HERE / "OutputData" / "sweep"

# As in electrical_contact.i.
SPOT_RADIUS = 1e-3
ROD_RADIUS = 5e-3
SIGMA = (5.8e7, 3.77e7)
RHO = 1 / SIGMA[0] + 1 / SIGMA[1]
APPLIED_VOLTAGE = 1e-3
DEFAULT_CONDUCTANCE = 1e11

# The finite element model stands in a conductance high enough for its film to
# be indistinguishable from none for a clean contact; it adds 3e-11 ohm.
CLEAN = 1e16
CONDUCTANCES = [10 ** (e / 2) for e in range(18, 31)] + [CLEAN]
PROFILE_CONDUCTANCES = [1e9, 1e10, 1e11, 1e12, 1e13, CLEAN]


def find_executable(explicit):
    for candidate in (explicit, os.environ.get("MOOSE_APP")):
        if candidate:
            found = shutil.which(str(candidate)) or (candidate if Path(candidate).is_file() else None)
            if found:
                return found
            sys.exit(f"error: executable not found: {candidate}")
    moose_dir = os.environ.get("MOOSE_DIR")
    if moose_dir and (Path(moose_dir) / "test" / "moose_test-opt").is_file():
        return str(Path(moose_dir) / "test" / "moose_test-opt")
    sys.exit("error: no MOOSE application given; pass --exe, or set MOOSE_APP or MOOSE_DIR")


def tag(conductance, label=None):
    return f"C{conductance:.3g}" + (f"_{label}" if label else "")


def run(exe, conductance, label=None, mesh=None, extra=()):
    """Run one case, one at a time; the solve is small and serial."""
    base = OUT / tag(conductance, label)
    cmd = [exe, "-i", INPUT.name, f"contact_conductance={conductance:.17g}",
           "Outputs/ParaViewDataCollection/enable=false",
           f"Outputs/ReportedPostprocessors/file_base={base}"]
    if mesh:
        cmd.append(f"Mesh/file={mesh}")
    cmd += list(extra)
    with open(f"{base}.log", "w") as log:
        subprocess.run(cmd, cwd=HERE, check=True, stdout=log, stderr=subprocess.STDOUT)


def results(conductance, label=None):
    with open(OUT / f"{tag(conductance, label)}.csv") as fh:
        return {k: float(v) for k, v in list(csv.DictReader(fh))[-1].items()}


def profile(conductance):
    with open(OUT / f"{tag(conductance)}_spot_current_density_0001.csv") as fh:
        rows = list(csv.DictReader(fh))
    return (np.array([float(r["x_0"]) for r in rows]),
            np.array([float(r["current_density_0"]) for r in rows]))


def axial_potential(conductance):
    with open(OUT / f"{tag(conductance)}_axial_potential_0001.csv") as fh:
        rows = list(csv.DictReader(fh))
    return (np.array([float(r["x_1"]) for r in rows]),
            np.array([float(r["potential_0"]) for r in rows]))


def reference_at(conductance):
    """R - R_bulk of the reference and of the upper bound at one conductance."""
    a, b = SPOT_RADIUS, ROD_RADIUS
    exact = None if conductance == CLEAN else ref.spot_resistance(a, b, RHO, conductance)
    return exact, ref.uniform_spot(a, b, RHO, conductance, 200000)


def reference(conductances, reuse):
    """R - R_bulk of the reference and of its bounds at each conductance, cached.
    The conductances missing from the cache are computed in parallel."""
    path = OUT / "reference.json"
    cache = json.loads(path.read_text()) if path.exists() else {}
    a = SPOT_RADIUS
    missing = [c for c in conductances if f"{c:.3g}" not in cache]
    if missing and reuse:
        sys.exit(f"error: no cached reference for C = {missing[0]:.3g}; run without --reuse")
    if "clean" not in cache:
        cache["clean"] = ref.spot_resistance(a, ROD_RADIUS, RHO)
    if missing:
        print(f"computing the reference at {len(missing)} conductances", flush=True)
        workers = min(len(missing), max(1, (os.cpu_count() or 2) // 2))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for c, (exact, upper) in zip(missing, pool.map(reference_at, missing)):
                film = 1 / (c * math.pi * a**2)
                cache[f"{c:.3g}"] = dict(exact=cache["clean"] + film if exact is None else exact,
                                         lower=cache["clean"] + film, upper=upper)
        path.write_text(json.dumps(cache, indent=1))
    return cache


def sweep_table(cache):
    clean = cache["clean"]
    print(f"\nClean constriction resistance {clean:.6e} ohm, psi = "
          f"{clean / (RHO / (4 * SPOT_RADIUS)):.6f}\n")
    print("| $C$ (S/m$^2$) | $R_f / R_c$ | $R - R_\\text{bulk}$, MFEM ($\\mu\\Omega$) "
          "| Reference ($\\mu\\Omega$) | Difference | Holm, lower bound ($\\mu\\Omega$) "
          "| Uniform, upper bound ($\\mu\\Omega$) | Film share of power |")
    print("| - | - | - | - | - | - | - | - |")
    for c in CONDUCTANCES:
        r = results(c)
        e = cache[f"{c:.3g}"]
        film = 1 / (c * math.pi * SPOT_RADIUS**2)
        label = "clean" if c == CLEAN else f"{c:.3g}"
        print(f"| {label} | {film / clean:.3g} | {r['ContactResistance'] * 1e6:.6f} "
              f"| {e['exact'] * 1e6:.6f} "
              f"| {(r['ContactResistance'] - e['exact']) / e['exact']:+.1e} "
              f"| {e['lower'] * 1e6:.6f} | {e['upper'] * 1e6:.6f} "
              f"| {r['FilmPowerFraction']:.4f} |")


def profile_table():
    print("\nCurrent density through the spot, max |q_MFEM / q_ref - 1| over r < 0.99 a\n")
    a, b = SPOT_RADIUS, ROD_RADIUS
    for c in PROFILE_CONDUCTANCES:
        x, j = profile(c)
        current = APPLIED_VOLTAGE / results(c)["Resistance"]
        r, q = ref.spot_current_density(a, b, RHO, math.inf if c == CLEAN else c)
        qi = np.interp(x, r, q * current)
        inner = x < 0.99 * a
        print(f"C = {c:.3g}: {np.max(np.abs(j[inner] / qi[inner] - 1)):.2e}")


def plot(cache):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    MEDIA.mkdir(parents=True, exist_ok=True)
    a, b = SPOT_RADIUS, ROD_RADIUS
    clean = cache["clean"]
    filmed = [c for c in CONDUCTANCES if c != CLEAN]
    ratio = np.array([1 / (c * math.pi * a**2) / clean for c in filmed])
    fe = np.array([results(c)["ContactResistance"] for c in filmed])
    exact = np.array([cache[f"{c:.3g}"]["exact"] for c in filmed])
    lower = np.array([cache[f"{c:.3g}"]["lower"] for c in filmed])
    upper = np.array([cache[f"{c:.3g}"]["upper"] for c in filmed])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.loglog(filmed, exact * 1e6, "k-", label="Reference")
    ax1.loglog(filmed, lower * 1e6, "C0--", label="Holm's formula (lower bound)")
    ax1.loglog(filmed, upper * 1e6, "C3:", label="Uniform current (upper bound)")
    ax1.loglog(filmed, fe * 1e6, "o", mfc="none", color="C1", label="MFEM-MOOSE")
    ax1.axhline(clean * 1e6, color="0.6", lw=0.8)
    ax1.text(filmed[-1], clean * 1e6 * 1.15, "clean spot", ha="right", color="0.4")
    ax1.set_xlabel("Contact conductance of the film, $C$ (S/m$^2$)")
    ax1.set_ylabel(r"$R - R_\mathrm{bulk}$ ($\mu\Omega$)")
    ax1.legend(fontsize=8)
    ax1.grid(True, which="both", alpha=0.3)

    holm = lower
    ax2.semilogx(ratio, exact / holm, "k-", label="Reference")
    ax2.semilogx(ratio, upper / holm, "C3:", label="Uniform current (upper bound)")
    ax2.semilogx(ratio, fe / holm, "o", mfc="none", color="C1", label="MFEM-MOOSE")
    ax2.axhline(1, color="C0", ls="--", label="Holm's formula (lower bound)")
    ax2.set_xlabel(r"Film over constriction resistance, $R_f / R_c$")
    ax2.set_ylabel("Contact resistance over Holm's formula")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(MEDIA / "electrical_contact_resistance.png", dpi=150)
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, c in enumerate(PROFILE_CONDUCTANCES):
        x, j = profile(c)
        current = APPLIED_VOLTAGE / results(c)["Resistance"]
        mean = current / (math.pi * a**2)
        r, q = ref.spot_current_density(a, b, RHO, math.inf if c == CLEAN else c)
        label = "clean" if c == CLEAN else f"$C$ = {c:.0e} S/m$^2$"
        ax1.plot(r / a, q * math.pi * a**2, "-", color=f"C{i}", lw=1)
        ax1.plot(x[::8] / a, j[::8] / mean, "o", ms=3, mfc="none", color=f"C{i}", label=label)
    ax1.set_xlabel("$r / a$")
    ax1.set_ylabel(r"Current density through the spot, $|J| / \bar{J}$")
    ax1.set_ylim(0, 5)
    ax1.legend(fontsize=8, title="MFEM-MOOSE (o), reference (-)", title_fontsize=8)
    ax1.grid(True, alpha=0.3)

    for i, c in enumerate(PROFILE_CONDUCTANCES):
        z, phi = axial_potential(c)
        label = "clean" if c == CLEAN else f"$C$ = {c:.0e} S/m$^2$"
        ax2.plot(z * 1e3, phi * 1e3, "-", color=f"C{i}", label=label)
    ax2.set_xlim(-4, 4)
    ax2.set_xlabel("Distance along the axis from the contact, $z$ (mm)")
    ax2.set_ylabel(r"Potential on the axis, $\phi$ (mV)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(MEDIA / "electrical_contact_spot.png", dpi=150)
    plt.close(fig)
    print(f"\nwrote figures to {MEDIA}")


# Meshes for the convergence study: the growth rate of the elements away from the
# edge of the spot, and the size of the smallest, at the edge, in film thicknesses.
CONVERGENCE_MESHES = [(0.4, 0.25), (0.2, 0.25), (0.1, 0.25), (0.05, 0.25),
                      (0.1, 1.0), (0.1, 0.0625)]


def convergence(exe, python, reuse):
    cache = reference([DEFAULT_CONDUCTANCE, CLEAN], reuse)
    rows = []
    for growth, edge in CONVERGENCE_MESHES:
        mesh = OUT / f"mesh_g{growth:g}_e{edge:g}.e"
        for order in ("FIRST", "SECOND"):
            extra = [f"FESpaces/H1FESpace/fec_order={order}",
                     f"FESpaces/HCurlFESpace/fec_order={order}"]
            label = f"{mesh.stem}_{order.lower()}"
            if not reuse:
                if not mesh.exists():
                    subprocess.run([python, str(HERE / "electrical_contact_mesh.py"),
                                    "--growth", str(growth), "--edge-size", str(edge),
                                    "--msh", str(mesh.with_suffix(".msh")), "--exodus",
                                    str(mesh), "--exe", exe],
                                   check=True, stdout=subprocess.DEVNULL)
                for c in (DEFAULT_CONDUCTANCE, CLEAN):
                    print(f"running C = {c:.3g}, {label}", flush=True)
                    run(exe, c, label, mesh, extra)
            with open(mesh.with_suffix(".msh")) as fh:
                lines = fh.read().split("$Elements\n")[1].split("\n")
            elements = sum(1 for ln in lines[1:int(lines[0]) + 1] if ln.split()[1] in ("2", "3"))
            errors = []
            for c in (DEFAULT_CONDUCTANCE, CLEAN):
                exact = cache[f"{c:.3g}"]["exact"]
                errors.append((results(c, label)["ContactResistance"] - exact) / exact)
            rows.append((growth, edge, order, elements, *errors))
    print("\n| Growth | Edge size / $t$ | Order | Elements | Error, $C$ = 10$^{11}$ S/m$^2$ "
          "| Error, clean |")
    print("| - | - | - | - | - | - |")
    for g, e, o, n, err_film, err_clean in rows:
        print(f"| {g:g} | {e:g} | {o.lower()} | {n} | {err_film:+.1e} | {err_clean:+.1e} |")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", help="MOOSE application (default: $MOOSE_APP or "
                                 "$MOOSE_DIR/test/moose_test-opt)")
    p.add_argument("--python", default=sys.executable,
                   help="python with gmsh, for the meshes of the convergence study")
    p.add_argument("--reuse", action="store_true",
                   help="re-analyse the results of an earlier sweep without re-running")
    p.add_argument("--plot", action="store_true", help="write the figures")
    p.add_argument("--convergence", action="store_true",
                   help="run the mesh convergence study instead of the sweep")
    args = p.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    exe = None if args.reuse else find_executable(args.exe)
    if args.convergence:
        convergence(exe, args.python, args.reuse)
        return
    if not args.reuse:
        for c in CONDUCTANCES:
            print(f"running C = {c:.3g}", flush=True)
            run(exe, c)
    cache = reference(CONDUCTANCES, args.reuse)
    sweep_table(cache)
    profile_table()
    if args.plot:
        plot(cache)


if __name__ == "__main__":
    main()
