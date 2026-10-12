"""
Compares FREEGSnke's vessel response during the VDE of run_vde.py with MFEM
solutions of torus_vessel.i driven by the exported source field.

Both models describe the same vessel, so its total toroidal current and the
vertical force on it should agree, up to FREEGSnke's filament and vessel-mode
truncation and MFEM's mesh and finite air region. The vertical force is the load
a VDE puts on the vessel.

Usage: python compare.py FREEGSNKE_CSV MFEM_CSV [MFEM_CSV ...]
                         [--labels LABEL ...] [--plot PNG]

Several MFEM runs, for example on meshes with different air regions, can be
compared at once; --labels names them in the output and the plot.
"""

import argparse
import csv

import numpy as np


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


# Quantity name, FREEGSnke column and MFEM postprocessor.
QUANTITIES = (
    ("vessel toroidal current [A]", "vessel_current", "toroidal_current"),
    ("vessel vertical force [N]", "vessel_vertical_force", "vertical_force"),
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("freegsnke_csv")
    parser.add_argument("mfem_csv", nargs="+")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--plot")
    args = parser.parse_args()
    labels = args.labels or [f"MFEM {i + 1}" for i in range(len(args.mfem_csv))]
    if len(labels) != len(args.mfem_csv):
        raise ValueError("Give one label per MFEM run.")

    freegsnke = read(args.freegsnke_csv)
    runs = []
    for label, path in zip(labels, args.mfem_csv):
        mfem = read(path)
        # Both histories start with a row at the initial time.
        n = min(len(freegsnke["time"]), len(mfem["time"]))
        if not np.allclose(freegsnke["time"][:n], mfem["time"][:n], atol=1e-9):
            raise ValueError(f"{path} is not at the FREEGSnke times.")
        runs.append((label, mfem, n))

    for label, mfem, n in runs:
        print(f"== {label}")
        for name, f_key, m_key in QUANTITIES:
            a, b = freegsnke[f_key][:n], mfem[m_key][:n]
            print(f"{name}: max |FREEGSnke - MFEM| / max |FREEGSnke| = "
                  f"{np.abs(a - b).max() / np.abs(a).max():.3g}")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
        axes[0].plot(freegsnke["time"], freegsnke["axis_z"] * 100, color="k")
        axes[0].set_ylabel("Magnetic axis Z [cm]")
        for ax, (name, f_key, m_key) in zip(axes[1:], QUANTITIES):
            ax.plot(freegsnke["time"], freegsnke[f_key], color="k", label="FREEGSnke")
            for label, mfem, n in runs:
                ax.plot(mfem["time"][:n], mfem[m_key][:n], linestyle="--", label=label)
            ax.set_ylabel(name)
            ax.legend()
        axes[-1].set_xlabel("Time [s]")
        for ax in axes:
            ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=120)


if __name__ == "__main__":
    main()
