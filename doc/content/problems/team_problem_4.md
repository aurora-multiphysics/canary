# TEAM Problem 4: The FELIX Brick

[TEAM](https://www.compumag.org/wp/team/) Problem 4 is an eddy current benchmark from the COMPUMAG workshop
series. A hollow aluminium brick sits in a uniform axial magnetic field which is then switched off with an
exponential decay. The decaying field drives eddy currents that circulate around the hole, and the benchmark
quantities are the current crossing the brick wall and the ohmic power those currents dissipate.

The hole makes the problem non-trivial by making the conducting region multiply connected, so the
formulation has to support a net current circulating around a loop that cannot be shrunk to a point.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

| Quantity | Value |
| :- | :- |
| Brick, outer | 0.1524 m ($x$) by 0.1016 m ($y$) by 0.0508 m ($z$) |
| Hole | 0.0889 m ($x$) by 0.0381 m ($y$), bored through the brick along $z$ |
| Wall thickness | 0.03175 m, uniform in both $x$ and $y$ |
| Conductivity | $2.538 \times 10^{7}$ S/m (resistivity $3.94 \times 10^{-8}\ \Omega$ m) |
| Relative permeability | 1.0, everywhere |
| Applied field | $\vec{B}_\text{ext}(t) = B_0 e^{-t/\tau}\hat{z}$, $B_0 = 0.1$ T, $\tau = 0.0119$ s |
| Simulated interval | 0 to 20 ms |

The applied field is parallel to the axis of the hole.

## Symmetry and Mesh id=mesh

The geometry and the applied field share the symmetry planes at $x = 0$, $y = 0$ and $z = 0$, so only one octant
is meshed: the brick corner, surrounded by a 0.508 m cube of free space standing in for the unbounded space
around it. The supplied mesh has 26760 tetrahedra.

Two blocks are defined, `Brick` for the conductor and `Vacuum` for the surrounding free space. The sidesets are
named for the role they play:

- `x_symmetry`, `y_symmetry`, `z_symmetry` (ids 4, 3, 2) — the three symmetry planes
- `x_high`, `y_high`, `z_high` (ids 6, 5, 1) — the outer faces of the free space box
- `semi_xsection_in` (id 7) — the cut through the brick wall on the $y = 0$ plane, between the hole and the
  outer face
- `semi_xsection_out` (id 8) — the corresponding cut on the $x = 0$ plane

The two `semi_xsection` sidesets are where the circulating current is measured. Each spans only half the brick
thickness, because the octant is cut at $z = 0$, so the current through a full cross-section of the brick is
twice the value reported. Likewise the reported dissipation is that of one octant; the whole brick dissipates
eight times as much.

## Formulation id=formulation

The total magnetic field is split into the applied background field and the field induced by the eddy currents,
$\vec{H} = \vec{H}_\text{ext} + \vec{H}_\text{ind}$, and the two parts are computed by two coupled apps.

### The Applied Field id=applied-field

The sub-app builds the spatial part of the applied field. Rather than projecting $B_0\hat{z}$ onto $H(\text{curl})$
directly, it solves a scalar potential problem

!equation
\nabla \cdot \left( \mu \nabla \phi \right) = 0, \qquad \phi = z \ \text{ where } \vec{B} \text{ is normal},

whose solution is $\phi = z$, and takes $\vec{H}_\text{ext} = (B_0/\mu_0)\nabla\phi$. The result is a discrete
gradient by construction, so its discrete curl vanishes to machine precision and contributes no spurious source
current to the main app. Projecting the analytic field directly would leave a small discrete curl behind, which
the main app would read as a real current.

The $x = 0$ and $y = 0$ symmetry planes are left with the natural condition $\mu\,\partial\phi/\partial n = 0$,
which $\phi = z$ satisfies exactly, so the discrete solution recovers the uniform field to round-off. Because the
applied field has no spatial time dependence, this solve runs once, at `INITIAL`, and the result is reused at
every step.

!listing canary/test/tests/electromagnetics/team/problem_4/team4_external_source_field.i

### The Induced Field id=induced-field

The main app solves for $\vec{H}_\text{ind}$ in $H(\text{curl})$:

!equation
\left( \mu \frac{\partial \vec{H}_\text{ind}}{\partial t}, \vec{H}' \right)
+ \left( \rho \nabla \times \vec{H}_\text{ind}, \nabla \times \vec{H}' \right)
= -\left( \frac{\partial \vec{B}_\text{ext}}{\partial t}, \vec{H}' \right)

Since $\nabla \times \vec{H}_\text{ext} = 0$, the entire current density is
$\vec{J} = \nabla \times \vec{H}_\text{ind}$, projected onto $H(\text{div})$ by an aux kernel.

Two details are worth knowing.

The curl-curl operator is singular wherever the resistivity is zero, so free space is given a large but finite
resistivity ($1\ \Omega$ m) rather than an infinite one. At eight orders of magnitude above the brick this leaves
a negligible current outside the conductor while keeping the system solvable, and it is why the curl-curl kernel
is applied over the whole domain rather than restricted to the `Brick` block.

The tangential induced field is set to zero on `z_symmetry` and on the three far-field faces. On `z_symmetry`
this is exact rather than an approximation: reflecting in $z$ makes $H_x$ and $H_y$ odd, and those are precisely
the tangential components there. The `x_symmetry` and `y_symmetry` planes are deliberately left unconstrained —
their tangential components include $H_z$, which is not zero — and pick up the natural condition of the weak
form, $\vec{n} \times \vec{E} = 0$, which is exactly the required statement that no current crosses a symmetry
plane.

!listing canary/test/tests/electromagnetics/team/problem_4/team4_induced_field.i

## Files id=files

| File | Role |
| :- | :- |
| `team4_induced_field.i` | Main app: the transient eddy current solve. This is the input to run. |
| `team4_external_source_field.i` | Sub-app: the applied background field. Launched automatically. |
| `team4_symmetrized.e` | Octant mesh, 26760 tetrahedra. |
| `plot_team4_current.py` | Draws the current history below from the postprocessor CSV. |
| `run_convergence_study.py` | Repeats the convergence study and extrapolates to zero error. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

All of them live in `test/tests/electromagnetics/team/problem_4`.

## Running id=running

The main app launches the sub-app itself, so only one input needs to be given:

```bash
<path to a MOOSE app built with MFEM support> -i team4_induced_field.i
```

or, from the repository root, through the test harness:

```bash
./run_tests --re TEAM_Problem_4
```

Results are written to `OutputData/`. `TEAM4CSV.csv` holds the postprocessor time histories and `TEAM4/` holds a
ParaView collection containing `induced_h_field`, `external_h_field` and `j_field`. The applied field has no
output of its own; it is transferred into the main app and written out from there.

Each field snapshot is around 18 MB, so the ParaView collection is written only every fortieth step, giving five
across the 20 ms. The postprocessor histories keep the full time resolution. Drop `time_step_interval` from the
output block to get every step.

The run takes about 90 seconds. Adding `Executioner/dt=0.001` on the command line finishes in around 12 seconds,
at the cost of roughly 6% in the peak current — see [Convergence](#convergence).

## Results id=results

Three postprocessors are reported at every step: `OhmicHeating`, the power
$\int \rho \left| \vec{J} \right|^2$ dissipated in the modelled octant, and `InnerLimbCurrent` and
`OuterLimbCurrent`, the currents through the two cross-section sidesets.

The limb currents are a correctness check as much as a result. The current circulates around the hole, so
whatever enters through one cross-section must leave through the other: the two should agree in magnitude and
differ in sign. They match to about 1 part in $10^5$ throughout the run.

The current rises from zero, peaks at 10.9 ms and then decays away with the applied field.

!media media/team4_current.png
       id=current-history
       caption=Total current circulating around the hole, for a full cross-section of the brick. The eddy current
               lags the applied field: it peaks at 10.9 ms, well after the field has begun to collapse, and then
               decays with it.
       style=width:80%;margin-left:auto;margin-right:auto;

The figure is produced by `plot_team4_current.py`, which reads the postprocessor CSV and doubles the reported
limb current to span the full brick thickness. Regenerate it after a run with:

```bash
cd test/tests/electromagnetics/team/problem_4
./plot_team4_current.py
```

It needs only `matplotlib`, and writes `doc/content/media/team4_current.png` by default. It also prints the worst
disagreement between the two limb currents over the run, which is the consistency check described above.

| Quantity | Modelled octant | Whole brick |
| :- | :- | :- |
| Peak current, at 10.9 ms | 1711 A | 3423 A ($\times 2$, full cross-section) |
| Peak ohmic power, at 9.9 ms | 14.17 W | 113.4 W ($\times 8$, all octants) |

## Convergence id=convergence

Time discretisation is first-order and the timestep dominates the error. `run_convergence_study.py` repeats the
study below and reports the extrapolated limit:

```bash
cd test/tests/electromagnetics/team/problem_4
./run_convergence_study.py --exe /path/to/app-opt        # or set MOOSE_APP / MOOSE_DIR
```

It halves the timestep from `--dt0` and prints the table straight into this page's format. Both figures of merit
are quoted for the whole brick, as in the [Results](#results) table's second column: the current is doubled to
span the full cross-section, and the dissipation is summed over all eight octants.

| $\Delta t$ (s) | Steps | Peak current (kA) | Peak power (W) | Error in peak current | Runtime |
| :- | -: | -: | -: | -: | -: |
| $4 \times 10^{-3}$ | 5 | 2.685 | 70.57 | 22.2% | 4 s |
| $2 \times 10^{-3}$ | 10 | 3.028 | 90.03 | 12.3% | 7 s |
| $1 \times 10^{-3}$ | 20 | 3.232 | 101.68 | 6.4% | 13 s |
| $5 \times 10^{-4}$ | 40 | 3.339 | 108.16 | 3.3% | 25 s |
| $2.5 \times 10^{-4}$ | 80 | 3.394 | 111.58 | 1.7% | 47 s |
| $1.25 \times 10^{-4}$ | 160 | 3.423 | 113.36 | 0.8% | 93 s |
| Richardson extrapolation | -- | 3.451 | 115.14 | -- | -- |

The observed order of convergence in the peak current, taken from each triple of successive halvings, is 0.75,
0.93, 0.95 and 0.96 — converging on the expected 1. Errors above are measured against the Richardson
extrapolation to $\Delta t \rightarrow 0$ from the two finest steps, which assumes that first order; pass
`--order` to assume another. The shipped $\Delta t = 1.25 \times 10^{-4}$ s is the coarsest that keeps the peak
current within 1% of the limit, at around 90 seconds. Adding `Executioner/dt=0.001` on the command line finishes
in about 13 seconds at the cost of roughly 6% in the peak.

Spatial error is smaller and does not drive the choice of timestep. The same script runs the mesh study, at a
fixed $\Delta t$:

```bash
./run_convergence_study.py --study mesh --levels 2 --exe /path/to/app-opt
```

| Refinements | Elements | Peak current (kA) | Peak power (W) | Runtime |
| :- | -: | -: | -: | -: |
| 0 | 26760 | 3.232 | 101.68 | 13 s |
| 1 | 214080 | 3.244 | 102.36 | 109 s |
| Richardson extrapolation | -- | 3.257 | 103.05 | -- |

One uniform refinement moves the peak current by 0.37% for eight times the element count, and a first-order
extrapolation from those two levels puts the supplied mesh about 0.8% short of the spatial limit — the same
ballpark as the timestep error at the shipped $\Delta t$, and in the same direction. Two levels fix a magnitude
but not an order; a third would be needed to establish that, and at 64 times the base element count it is a much
larger job than anything else here.

The script refines both apps together. That matters: `MultiAppMFEMCopyTransfer` copies degrees of freedom
directly, so a run that refined only the main app would fail rather than quietly give a wrong answer.

## References id=references

- A. Kameari, *Results for benchmark calculations of problem 4 (the FELIX brick)*, COMPEL, 1988.
- [Problem statement summary](https://ceae-server.colorado.edu/v2016/books/bmk/ch01s08ach67.html) — eddy current
  simulation of a conducting brick in a decaying magnetic field.
