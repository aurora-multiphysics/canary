# TEAM Problem 4: The FELIX Brick

[TEAM](https://www.compumag.org/wp/team/) Problem 4 is an eddy current benchmark from the COMPUMAG workshop
series. A hollow aluminium brick sits in a uniform axial magnetic field which is then switched off with an
exponential decay. The decaying field drives eddy currents that circulate around the hole, and the benchmark
quantities are the current crossing the brick wall and the ohmic power those currents dissipate.

The hole is what gives the problem its teeth. It makes the conducting region multiply connected, so the
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

- `x_symmetry`, `y_symmetry`, `z_symmetry` — the three symmetry planes
- `x_high`, `y_high`, `z_high` — the outer faces of the free space box
- `semi_xsection_in` — the cut through the brick wall on the $y = 0$ plane, between the hole and the outer face
- `semi_xsection_out` — the corresponding cut on the $x = 0$ plane

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

## Running id=running

The main app launches the sub-app itself, so only one input needs to be given:

```bash
<path to a MOOSE app built with MFEM support> -i team4_induced_field.i
```

or, from the repository root, through the test harness:

```bash
./run_tests --re team4
```

Results are written to `OutputData/`. `TEAM4CSV.csv` holds the postprocessor time histories and `TEAM4/` holds a
ParaView collection containing `induced_h_field`, `external_h_field` and `j_field`. The applied field has no
output of its own; it is transferred into the main app and written out from there.

## Results id=results

Three postprocessors are reported at every step: `OhmicHeating`, the power
$\int \rho \left| \vec{J} \right|^2$ dissipated in the modelled octant, and `InnerLimbCurrent` and
`OuterLimbCurrent`, the currents through the two cross-section sidesets.

The limb currents are a correctness check as much as a result. The current circulates around the hole, so
whatever enters through one cross-section must leave through the other: the two should agree in magnitude and
differ in sign. They match to about 1 part in $10^5$ throughout the run.

The current rises from zero, peaks at 10.9 ms and then decays away with the applied field.

| Quantity | Modelled octant | Whole brick |
| :- | :- | :- |
| Peak current, at 10.9 ms | 1711 A | 3423 A ($\times 2$, full cross-section) |
| Peak ohmic power, at 9.9 ms | 14.17 W | 113.4 W ($\times 8$, all octants) |

## Convergence id=convergence

`MFEMTransient` currently offers only implicit Euler — `bdf2` and `crank-nicolson` are rejected — so the time
discretisation is first-order and the timestep dominates the error. Refining it gives textbook first-order
behaviour:

| $\Delta t$ (s) | Steps | Peak current (A) | Peak power (W) | Error in peak current | Runtime |
| :- | -: | -: | -: | -: | -: |
| $5 \times 10^{-3}$ | 4 | 1245.7 | 7.911 | 27.8% | 6 s |
| $2 \times 10^{-3}$ | 10 | 1513.9 | 11.254 | 12.3% | 6 s |
| $1 \times 10^{-3}$ | 20 | 1615.8 | 12.710 | 6.4% | 12 s |
| $5 \times 10^{-4}$ | 40 | 1669.3 | 13.520 | 3.3% | 24 s |
| $2.5 \times 10^{-4}$ | 80 | 1697.0 | 13.947 | 1.7% | 46 s |
| $1.25 \times 10^{-4}$ | 160 | 1711.3 | 14.170 | 0.8% | 92 s |
| Richardson extrapolation | | 1725.6 | 14.392 | | |

Successive differences fall by factors of 1.90, 1.93 and 1.94 as the timestep is halved, giving observed orders
of 0.93, 0.95 and 0.96, converging on the expected 1. Errors above are measured against the Richardson
extrapolation to $\Delta t \rightarrow 0$ taken from the two finest steps. The shipped
$\Delta t = 1.25 \times 10^{-4}$ s is the coarsest that keeps the peak current within 1% of that limit; it takes
around 90 seconds. Adding `Executioner/dt=0.001` on the command line finishes in about 12 seconds at the cost of
roughly 6% in the peak.

Spatial error is much smaller and does not drive the choice. One uniform refinement of the mesh, at fixed
$\Delta t = 10^{-3}$ s, moves the peak current by 0.39% (1615.8 A to 1622.1 A) at eight times the element count,
so the supplied mesh contributes a few tenths of a percent at most. To check this:

```bash
<moose app> -i team4_induced_field.i Mesh/uniform_refine=1 \
  MultiApps/external_source_field/cli_args=Mesh/uniform_refine=1
```

Both meshes have to be refined together, because `MultiAppMFEMCopyTransfer` copies degrees of freedom directly
and requires the two apps to agree exactly.

## References id=references

- A. Kameari, *Results for benchmark calculations of problem 4 (the FELIX brick)*, COMPEL, 1988.
- [Problem statement summary](https://ceae-server.colorado.edu/v2016/books/bmk/ch01s08ach67.html) — eddy current
  simulation of a conducting brick in a decaying magnetic field.
