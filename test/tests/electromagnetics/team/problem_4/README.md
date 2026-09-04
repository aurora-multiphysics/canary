# TEAM Problem 4 — the FELIX brick

TEAM Problem 4 is an eddy current benchmark from the COMPUMAG TEAM (Testing Electromagnetic Analysis Methods)
workshop series. A hollow aluminium brick sits in a uniform axial magnetic field which is then switched off with
an exponential decay. The decaying field drives eddy currents that circulate around the hole, and the benchmark
quantities are the current crossing the brick wall and the ohmic power the currents dissipate.

The hole makes the conducting region multiply connected, which is what gives the problem its teeth: the
formulation has to support a net current circulating around a loop that cannot be shrunk to a point.

## Problem specification

| Quantity | Value |
| --- | --- |
| Brick, outer | 0.1524 m (x) × 0.1016 m (y) × 0.0508 m (z) |
| Hole | 0.0889 m (x) × 0.0381 m (y), bored through the brick along z |
| Wall thickness | 0.03175 m, uniform in both x and y |
| Conductivity | 2.538 × 10⁷ S/m (resistivity 3.94 × 10⁻⁸ Ω m) |
| Relative permeability | 1.0, everywhere |
| Applied field | `B_ext(t) = B0 exp(-t/tau) z_hat`, `B0` = 0.1 T, `tau` = 0.0119 s |
| Simulated interval | 0 to 20 ms |

The applied field is parallel to the axis of the hole.

## Symmetry and mesh

The geometry and the applied field share the symmetry planes at `x = 0`, `y = 0` and `z = 0`, so `team4_symmetrized.e`
meshes one octant of the problem: the brick corner, surrounded by a 0.508 m cube of free space standing in for the
unbounded space around it.

Blocks:

- `Brick` — the conductor
- `Vacuum` — the surrounding free space

Sidesets:

- `x_symmetry`, `y_symmetry`, `z_symmetry` (ids 4, 3, 2) — the three symmetry planes
- `x_high`, `y_high`, `z_high` (ids 6, 5, 1) — the outer faces of the free space box
- `semi_xsection_in` (id 7) — the cut through the brick wall on the `y = 0` plane, between the hole and the outer face
- `semi_xsection_out` (id 8) — the corresponding cut on the `x = 0` plane

The two `semi_xsection` sidesets are where the circulating current is measured. Each spans only half the brick
thickness, because the octant is cut at `z = 0`, so the current through a full cross-section of the brick is twice
the value reported. Likewise the reported dissipation is that of one octant; the whole brick dissipates eight times
as much.

## Formulation

The total magnetic field is split into the applied background field and the field induced by the eddy currents,
`H = H_ext + H_ind`, and the two parts are computed by two coupled apps.

`team4_external_source_field.i` (sub-app) builds the spatial part of the applied field. Rather than projecting
`B0 z_hat` onto H(curl) directly, it solves a scalar potential problem `div(mu grad(phi)) = 0` with `phi = z`
imposed on the boundaries where `B` is normal, and takes `H_ext = (B0/mu0) grad(phi)`. The result is a discrete
gradient by construction, so its discrete curl vanishes to machine precision and contributes no spurious source
current to the main app. Because the applied field has no spatial time dependence, this solve runs once, at
`INITIAL`, and the result is reused at every step.

`team4_induced_field.i` (main app) solves for `H_ind` in H(curl):

```text
(mu dH_ind/dt, H') + (rho curl(H_ind), curl(H')) = -(dB_ext/dt, H')
```

Since `curl(H_ext) = 0`, the entire current density is `J = curl(H_ind)`, projected onto H(div) by an aux kernel.

Two details are worth knowing:

- **Free space resistivity.** The curl-curl operator is singular wherever the resistivity is zero, so free space is
  given a large but finite resistivity (1 Ω m) rather than an infinite one. At eight orders of magnitude above the
  brick this leaves a negligible current outside the conductor while keeping the system solvable. This is why the
  curl-curl kernel is applied over the whole domain rather than restricted to the `Brick` block.
- **Boundary conditions.** The tangential induced field is set to zero on `z_symmetry` and on the three far-field
  faces. On `z_symmetry` this is exact rather than an approximation: reflecting in `z` makes `H_x` and `H_y` odd, and
  those are precisely the tangential components there. The `x_symmetry` and `y_symmetry` planes are deliberately left
  unconstrained — their tangential components include `H_z`, which is not zero — and pick up the natural condition of
  the weak form, `n × E = 0`, which is exactly the required statement that no current crosses a symmetry plane.

## Files

| File | Role |
| --- | --- |
| `team4_induced_field.i` | Main app: transient eddy current solve. This is the input to run. |
| `team4_external_source_field.i` | Sub-app: the applied background field. Launched automatically. |
| `team4_symmetrized.e` | Octant mesh, 26760 tetrahedra. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

## Running

```sh
<path to a MOOSE app built with MFEM support> -i team4_induced_field.i
```

or, from the repository root, through the test harness:

```sh
./run_tests --re team4
```

Results are written to `OutputData/`: `TEAM4CSV.csv` holds the postprocessor time histories, and `TEAM4/` holds a
ParaView collection containing `induced_h_field`, `external_h_field` and `j_field`. The applied field has no output
of its own; it is transferred into the main app and written out from there.

Each field snapshot is around 18 MB, so the ParaView collection is written only every fifth step, giving five
snapshots across the 20 ms. The postprocessor histories keep the full time resolution. Drop
`time_step_interval` from the output block if you want every step.

## Output

Three postprocessors are reported at each step:

- `OhmicHeating` — `integral of rho |J|^2` over the brick, the power dissipated in the modelled octant, in W
- `InnerLimbCurrent` — current through `semi_xsection_in`, in A
- `OuterLimbCurrent` — current through `semi_xsection_out`, in A

The two limb currents are a useful correctness check as much as a result: the current circulates around the hole, so
whatever enters through one cross-section must leave through the other. They should agree in magnitude and differ in
sign. With the settings shipped here they match to seven significant figures.

The current rises from zero, peaks at around 11 ms and then decays with the applied field.

## Caveats

The mesh and the timestep are chosen to keep the example runnable in well under a minute, not to be converged. The
implicit Euler integrator is noticeably diffusive at this timestep: coarsening `dt` from 0.001 s to 0.005 s lowers
the peak current by over 20%. Refine both before comparing against published benchmark data.

## References

- A. Kameari, *Results for benchmark calculations of problem 4 (the FELIX brick)*, COMPEL, 1988.
- [Problem statement summary][team4-spec] — eddy current simulation of a conducting brick in a decaying
  magnetic field.

[team4-spec]: https://ceae-server.colorado.edu/v2016/books/bmk/ch01s08ach67.html
