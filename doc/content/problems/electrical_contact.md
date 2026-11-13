# Electrical Contact: Constriction and Film Resistance of a Contact Spot

Two conductors pressed together do not touch over the whole of their faces. Real surfaces are rough, and
they meet only where the asperities on one meet those on the other, in small spots that Holm called
a-spots. The current between them has to crowd through those spots, which adds a *constriction
resistance*, and the spots are usually covered by oxide or tarnish, which adds a *film resistance*. Together
they are the contact resistance, which sets the heating at every bolted busbar joint, demountable magnet
joint and current lead.

This example models the classical problem of Holm's theory of electric contacts: two coaxial rods, here of
copper and aluminium, touching over a single circular spot on their axis, with a film of given contact
conductance over it. It shows how to set up an electrical contact problem with the MFEM backend of MOOSE:

- an axisymmetric electrostatic solve, with the $2 \pi r$ of the cylindrical volume element in the
  conductivity, as MFEM has no axisymmetric coordinate system of its own;
- a contact conductance at an interface, modelled by a thin conducting layer, as MFEM-MOOSE has no interface
  kernel to impose the jump in potential across it;
- the resistance from the power dissipated, which converges from below and faster than the current through
  the electrodes.

It is verified against a semi-analytic solution of the same problem, independent of the finite element
model, to better than $10^{-6}$ for the default film and $3 \times 10^{-5}$ for a clean contact, and against
the exact one-dimensional solution when the rods touch over their whole faces. It also shows by how much
Holm's formula, the sum of the constriction and film resistances, underestimates the resistance of the
contact: by up to 5% when the two are comparable.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

| Quantity | Value |
| :- | :- |
| Rods | radius $b = 5$ mm, 20 mm from the contact plane to each electrode |
| Upper rod | copper, $\sigma_1 = 5.8 \times 10^7$ S/m |
| Lower rod | aluminium, $\sigma_2 = 3.77 \times 10^7$ S/m |
| Contact spot | radius $a = 1$ mm, on the axis; $a / b = 0.2$ |
| Film over the spot | contact conductance $C = 10^{11}$ S/m$^2$; swept from $10^9$ S/m$^2$ to a clean contact |
| Outside the spot | an insulating gap |
| Electrodes | the end faces of the rods, at 1 mV (copper) and 0 V (aluminium) |
| Reported | resistance between the electrodes; what the contact adds to that of the rods; current density through the spot; potential along the axis |

The contact conductance relates the current density $q$ through the film to the jump in potential across it,
$q = C (\phi_1 - \phi_2)$; its inverse, $1/C$, is the film resistivity in $\Omega$ m$^2$. The sweep spans
$1/C$ from $10^{-16}$ to $10^{-9}$ $\Omega$ m$^2$, from a film that barely registers to one that dominates the
resistance of the contact.

## Theory id=theory

### Constriction and Film Resistance id=holm

Holm's analysis of an a-spot rests on an observation of Maxwell's: an equipotential disk of radius $a$ on
the surface of a half-space of resistivity $\rho$ has a spreading resistance into it of $\rho / (4a)$. Two
half-spaces in contact over the disk, with no film, then have a constriction resistance of

!equation
R_c = \frac{\rho_1 + \rho_2}{4 a},

since the field on each side is that of an equipotential disk, and the current densities through the disk
match. In a rod of finite radius the current has less far to spread, and the constriction resistance is
reduced by a factor $\psi(a/b)$ that depends only on the ratio of the radii,

!equation
R_c = \psi(a/b) \, \frac{\rho_1 + \rho_2}{4 a}, \qquad \psi(0.2) = 0.720876,

the latter from the reference solution below. The film adds its own resistance; if the current crossed the
spot uniformly, that would be $R_f = 1/(C \pi a^2)$. Holm's formula for the resistance of the contact is the
sum,

!equation id=holm-formula
R_\text{Holm} = R_c + R_f,

and the resistance between the electrodes is that of the rods, carrying a uniform current over their length
$l$ from electrode to film, plus the contact's:

!equation
R = R_\text{bulk} + R_\text{contact}, \qquad R_\text{bulk} = (\rho_1 + \rho_2) \frac{l}{\pi b^2}.

### Holm's Formula Is a Lower Bound id=bounds

The two terms of [!eqref](holm-formula) are for different current distributions. The constriction
resistance is that of an equipotential spot, through which the current density rises as
$1/\sqrt{a^2 - r^2}$ towards the edge; the film resistance is that of a uniform current density. A real
contact has one current distribution, which is a compromise between them, and that makes Holm's sum a
lower bound.

The current distribution through a network of resistances is the one that dissipates least power for the
total current (Thomson's theorem). Of all distributions through the spot, the equipotential one dissipates
least in the rods, and the uniform one least in the film. So whatever the actual distribution, its power in
the rods is at least that of the equipotential spot and its power in the film at least that of a uniform
current, and $R_\text{contact} \ge R_c + R_f$. Equally, the uniform current is one possible distribution, so
the resistance it gives is an upper bound. Between the two, the film evens out the current over the spot,
which raises the constriction resistance above $R_c$: the thicker the film, the closer the current comes to
uniform.

### Reference Solution id=reference

`contact_reference.py` solves the problem without finite elements. A current density $q(r)$ through the spot
of a rod of resistivity $\rho$ sets up the potential

!equation
\phi(r, z) = \rho \left( C_0 - \frac{I z}{\pi b^2} + \sum_n c_n J_0(k_n r) e^{-k_n z} \right),
\qquad J_1(k_n b) = 0,

with $c_n$ given by the Hankel transform of $q$ on the zeros of $J_1$, which keep current from crossing the
side of the rod. The power that $q$ dissipates in both rods, over and above a uniform current, and in the
film, is a quadratic form in $q$, and the true $q$ minimises it for the given total current. Minimising it
over a finite family of distributions gives an upper bound on the resistance that converges from above as
the family grows:

- for a clean spot, $q = \sum_m \alpha_m (1 - r^2/a^2)^{m - 1/2}$, which has the singularity of an
  equipotential spot at its edge; eight terms and twelve agree to $10^{-8}$;
- with a film, which that singularity would make dissipate infinite power, $q$ is constant on each of 200
  and then 400 annuli graded towards the edge of the spot, and extrapolated in the width of the annuli.

The Hankel transforms of both are closed form. The Bessel series is summed to 200,000 terms, and its tail is
extrapolated too. The reference is good to about $10^{-7}$.

## Formulation id=formulation

The problem is solved for the electric potential $\phi$ in $H^1$, with $\vec{J} = -\sigma \nabla \phi$ and
$\nabla \cdot \vec{J} = 0$. It is axisymmetric, and is meshed in the $(r, z)$ half-plane with $r$ along $x$ and
$z$ along $y$. MFEM has no axisymmetric coordinate system, so the weak form is written out in full,

!equation
\left( 2 \pi r \, \sigma \nabla \phi, \nabla \phi' \right) = 0,

with the $2 \pi r$ in the conductivity coefficient of
[MFEMDiffusionKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html). The
weighted conductivity of each block is a `ParsedFunction` of $x$, collected into one coefficient by
[MFEMGenericFunctorMaterial](https://mooseframework.inl.gov/source/mfem/functormaterials/MFEMGenericFunctorMaterial.html).
With the weight in it, every integral over the half-plane is an integral over the volume of the rods: the
power reported is that dissipated in the rods, not per radian.

### The Film as a Thin Layer id=thin-layer

Across the film the potential jumps by $q / C$. MFEM-MOOSE has no interface kernel with which to impose a
jump of that kind between the two sides of a shared boundary, like the libMesh
[ElectrostaticContactCondition](https://mooseframework.inl.gov/source/interfacekernels/ElectrostaticContactCondition.html)
does, and a continuous $H^1$ potential could not represent one anyway. Instead the film is meshed as a layer
of elements of thickness $t = 1$ µm over the spot, with conductivity

!equation
\sigma_f = C t.

Its resistance through its thickness is then $t / \sigma_f = 1/C$ per unit area, whatever $t$; the current it
can carry along its length, which a true interface cannot, is $\sigma_f t = C t^2$ per unit width, and vanishes
as $t$ does. At $t / a = 10^{-3}$ the layer stands in for the interface to a few parts in $10^{6}$
([see below](#mesh-convergence)), and the same layer models a clean contact when $C$ is made large enough.
The gap outside the spot, of the same thickness, is not meshed at all: the faces of the rods either side of it
are boundaries of the mesh, and carry no current.

### Resistance from the Power id=power

The resistance is computed from the power dissipated,

!equation
R = \frac{V^2}{P}, \qquad P = \int 2 \pi r \, \sigma \lvert \nabla \phi \rvert^2 \, dA,

integrated by
[MFEMVectorFEInnerProductIntegralPostprocessor](https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorFEInnerProductIntegralPostprocessor.html)
from $\nabla \phi$ held in the $H(\text{curl})$ space of the same order, which holds it exactly. The finite
element solution is the discrete potential that dissipates least power, but the discrete space cannot carry
the current quite as the true potential would, so it dissipates more: the resistance converges from below,
and as the square of the error in the potential. `ContactResistance` subtracts the analytic $R_\text{bulk}$,
leaving what the contact adds. The current through each electrode is also reported, from the flux of
$\nabla \phi$ through it by
[MFEMVectorBoundaryFluxIntegralPostprocessor](https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorBoundaryFluxIntegralPostprocessor.html),
as a check: both match $V / R$ to $3 \times 10^{-8}$.

The flux postprocessors take the weighted conductivity of their block as a function rather than the
material property: the property is defined block by block, and the boundary elements on which the flux is
integrated belong to no block.

## Boundary Conditions and Mesh id=mesh

The electrodes are held at 1 mV on the copper and 0 V on the aluminium. They are four rod radii from the
contact, by which distance the slowest decaying mode of the series above has fallen to $e^{-3.83 \times 4}$
of its value at the contact, so the current through them is uniform. Everywhere else, the axis, the sides of
the rods and their faces either side of the gap, the natural boundary condition applies: no current crosses.

`electrical_contact_mesh.py` builds the mesh with [gmsh](https://gmsh.info). The film is a single row of
quadrilaterals, graded towards the edge of the spot, where the current crowds and, for a clean contact, the
field is singular. The rods are meshed with triangles, which grow from the edge of the spot at a rate of 0.1,
from 0.25 µm, a quarter of the film thickness, to 1 mm. The mesh has 10037 elements, and with second order
elements 20495 unknowns. The blocks are `Copper`, `Aluminium` and `Film`, and the boundaries
`copper_electrode`, `aluminium_electrode`, `axis`, `outer` and `gap`. `electrical_contact_full.e` is the same
mesh with the spot as wide as the rods, for the one-dimensional check.

## Solver id=solver

The system is symmetric positive definite, and is solved with conjugate gradients preconditioned by
BoomerAMG to a relative residual of $10^{-12}$. The conductivity is weighted by $r$, which vanishes on the axis,
and the film's, $C t$, is anywhere from $10^3$ S/m to $10^{10}$ S/m over the sweep, against the rods'
$4 \times 10^7$ S/m to $6 \times 10^7$ S/m. Neither troubles the multigrid much: the solve takes 23 iterations
for films up to $C = 10^{11}$ S/m$^2$, 27 at $10^{13}$ S/m$^2$, and 150 for the clean contact, whose layer is
170 times as conductive as copper. Each takes about a second.

!listing canary/test/tests/electromagnetics/electrical_contact/electrical_contact.i

## Files id=files

| File | Role |
| :- | :- |
| `electrical_contact.i` | The solve. This is the input to run. |
| `electrical_contact_mesh.py` | Builds `electrical_contact.e` with gmsh. Run with `--help` for the meshing controls. |
| `electrical_contact.e` | The mesh, as built by `electrical_contact_mesh.py` with its defaults. |
| `electrical_contact_full.e` | The mesh for rods in contact over their whole faces, built with `--spot-radius 5`. |
| `contact_reference.py` | The semi-analytic reference solution, and the bounds of Holm's formula and a uniform current. |
| `run_contact_sweep.py` | Sweeps the contact conductance, runs the mesh convergence study, and writes the tables and figures below. |

## Running id=running

```bash
cd test/tests/electromagnetics/electrical_contact
moose_test-opt -i electrical_contact.i                              # C = 1e11 S/m^2
moose_test-opt -i electrical_contact.i contact_conductance=1e16     # a clean contact
moose_test-opt -i electrical_contact.i Mesh/file=electrical_contact_full.e spot_radius=5e-3
./contact_reference.py --conductance 1e11                           # the reference
./run_contact_sweep.py --exe moose_test-opt --plot                  # the sweep and figures
```

Each run takes about a second on one core. `OutputData/ElectricalContactCSV.csv` holds the resistance, the
contact resistance, the power and the share of it dissipated in the film, and the electrode currents;
`OutputData/ElectricalContactCSV_spot_current_density_0001.csv` the current density across the middle of the
film; and `OutputData/ElectricalContactCSV_axial_potential_0001.csv` the potential along the axis. The
potential and current density are also written for ParaView to `OutputData/ElectricalContact`.

The reference takes about thirty seconds a conductance; the sweep computes them in parallel and caches them.
To rebuild the mesh, which needs the gmsh Python module and a MOOSE executable to convert it to Exodus:

```bash
./electrical_contact_mesh.py --exe /path/to/moose_test-opt
./electrical_contact_mesh.py --exe /path/to/moose_test-opt --spot-radius 5 \
    --msh electrical_contact_full.msh --exodus electrical_contact_full.e
```

## Results id=results

### Default Contact id=default-results

With $C = 10^{11}$ S/m$^2$ the film's resistance $R_f = 3.18$ µΩ is 0.4 of the constriction resistance
$R_c = 7.89$ µΩ, and the two interact about as strongly as they can. The resistance between the electrodes is
22.7726 µΩ, of which the rods account for 11.1448 µΩ and the contact 11.6278 µΩ, against 11.6278 µΩ from the
reference: they differ by $7.6 \times 10^{-7}$. Holm's formula gives 11.0707 µΩ, which the contact resistance
is 5.0% above. The film
dissipates 14.8% of the power, against the 14.0% its share of the resistance would give it with a uniform
current through it, since more of the current crosses it near the edge.

### Contact Conductance id=sweep

The resistance the contact adds to that of the rods, from a film that dominates it to a clean contact. $R_f /
R_c$ is the film resistance over the clean constriction resistance, 7.8876 µΩ. The clean contact is run with
$C = 10^{16}$ S/m$^2$, whose film adds $3 \times 10^{-11}$ Ω, and is compared with the reference for a clean
spot plus that.

| $C$ (S/m$^2$) | $R_f / R_c$ | MFEM-MOOSE (µΩ) | Reference (µΩ) | Difference | Holm, lower bound (µΩ) | Uniform, upper bound (µΩ) | Film share of power |
| :- | :- | :- | :- | :- | :- | :- | :- |
| $10^{9}$ | 40.4 | 327.06777 | 327.06778 | $-2.2 \times 10^{-8}$ | 326.19745 | 327.07361 | 0.941 |
| $10^{9.5}$ | 12.8 | 109.40405 | 109.40406 | $-6.5 \times 10^{-8}$ | 108.54599 | 109.42214 | 0.835 |
| $10^{10}$ | 4.04 | 40.540693 | 40.540700 | $-1.8 \times 10^{-7}$ | 39.718553 | 40.594708 | 0.617 |
| $10^{10.5}$ | 1.28 | 18.683877 | 18.683885 | $-4.1 \times 10^{-7}$ | 17.953407 | 18.829562 | 0.341 |
| $10^{11}$ | 0.404 | 11.627763 | 11.627771 | $-7.6 \times 10^{-7}$ | 11.070663 | 11.946818 | 0.148 |
| $10^{11.5}$ | 0.128 | 9.240884 | 9.240896 | $-1.3 \times 10^{-6}$ | 8.894149 | 9.770304 | 0.058 |
| $10^{12}$ | 0.0404 | 8.386196 | 8.386218 | $-2.7 \times 10^{-6}$ | 8.205874 | 9.082029 | 0.022 |
| $10^{12.5}$ | 0.0128 | 8.070781 | 8.070833 | $-6.5 \times 10^{-6}$ | 7.988223 | 8.864378 | 0.008 |
| $10^{13}$ | 0.00404 | 7.954018 | 7.954153 | $-1.7 \times 10^{-5}$ | 7.919395 | 8.795550 | 0.003 |
| $10^{13.5}$ | 0.00128 | 7.911135 | 7.911436 | $-3.8 \times 10^{-5}$ | 7.897630 | 8.773785 | 0.001 |
| $10^{14}$ | 0.000404 | 7.895553 | 7.896015 | $-5.9 \times 10^{-5}$ | 7.890748 | 8.766902 | 0.000 |
| $10^{14.5}$ | 0.000128 | 7.890046 | 7.890523 | $-6.0 \times 10^{-5}$ | 7.888571 | 8.764726 | 0.000 |
| $10^{15}$ | 0.0000404 | 7.888189 | 7.888587 | $-5.0 \times 10^{-5}$ | 7.887883 | 8.764038 | 0.000 |
| clean | 0.00000404 | 7.887393 | 7.887596 | $-2.6 \times 10^{-5}$ | 7.887596 | 8.763751 | 0.000 |

!media media/electrical_contact_resistance.png
       id=resistance
       caption=Left: the resistance of the contact against the contact conductance of the film, with the
               reference and the bounds of Holm's formula and a uniform current. Right: the same, over Holm's
               formula, against the ratio of the film resistance to the clean constriction resistance.

MFEM-MOOSE matches the reference to $10^{-6}$ or better wherever the film resistance is at least a tenth of
the constriction resistance, and to $6 \times 10^{-5}$ everywhere. Every difference is negative, as it must
be for a resistance [computed from the power](#power). The largest are where the film is thin enough to let
the current crowd towards the edge of the spot, into a layer of width $1 / (C (\rho_1 + \rho_2))$: 0.2 µm at
$C = 10^{14}$ S/m$^2$, the size of the smallest elements. Thinner films narrow the layer further, until the clean
singularity takes over.

Holm's formula is exact in both limits, where one of the resistances is negligible, and low in between: the
resistance of the contact is up to 5.0% above it where $R_f / R_c \approx 0.4$, close to the default film. The uniform current bound is close to
the reference for films that dominate, but approaches 11% too high for a clean contact, the ratio of the
constriction resistance of a uniform current to that of an equipotential spot in a rod with $a / b = 0.2$.

!media media/electrical_contact_spot.png
       id=spot
       caption=Left: the current density through the middle of the film over its mean, from MFEM-MOOSE
               (circles, every eighth sample) and the reference (lines). Right: the potential along the axis,
               which jumps across the film; the jump looks sloped only because the samples are 0.1 mm apart.

The current density through the spot shows the film at work. With a clean contact it has the
$1/\sqrt{1 - r^2/a^2}$ shape of an equipotential spot, half its mean value at the centre (0.501, against 0.5
for a half-space) and singular at the edge. The film evens it out: at $C = 10^{11}$ S/m$^2$ it is 0.71 of the
mean at the centre and 1.68 at the edge, and at $10^9$ S/m$^2$ uniform to 1%. MFEM-MOOSE reproduces the
reference profiles to within $2 \times 10^{-5}$ for the thickest film, $9 \times 10^{-4}$ for the default, and
$6 \times 10^{-3}$ for a clean contact, over $r < 0.99a$; the current density is projected onto discontinuous
linear elements, which limits its accuracy away from the edge more than the potential does. On the axis, the
potential falls steeply as the current converges on the spot, and jumps across the film by $q(0)/C$.

### Full Contact id=full-contact

With the spot as wide as the rods, $a = b$, there is no constriction, and the problem is one-dimensional: the
rods and the film in series, $R = R_\text{bulk} + 1 / (C \pi b^2)$ exactly, and the potential linear in $z$ in
each rod. At $C = 10^{9}$ S/m$^2$, which makes the film resistance, 12.732 µΩ, comparable to the rods',
MFEM-MOOSE gives a contact resistance of 12.732395447345 µΩ against 12.732395447352 µΩ, and the potential
along the axis to $5 \times 10^{-13}$ of the applied voltage. The exact solution lies in the finite element
space, so this is a check of the thin layer and the axisymmetric weighting rather than of the discretisation.

## Convergence id=convergence

### Mesh id=mesh-convergence

The error in the contact resistance against the reference, varying the growth rate of the elements away from
the edge of the spot and the size of the smallest, at the edge, as a multiple of the film thickness. The
shipped mesh has a growth rate of 0.1 and edge elements a quarter of the film thickness.

| Growth | Edge size / $t$ | Order | Elements | Error, $C = 10^{11}$ S/m$^2$ | Error, clean |
| :- | :- | :- | :- | :- | :- |
| 0.4 | 0.25 | first | 1340 | $-1.3 \times 10^{-2}$ | $-1.9 \times 10^{-2}$ |
| 0.4 | 0.25 | second | 1340 | $-9.1 \times 10^{-5}$ | $-1.7 \times 10^{-4}$ |
| 0.2 | 0.25 | first | 3187 | $-3.8 \times 10^{-3}$ | $-5.7 \times 10^{-3}$ |
| 0.2 | 0.25 | second | 3187 | $-8.5 \times 10^{-6}$ | $-3.8 \times 10^{-5}$ |
| 0.1 | 0.25 | first | 10037 | $-1.0 \times 10^{-3}$ | $-1.6 \times 10^{-3}$ |
| 0.1 | 0.25 | second (shipped) | 10037 | $-7.6 \times 10^{-7}$ | $-2.6 \times 10^{-5}$ |
| 0.05 | 0.25 | first | 34279 | $-2.5 \times 10^{-4}$ | $-4.7 \times 10^{-4}$ |
| 0.05 | 0.25 | second | 34279 | $-2.0 \times 10^{-7}$ | $-2.4 \times 10^{-5}$ |
| 0.1 | 1 | first | 7871 | $-1.0 \times 10^{-3}$ | $-1.9 \times 10^{-3}$ |
| 0.1 | 1 | second | 7871 | $-8.6 \times 10^{-7}$ | $-1.2 \times 10^{-4}$ |
| 0.1 | 0.0625 | first | 12294 | $-9.9 \times 10^{-4}$ | $-1.5 \times 10^{-3}$ |
| 0.1 | 0.0625 | second | 12294 | $-7.4 \times 10^{-7}$ | $-1.7 \times 10^{-6}$ |

Second order elements are three to four orders of magnitude more accurate than first on the same mesh, and
the shipped mesh uses them. With the default film the current density is bounded, and the error falls with
the growth rate and not at all with the edge size. With a clean contact it is set by the edge size, which has
to resolve the singularity: once the growth rate is 0.2 or less, the error stops falling with it and falls
instead with the size of the edge elements, to $2 \times 10^{-6}$ at a sixteenth of the film thickness. That
also bounds the difference between the 1 µm layer and the interface it stands in for, which the reference
models: it is no more than a few parts in $10^6$.

## Notes id=notes

The test suite runs three cases, each in about a second: the default film, a clean contact, and full
contact. It compares the resistance, the power and the electrode currents in each, the current density
through the spot for the default film, and the potential along the axis for full contact. The rest is
`run_contact_sweep.py`, which takes about five minutes, most of it in the reference.

The thin layer is the simplest way to impose a contact conductance with an $H^1$ potential in MFEM, and it
needs nothing beyond the framework: only a mesh with the layer in it. It extends directly to a conductance
that varies over the contact, through a coefficient for $\sigma_f$ that depends on position; for instance on
the contact pressure, through a correlation of the kind used by
[ElectrostaticContactCondition](https://mooseframework.inl.gov/source/interfacekernels/ElectrostaticContactCondition.html).
For a coupled electrothermal problem the power density in the film, $\sigma_f \lvert \nabla \phi \rvert^2$, is
the heat generated at the interface.

## References id=references

- R. Holm, *Electric Contacts: Theory and Application*, 4th ed., Springer (1967) - the theory of
  constriction and film resistance, and the formula they are summed by.
- J. C. Maxwell, *A Treatise on Electricity and Magnetism*, Clarendon Press (1873) - the resistance of an
  equipotential disk on a half-space.
- R. S. Timsit, *Electrical contact resistance: properties of stationary interfaces*, IEEE Trans. Compon.
  Packag. Technol., 22, 85 (1999) - a review of constriction and film resistance.
- G. Cincotti et al., *Modeling of SPS apparatus: temperature, current and strain distribution with no
  powders*, AIChE J., 53, 703 (2007) - contact conductances between electrodes, used by
  ElectrostaticContactCondition.
