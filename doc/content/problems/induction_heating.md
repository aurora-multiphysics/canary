# Induction Heating of a Water-Cooled Billet

A three-way coupled problem, run as a MOOSE MultiApp, in which each physics is solved by the code
best suited to it and the fields that link them are exchanged every timestep.

A solenoid coil carrying alternating current induces eddy currents in a cylindrical steel billet
sitting inside it. The eddy currents dissipate a few kilowatts in a skin a centimetre thick at the
billet surface, which heats it; an annulus of water flowing along the billet carries that heat away.
The two couplings run in opposite directions and close a loop: the electromagnetic problem sets the
heat source for the thermal problem, and the billet temperature sets the electrical conductivity the
electromagnetic problem sees, because the conductivity of steel falls as it heats.

| Physics | Solved by | Input |
| :- | :- | :- |
| Eddy currents in the billet, the coil and the air around them | MFEM, $H(\mathrm{curl})$ | `induction_em.i` |
| Heat conduction in the billet | libMesh, [heat transfer module](https://mooseframework.inl.gov/modules/heat_transfer/index.html) | `induction_heating.i` |
| The coolant channel | libMesh, [thermal hydraulics module](https://mooseframework.inl.gov/modules/thermal_hydraulics/index.html) | `induction_coolant.i` |

The thermal problem is the parent app and the other two are its sub-apps:

```
            +-----------------------------------------+
            |          induction_heating.i            |
            |      billet temperature, libMesh        |
            +-----------------------------------------+
                |   ^                     |   ^
     temperature|   |joule_heating T_wall |   |T_fluid, htc
                v   |                     v   |
     +----------------------+   +--------------------------+
     |    induction_em.i    |   |    induction_coolant.i   |
     |    eddy currents,    |   |    coolant channel,      |
     |         MFEM         |   |   thermal hydraulics     |
     +----------------------+   +--------------------------+
```

!alert note title=This example needs more than the framework
Unlike the rest of Canary, this problem uses the heat transfer and thermal hydraulics modules as
well as MFEM, so it cannot run against the `moose_test` executable, which carries MFEM but no
modules. The directory carries a `testroot` naming the `thermal_hydraulics` application, the
smallest one that has everything: it pulls in heat transfer and fluid properties as dependencies,
and MFEM comes from the framework. Build it with

```bash
cd $MOOSE_DIR/modules/thermal_hydraulics && make -j 8
```

and run the suite as `./run_tests --use-subdir-exe`; everything else in Canary continues to use
`moose_test`. `combined-opt` works equally well if that is the build you have.

## Geometry and Mesh id=mesh

Everything is axisymmetric about $z$, although it is solved in full three dimensions because the
$H(\mathrm{curl})$ discretisation is Cartesian.

| Region | Extent |
| :- | :- |
| Billet | $r \le 20$ mm, $\left|z\right| \le 100$ mm |
| Coil | $30 \le r \le 45$ mm, $\left|z\right| \le 50$ mm |
| Air | out to $r = 100$ mm and $\left|z\right| \le 150$ mm |

`mesh.i` builds both meshes the example needs from one common cross section, so that the billet
occupies exactly the same space in each. A
[ConcentricCircleMeshGenerator](https://mooseframework.inl.gov/source/meshgenerators/ConcentricCircleMeshGenerator.html)
lays down four concentric rings, which an
[AdvancedExtruderGenerator](https://mooseframework.inl.gov/source/meshgenerators/AdvancedExtruderGenerator.html)
turns into four nested cylinders; the billet and the coil are then carved out of their own rings
over the axial extent each occupies, and everything left over is air. `preserve_volumes = false`
puts the nodes exactly on the nominal radii rather than pushing them out to preserve the area of the
circle they approximate, which is what makes the two meshes agree on where the billet surface is.

!listing canary/test/tests/multiphysics/induction_heating/mesh.i

The electromagnetic mesh is the whole thing, 3000 hexes; the thermal mesh is the billet alone, 1056
of them, with its curved side and its two flat ends in sidesets. Neither is committed, because
Canary ignores `*.e`; run

```bash
thermal_hydraulics-opt -i mesh.i --mesh-only induction_em.e
thermal_hydraulics-opt -i mesh.i Mesh/final_generator=billet --mesh-only induction_billet.e
```

before the first solve, or let the test suite do it.

The coil is a true annulus rather than a square one, which matters: the source current follows the
azimuthal unit vector $\hat{\varphi} = (-y, x, 0)/r$, and that field is divergence free on an
annulus but not on a square ring. A source current with a non-zero divergence is not consistent with
the vector potential formulation below.

## The Electromagnetic Problem id=electromagnetics

### Formulation

Writing the magnetic vector potential as $\mathbf{A}(\mathbf{x}, t) = \mathrm{Re}\{(\mathbf{A}_r + j
\mathbf{A}_i) e^{j \omega t}\}$, the eddy current A-formulation is

!equation
\nabla \times \left( \nu \nabla \times \mathbf{A} \right) + j \omega \sigma \mathbf{A} = \mathbf{J}_s

with $\nu = 1/\mu$ the magnetic reluctivity and $\mathbf{J}_s$ the source current density prescribed
in the coil. This splits into a pair of real equations for the in-phase and quadrature parts,

!equation
\nabla \times \left( \nu \nabla \times \mathbf{A}_r \right) - \omega \sigma \mathbf{A}_i = \mathbf{J}_s,
\qquad
\nabla \times \left( \nu \nabla \times \mathbf{A}_i \right) + \omega \sigma \mathbf{A}_r = 0 ,

each discretised in first order Nédélec elements. The diagonal blocks are
`MFEMCurlCurlKernel`s and the off-diagonal loss terms are `MFEMMixedVectorMassKernel`s, which
couple the two fields into one block system.

The electric field is $\mathbf{E} = -j \omega \mathbf{A}$, so the current density induced in the
billet is $\mathbf{J} = \sigma \mathbf{E}$ and the dissipation averaged over a cycle is

!equation
q = \frac{\left|\mathbf{J}\right|^2}{2\sigma}
  = \frac{\sigma \left|\mathbf{E}\right|^2}{2}
  = \frac{\sigma \omega^2}{2} \left( \left|\mathbf{A}_r\right|^2 + \left|\mathbf{A}_i\right|^2 \right) ,

which an `MFEMParsedFunction` assembles from the `A_real_mag` and `A_imag_mag` coefficients that
MFEM declares alongside each vector variable, and an `MFEMScalarProjectionAux` projects onto a
constant $L^2$ variable. That variable is what is handed to the thermal app.

### Why the split form id=split

MFEM can express this directly with an `MFEMComplexVariable` and `MFEMComplexKernel`, which is the
more natural way to write a time-harmonic problem and is what the framework's own complex tests do.
It cannot be used here, for two reasons that compound:

- A complex `MFEMProblem` may not hold real variables at all. `ComplexEquationSystem::Init` rejects
  the mixture outright, so the temperature coming in and the dissipation going out would both have
  to be complex-valued as well.
- `MultiAppMFEMTolibMeshShapeEvaluationTransfer` and
  `MultiApplibMeshToMFEMShapeEvaluationTransfer` reach for real grid functions by name and do not
  know about complex ones. (The MFEM-to-MFEM transfers do.)

Splitting the problem by hand into real and imaginary parts sidesteps both: every variable in the
app is then an ordinary real `MFEMVariable` and the transfers work. It solves the same equations
either way, at the same cost, since the complex path assembles the same $2 \times 2$ real block
operator internally.

### Gauge regularisation id=gauge

The curl-curl operator has a large gradient nullspace wherever $\sigma = 0$, which would leave the
system singular in the air. A small artificial conductivity there regularises it: 1 S/m, against
1.4 MS/m in the billet. Its skin depth is several metres, so it does not perturb the solution.

The coil carries the same regularisation conductivity, which is to say it is treated purely as a
prescribed current source: no eddy currents are induced in it and none of its own resistive loss is
modelled. That is the usual treatment for a driven coil, and it is why the coil is not part of the
thermal mesh.

The block operator the two Nédélec fields assemble into is not one the $H(\mathrm{curl})$ algebraic
multigrid preconditioners (AMS) are set up for, so it is solved directly, with `MFEMMUMPS`. On this
mesh that takes about two seconds and a gigabyte of host memory.

`MFEMSuperLU` gives a bit-identical answer and is equally quick, but it is the wrong choice on a
machine with a GPU. Where PETSc is configured `--with-cuda`, its SuperLU_DIST puts the factorisation
on the device and takes roughly 700 MB of GPU memory per MPI rank; ten ranks then exhaust an 8 GB
card and the run dies in PETSc's device initialisation with `cudaErrorMemoryAllocation`. MUMPS stays
on the host and takes none, and runs unchanged at sixteen ranks. This is worth knowing more
generally: a direct solver that silently offloads is easy to mistake for a solver that simply needs
less memory than the alternative.

!listing canary/test/tests/multiphysics/induction_heating/induction_em.i

### Skin effect id=skin

At $f = 3$ kHz with $\sigma = 1.4$ MS/m the plane-wave skin depth is

!equation
\delta = \sqrt{\frac{2}{\omega \mu_0 \sigma}} = 7.8\ \text{mm} ,

against a billet radius of 20 mm, so the heating is strongly concentrated near the surface. Sampling
the dissipation along a radius of the standalone solve at $T = 300$ K gives

| $r$ (mm) | $q$ (W/m³) |
| -: | -: |
| 0.0 | 5.92e5 |
| 4.0 | 3.76e6 |
| 8.0 | 1.26e7 |
| 12.0 | 3.12e7 |
| 16.0 | 6.96e7 |

Fitting $q \propto e^{2(r-a)/\delta}$ to the outer points gives $\delta = 9.7$ mm rather than 7.8 mm.
Two things account for the difference, and both are worth knowing about before reading anything into
a number from this mesh: the plane-wave result does not hold on a cylinder only 2.6 skin depths
across, and the constant $L^2$ elements the dissipation is projected onto are 4 to 5 mm wide, which
is comparable to $\delta$ itself. Refining the billet is the first thing to do with this example if
the profile, rather than the coupling, is what is of interest.

Over the billet the mean dissipation is $2.20 \times 10^7$ W/m³, which over its volume of
$2.485 \times 10^{-4}$ m³ is 5.47 kW.

### The fields id=fields

The app also carries E and B as auxvariables, and reports the cycle root-mean-square magnitude of
each, averaged over the billet. Since $\mathbf{E} = -j \omega \mathbf{A}$, the in-phase part of E
is built from the quadrature part of A and vice versa, which is a scaled sum
(`MFEMSumAux`); $\mathbf{B} = \nabla \times \mathbf{A}$ takes the potential from $H(\mathrm{curl})$
into $H(\mathrm{div})$ (`MFEMCurlAux`).

The two RMS scalars the postprocessors average are not built from those auxvariables, though. For a
phasor field the mean of $|\mathbf{F}|^2$ over a cycle is $(|\mathbf{F}_r|^2 + |\mathbf{F}_i|^2)/2$,
and MFEM already declares `<name>_mag` for a vector variable and `<name>_curl_mag` for the magnitude
of its curl, so both read directly off the solution. That matters: the `coefficient` parameter of
`MFEMScalarProjectionAux` is not registered as a dependency, so an auxvariable reached through it
would not be ordered against the AuxKernel that fills it, and could be a step stale.

Sampling $|\mathbf{B}|_{\rm rms}$ along a radius at the axial midplane:

| $r$ (mm) | 0 | 8 | 12 | 16 | 19.9 | 21 | 25 | 29 |
| :- | -: | -: | -: | -: | -: | -: | -: | -: |
| $\left|\mathbf{B}\right|_{\rm rms}$ (T) | 0.038 | 0.041 | 0.051 | 0.072 | 0.096 | 0.096 | 0.094 | 0.093 |

Two things are worth checking in that profile. It is continuous across the billet surface at 20 mm,
0.0956 T just inside against 0.0958 T just outside, as it must be for a non-magnetic conductor. And
the field does not fall to zero on the axis, because $a/\delta = 2.6$ is not large enough to screen
the billet completely.

The bore field is also the one number in this problem with a clean analytic counterpart. An infinite
solenoid of the same ampere-turns per unit length would give
$\mu_0 N I / L = 0.133$ T rms. This coil is 0.1 m long and 0.09 m across, and a finite solenoid
gives a fraction $L / \sqrt{L^2 + D^2} = 0.74$ of the infinite result, so 0.098 T. The solve gives
0.096 T, within 3%, which is a reasonable check that the coil drive, the gauge regularisation and
the truncation boundary are all doing what they should.

The mean electric field over the billet is 3.09 V/m rms. That is consistent with the dissipation
through $q = \sigma |\mathbf{E}|^2_{\rm rms}$: the mean dissipation implies
$\sqrt{\bar{q}/\sigma} = 3.97$ V/m, and the mean of $|\mathbf{E}|$ is necessarily below the root of
the mean of $|\mathbf{E}|^2$ for a field that varies as strongly across the billet as this one does.

## The Thermal Problem id=thermal

The parent app conducts heat through the billet with
[HeatConduction](https://mooseframework.inl.gov/source/kernels/HeatConduction.html) and
[HeatConductionTimeDerivative](https://mooseframework.inl.gov/source/kernels/HeatConductionTimeDerivative.html),
taking the eddy current dissipation as a volumetric source through a
[CoupledForce](https://mooseframework.inl.gov/source/kernels/CoupledForce.html). The billet is
stainless steel: $k = 16$ W/m/K, $c_p = 500$ J/kg/K, $\rho = 7900$ kg/m³.

The two flat ends are left adiabatic. The curved side is cooled by the channel, and all of the solid
side of that coupling is set up by one
[CoupledHeatTransfers](https://mooseframework.inl.gov/source/actions/CoupledHeatTransferAction.html)
block, which adds the convective flux boundary condition, the layered side average that reduces the
surface temperature to one value per channel element, and the three transfers that move the wall
temperature out and the fluid temperature and heat transfer coefficient back.

!listing canary/test/tests/multiphysics/induction_heating/induction_heating.i

## The Coolant Channel id=coolant

A single phase water annulus wrapped around the billet, modelled as a one-dimensional
`FlowChannel1Phase` with a `HeatTransferFromExternalAppTemperature1Phase` component, which is the
component built for exactly this arrangement: the solid lives in another app and sends this one its
wall temperature.

| Quantity | Value |
| :- | :- |
| Flow area | $\pi (25^2 - 20^2)\ \text{mm}^2 = 7.069 \times 10^{-4}$ m² |
| Hydraulic diameter | 10 mm |
| Heated perimeter | $2 \pi \times 20\ \text{mm} = 125.7$ mm |
| Mass flow rate | 0.03 kg/s |
| Inlet temperature | 300 K |
| Outlet pressure | 5 bar |
| Heat transfer coefficient | 2000 W/m²/K, held fixed |

The heat transfer coefficient is prescribed rather than taken from a correlation so that the example
does not turn on the choice of closure relation. The pressure is 5 bar so that the saturation
temperature, about 425 K, stays well above anything the channel reaches and the flow stays single
phase.

!listing canary/test/tests/multiphysics/induction_heating/induction_coolant.i

## Coupling id=coupling

Both sub-apps are `TransientMultiApp`s executed at `TIMESTEP_END`. Within a timestep each app sees
the other apps' state from the end of the previous step, so the loop is explicit rather than
Picard-iterated. The evidence that this is good enough here is that halving the timestep moves the
converged temperatures by 2 parts in 100,000: if the lag were doing damage, the answer would depend
on how long the lag was. It is not a free choice in general, though. Early in the transient the
billet gains 20 K in a step, which moves the conductivity by 6% and the dissipation by 3%, and a
problem run closer to thermal runaway, where that feedback is the dominant term rather than a
correction, would want `fixed_point_max_its` on the parent executioner instead.

The coolant channel takes several of its own smaller steps per parent step (`sub_cycling = true`),
because the flow needs a finer timestep than thermal diffusion in the solid does.

The electromagnetic app deliberately declares a very large `dt`. A `TransientMultiApp` that is not
sub-cycling constrains its parent to the smallest timestep among itself and its sub-apps, so a small
`dt` in a sub-app silently overrides whatever the parent asked for. This one has no timescale of its
own to impose, so it names one no parent would ever exceed and leaves the choice to the parent.

One consequence of the explicit loop is that the first thermal step sees no heating at all: the
electromagnetic app has not run yet when it is taken. The dissipation arrives at the end of that
step, and every step after it has a source. Over the full run that startup lag is invisible, but it
is worth knowing when reading the first row or two of the output.

One constraint worth noting on the `joule_heating` transfer: `MultiAppMFEMTolibMeshShapeEvaluationTransfer`
accepts only nodal targets or constant monomials, so the dissipation is computed in a constant $L^2$
space on the MFEM side and received into a `MONOMIAL`, `CONSTANT` auxiliary variable on the libMesh
side. The temperature going the other way is first order Lagrange on both sides.

## Results id=results

Running `induction_heating.i` as it stands takes about five minutes and reaches a steady state in
every quantity by 300 s.

| $t$ (s) | coolant (K) | peak billet (K) | mean surface (K) | power (W) |
| -: | -: | -: | -: | -: |
| 0 | 300.0 | 300.0 | 300.0 | 0 |
| 10 | 304.9 | 412.6 | 339.4 | 5737 |
| 20 | 310.6 | 493.6 | 366.7 | 5982 |
| 40 | 317.1 | 602.5 | 400.7 | 6182 |
| 60 | 321.0 | 670.2 | 421.2 | 6244 |
| 100 | 324.6 | 734.4 | 440.9 | 6274 |
| 200 | 326.4 | 763.9 | 451.1 | 6281 |
| 300 | 326.6 | 765.8 | 451.8 | 6282 |

The two timescales are clearly separated. The power is within 0.6% of its final value by 60 s,
because it is set by the surface temperature, which the coolant pins on a time constant of
$\rho V c_p / h A_{hf} \approx 20$ s. The peak temperature takes five times longer, because it is
set by conduction across the radius, $r^2 \rho c_p / k \approx 100$ s. Stopping at 60 s, as an
earlier version of this example did, leaves the interior 12% short of where it is going while every
integrated quantity looks converged.

### The conductivity feedback id=feedback

At a uniform 300 K the eddy currents dissipate 5473 W. Converged, with the billet between 600 and
770 K, they dissipate 6282 W, a rise of 14.8%. That is the loop closing: the conductivity falls by a
factor of 2.4 over that temperature range, and in the skin-effect regime the dissipated power goes
as $\sigma^{-1/2}$, so less conductive steel absorbs *more* power. The steady state is stable rather
than a runaway because the power grows only as the square root of the resistivity while the heat
removal grows linearly with the temperature difference.

The standalone `induction_em.i` solve is the 5473 W end of that: it runs at the 300 K initial
condition and never receives a temperature, which is what makes it a useful fixed point to test
against on its own.

### Where the heat goes id=profiles

Radially, at the axial midplane:

| $r$ (mm) | 0 | 4 | 8 | 12 | 16 | 19.9 |
| :- | -: | -: | -: | -: | -: | -: |
| $T$ (K) | 764.4 | 765.3 | 765.2 | 752.8 | 713.8 | 616.6 |

The inner half of the billet is isothermal to within a kelvin, and essentially all of the 149 K
radial drop is across the outer 8 mm. That is the shape to expect when the heat is generated in a
skin and has to conduct *outward* through the rest of that skin to reach the coolant: the core has
nowhere to send heat, so it floats up to the temperature of the region that is generating it, while
the generating region itself has to sit above the cooled surface by enough to drive the flux out
through stainless steel at $k = 16$ W/m/K.

The peak is therefore an interior temperature, 765.8 K, sitting 149 K above the surface below it and
314 K above the surface averaged over the whole billet. Watching only `surface_temperature` would
understate the hottest point in the billet by more than 300 K.

Axially, the surface runs from 305 K at the upstream end, through 617 K under the middle of the coil,
to 352 K at the downstream end. The asymmetry is the coolant doing its job: it enters at 300 K and
leaves at 348 K, so the downstream half of the billet is cooled by warmer water and runs hotter. At $z = \pm 50$ mm, symmetric points in both the geometry and the electromagnetic solution,
the surface sits at 416 K upstream and 451 K downstream, a 35 K difference that exists only because
the coolant is being solved for rather than prescribed. A fixed-temperature or fixed-htc boundary
condition on the billet would put those two points at the same temperature.

### Checking it by hand id=hand-check

A closed-form estimate is worth doing before reading a coupled result, and worth being clear about
where it stops working. Taking the converged 6282 W as given:

| Quantity | Hand estimate | Computed |
| :- | -: | -: |
| Coolant rise, $P / \dot{m} c_p$ | 49.1 K | 47.5 K |
| Mean surface, fluid $+\ P / h A_{hf}$ | 450 K | 451.8 K |
| Peak billet | 455 K | 765.8 K |

The first two are close because averaging is exactly what the estimate can legitimately do: total
power over total heated area gives a mean flux, and a mean flux against a mean fluid temperature
gives a mean surface temperature. (The 1.6 K shortfall in the coolant rise is not a modelling
discrepancy but a sampling one: the two postprocessors read the first and last element centroids, at
$z = \pm 90$ mm, and so miss a half element of channel at each end.)

The peak is where the estimate fails, and it fails by 300 K, for two reasons that compound. The coil spans only the middle half of the billet while the whole length is cooled,
so the local flux at the midplane is about twice the length-averaged value and the film drop there
is twice as large. And the interior sits another 149 K above that surface for the conduction reason
above. Neither effect is visible to an argument made only out of totals; both need the solve.
