# Electromagnetic Loads on a Tokamak Vessel During a VDE

During a vertical displacement event (VDE) a tokamak plasma drifts up or down, and the changing
poloidal field induces eddy currents in the vacuum vessel. Crossed with the magnetic field, these
currents load the vessel. Free-boundary equilibrium codes such as
[FREEGSnke](https://github.com/FusionComputingLab/freegsnke) compute the plasma's evolution with an
axisymmetric filament model of the vessel, which cannot represent a three-dimensional structure.
This example computes the eddy currents and Lorentz forces in a 3D vessel with the MFEM backend of
MOOSE, driven by the field of the plasma and coils exported from FREEGSnke.

The coupling is one-way: FREEGSnke's own vessel model determines the plasma's motion, and the MFEM
model computes the vessel's currents and loads in response. Both models here describe the same
axisymmetric vessel, so their vessel currents and forces can be compared.

The example uses only objects in the base framework, so it runs with `moose_test`. It requires an
application built with MFEM support, and is skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

The geometry is the torus of MOOSE's `embedded_concentric_torus.e` test mesh, symmetric about the
$z$ axis:

| Region | Mesh block | Geometry | Material |
| :- | :- | :- | :- |
| Plasma region | `TorusCore` | Major radius 5 m, minor radius 0.5 m | Non-conducting |
| Vessel | `TorusSheath` | Shell between minor radii 0.5 m and 1 m | 316L stainless steel, $\sigma = 1.35 \times 10^6$ S/m |
| Air | `Exterior` | Cube of half-width 7.5 m | Non-conducting |

The source field is axisymmetric: a poloidal flux $\psi_s(R, Z, t)$ per radian, following the FreeGS
conventions $B_R = -\frac{1}{R} \frac{\partial \psi_s}{\partial Z}$ and
$B_Z = \frac{1}{R} \frac{\partial \psi_s}{\partial R}$, together with a toroidal field
$B_\phi(R, Z, t)$. It comprises everything except the conductors in the MFEM model, here the plasma
and the poloidal field coils.

## Formulation id=formulation

The source field enters only through its electric field,

!equation
\vec E_s = -\frac{1}{R} \frac{\partial \psi_s}{\partial t} \hat \phi
= \frac{\partial \psi_s}{\partial t} \frac{(y, -x, 0)}{R^2},

and MFEM solves for the vector potential $\vec A$ of the eddy currents alone:

!equation
\sigma \frac{\partial \vec A}{\partial t} + \nabla \times \left( \nu \nabla \times \vec A \right)
= \sigma \vec E_s \quad \mathrm{in\ the\ vessel}, \qquad
\nabla \times \left( \nu \nabla \times \vec A \right) = 0 \quad \mathrm{elsewhere},

with $\nu = 1 / \mu_0$ and $\vec A \times \hat n = 0$ on the outer boundary. A mass term of
$10^{-6}$ times the curl-curl term, for fields varying over 1 m, makes the operator definite in the
non-conducting regions. $\vec A$ is a lowest-order Nédélec field, solved by preconditioned conjugate
gradients with Hypre's AMS, and time is integrated by implicit Euler. The eddy current density, the
total field and the Lorentz force density are then

!equation
\vec J = \sigma \left( \vec E_s - \frac{\partial \vec A}{\partial t} \right), \qquad
\vec B = \vec B_s + \nabla \times \vec A, \qquad
\vec f = \vec J \times \vec B.

These are evaluated on an `MFEMDomainSubMesh` of the vessel, so the source data only has to cover the
conductors. The time derivative and the curl of $\vec A$ are transferred to it with
`MFEMSubMeshTransfer`, $\vec f$ is projected onto an $[L^2]^3$ variable with `MFEMCrossProductAux`, and
`MFEMVectorFEInnerProductIntegralPostprocessor` integrates

!equation
I_\phi = \int_V \vec J \cdot \frac{\hat \phi}{2 \pi R} \, dV, \qquad
F_R = \int_V \vec f \cdot \hat R \, dV, \qquad
F_Z = \int_V \vec f \cdot \hat z \, dV.

$I_\phi$ is the toroidally averaged current of the vessel, which for an axisymmetric current is the
current through a poloidal cross-section.

### Reading the source field id=source

The source fields are written on a uniform $(R, Z)$ grid as `PiecewiseMultilinear` data over
$(R, Z, t)$, with the quantities $\partial \psi_s / \partial t$, $B_R$, $B_Z$ and $B_\phi$. A generated
MOOSE input converts them to a 2D ExodusII time series. In the 3D model a `SolutionUserObject` reads
that series and interpolates it in time, `Axisymmetric2D3DSolutionFunction` maps it from the $(R, Z)$
half-plane onto the 3D mesh, and `ParsedVectorFunction` assembles $\vec E_s$ and $\vec B_s$ from it.
All three are MOOSE functions that the MFEM backend uses as coefficients.

!listing canary/test/tests/electromagnetics/torus_vessel/torus_vessel.i

## Verification id=verification

The regression test drives the vessel with a uniform vertical field ramped at
$\dot B = 0.1$ T/s, whose flux $\psi_s = B_Z R^2 / 2$ satisfies the vacuum Grad-Shafranov equation. Once
the ramp has run for several L/R times of the vessel, the eddy currents stop changing and
$\vec J = \sigma \vec E_s$. After the 10 s of the regression test, with 1 s steps, $I_\phi$ is within
$9 \times 10^{-4}$ of the current that $\sigma \vec E_s$ carries on the same mesh, the remainder of the
initial transient; run to 20 s with 0.5 s steps, it settles to within $6 \times 10^{-6}$. The vertical
force is $7 \times 10^{-6}$ of the radial force, as the up-down symmetry of the problem requires.

$\partial \vec A / \partial t$ in the vessel settles at 0.3% of $\vec E_s$ rather than at zero. On the
faceted surface of the tetrahedral mesh the discrete $\vec E_s$ has a small normal component, and the
gradient part of $\partial \vec A / \partial t$ acts as the electrostatic potential that keeps
$\vec J \cdot \hat n = 0$ there.

## Coupling to FREEGSnke id=freegsnke

The scripts in `vde/` build a FREEGSnke machine whose vessel is the `TorusSheath` shell, split into
poloidal sectors and radial layers of passive structures, with the plasma inside `TorusCore` and a
circular limiter of radius 0.45 m. Up-down symmetric coil pairs hold a 150 kA plasma, limited at the
top and bottom of the limiter with an elongation of 1.38, and a radial field coil pair is wired
antisymmetrically to displace it.

The elongated plasma is vertically unstable, with a growth rate of 3.9 s$^{-1}$ in the linearised
FREEGSnke model. A constant voltage on the radial field coils starts the plasma moving, and over
0.68 s its magnetic axis rises 14 cm. FREEGSnke's time step is 0.02 s.

At each step `freegsnke_source.py` evaluates the source field of the coils and plasma on a grid
covering the vessel, excluding the vessel's own passive structures so their field is not counted
twice. Inside the FREEGSnke grid the plasma's flux is FREEGSnke's own solution; outside it, the plasma
current is summed with Green's functions. $\partial \psi_s / \partial t$ is the backward difference
over each step, matching an implicit Euler step of the same length, so the MFEM time step must equal
FREEGSnke's.

`run_vde.py` also records FREEGSnke's own vessel current and vertical force for comparison. The
force is the sum over the vessel's filaments of their currents times the radial field of the plasma
and coils, since the forces the vessel exerts on itself sum to zero.

## Results id=results

### The air region id=air-region

The air region of `embedded_concentric_torus.e` ends 1.5 m from the vessel. With the outer boundary
so close, the vessel's inductance is lower than in free space, where FREEGSnke treats it. Its slowest
L/R time is about 0.97 s in MFEM against 1.18 s in FREEGSnke. For times short compared with L/R the eddy
current scales inversely with the inductance, so MFEM's vessel current comes out about 20% higher than
FREEGSnke's.

`vde/make_mesh.py` meshes the same torus inside a larger air cube with Gmsh, meshing the torus the same
way whatever the cube's size. Growing the cube's half-width from 7.5 m to 15 m and 30 m reduces the
largest difference in vessel current over the VDE from 19.5% to 2.4% and 0.3%.

!media media/torus_vessel_air_region.png
       style=width:70%;margin-left:auto;margin-right:auto;display:block;
       caption=The VDE and the vessel's response, from FREEGSnke and from MFEM with three sizes of
               air region. The vessel current converges to FREEGSnke's as the air region grows; the
               vertical force does not.

### Resolving the vertical force id=resolution

The vertical force barely changes with the air region. It comes mainly from the up-down antisymmetric
part of the eddy currents, which opposes the plasma's motion, and depends on how that current is
distributed through the 0.5 m thickness of the vessel. Both models must resolve it. Refining each
moves its force towards the other's:

| FREEGSnke vessel (sectors x layers) | MFEM element size in the vessel | Largest vertical force difference |
| :- | :- | :- |
| 24 x 2 | 0.25 m | 30% |
| 24 x 2 | 0.15 m | 15% |
| 48 x 4 | 0.25 m | 19% |
| 48 x 4 | 0.15 m | 6.4% |

All four use the 30 m air cube. Refining the FREEGSnke vessel also changes the VDE itself slightly: the
growth rate falls to 3.8 s$^{-1}$ and the plasma moves 11 cm by 0.66 s.

!media media/torus_vessel_vde.png
       style=width:70%;margin-left:auto;margin-right:auto;display:block;
       caption=The VDE and the vessel's response with the 48 x 4 FREEGSnke vessel and MFEM elements
               of 0.15 m in a 30 m air cube. The vessel currents agree to 0.35% and the vertical
               forces to 6.4%.

## Files id=files

| File | Role |
| :- | :- |
| `torus_vessel.i` | The MFEM model of the vessel. This is the input to run. |
| `ramp_source.i` | Converts the ramp source of the regression test to ExodusII. Generated by `generate_ramp_source.py`. |
| `ramp_source_*.txt`, `ramp_source_times.csv` | The ramp source fields and their times. |
| `generate_ramp_source.py` | Writes the ramp source. |
| `gridded_fields.py` | Writes source fields and their conversion input. Used by the ramp and the FREEGSnke export. |
| `embedded_concentric_torus.e` | The mesh of the regression test, copied from MOOSE's `test/tests/mfem/mesh`. |
| `vde/torus_machine.py` | The FREEGSnke machine. |
| `vde/run_vde.py` | Runs the VDE in FREEGSnke and exports its source field. |
| `vde/freegsnke_source.py` | Evaluates the source field of a FREEGSnke equilibrium. |
| `vde/make_mesh.py` | Meshes the torus in a larger air region with Gmsh. |
| `vde/compare.py` | Compares FREEGSnke's and MFEM's vessel response and draws the figures above. |
| `tests` | Test harness specification. |
| `gold/` | Reference output. |

All of them live in `test/tests/electromagnetics/torus_vessel`.

## Running id=running

From the repository root, the regression test runs through the test harness in about 30 seconds:

```bash
./run_tests --re torus_vessel
```

The VDE needs FREEGSnke and Gmsh, which the test harness does not use. With `OUT` an output
directory, `EXAMPLE` this example's directory, and `<app>` a MOOSE application built with MFEM
support:

```bash
python $EXAMPLE/vde/run_vde.py $OUT --n-poloidal 48 --n-radial 4
python $EXAMPLE/vde/make_mesh.py 30 $OUT/torus.msh 0.15
cd $OUT
<app> -i vde_source.i
<app> -i $EXAMPLE/torus_vessel.i Mesh/file=$OUT/torus.msh \
  UserObjects/source/mesh=$OUT/vde_source_out.e \
  Executioner/dt=0.02 Executioner/end_time=0.66 \
  Outputs/active=csv Outputs/csv/file_base=torus_vde
python $EXAMPLE/vde/compare.py freegsnke_vde.csv torus_vde.csv --plot comparison.png
```

The FREEGSnke run takes about 5 minutes and the MFEM run about 9 minutes on one process. The end
time must be the last time in `vde_source_times.csv`. FREEGSnke stops once the magnetic axis has moved
0.12 m (`--max-axis-shift`): the plasma touches the limiter at top and bottom, so as it moves it is
scraped off, and FREEGSnke stalls once the plasma region vanishes. Adding `vessel` back to
`Outputs/active` writes the vessel's fields for ParaView.
