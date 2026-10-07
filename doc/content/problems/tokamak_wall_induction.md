# Eddy Currents Induced in a Tokamak Wall

Whenever the current in a tokamak plasma changes, so does the magnetic flux it links, and the change induces
currents in the conducting structures around it. Those eddy currents oppose the change, shield the plasma's
field from the outside, and then decay resistively. They matter most when the plasma is lost: in a disruption
the plasma current falls to zero within milliseconds, and the wall carries a comparable current for tens of
milliseconds afterwards, dissipating the plasma's magnetic energy and feeling the forces that come with it.

This example is a crude model of that: a torus of plasma current, inside a toroidal stainless steel wall, in
free space. The plasma current density is an input: any vector function $\vec{J}(x, y, z, t)$ of position and
time, with no assumption that it is time-harmonic, separable or axisymmetric. The low-frequency
(magnetoquasistatic) eddy current problem is solved in the time domain in 3-D, and reports the current induced
in the wall, the power it dissipates and the field at two magnetic probes. It is verified against an independent
axisymmetric filament (coupled circuit) model of the same device.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent. It also relies on `MFEMTransient` evaluating
time-dependent coefficients, and stepping its MultiApps, at the end-of-step time, as MOOSE's other transient
executioners do.

## Problem Specification id=specification

!media media/tokamak_wall_geometry.png
       id=geometry-figure
       caption=The poloidal cross-section of the device, which is revolved about the $z$ axis. The plasma is
               shaded by its current density. The dashed circle is where the plasma ends up in the
               vertical displacement event below; the crosses are the two magnetic probes. The air cylinder
               around the device, of radius and half-height 6 m, is not shown.
       style=width:70%;margin-left:auto;margin-right:auto;

| Quantity | Value |
| :- | :- |
| Major radius $R_0$ of the plasma and the wall | 1.0 m |
| Plasma minor radius $a$ | 0.3 m |
| Wall minor radius (inner surface) $b$ | 0.4 m |
| Wall thickness $d$ | 0.03 m |
| Wall conductivity $\sigma$ | $1.351 \times 10^{6}$ S/m: 316L stainless steel, resistivity 74 $\mu\Omega$ cm |
| Relative permeability | 1, everywhere |
| Free space | a cylinder of radius 6 m and height 12 m, centred on the device |
| Probes | $B_z$ at the centre, $R = 0$, $z = 0$, and on the outboard midplane at $R = 1.35$ m, $z = 0$ |

The wall is a complete toroidal shell, with no toroidal break, so it can carry a net toroidal current.

### The Plasma Current id=plasma-current

By default the plasma current $I_p(t)$ is ramped up linearly to 1 MA over 10 ms, held for 20 ms, and then lost
in a 2 ms linear current quench, as in a disruption. It flows toroidally, with a parabolic profile peaked on the
magnetic axis:

!equation
\vec{J} = \frac{2 I_p(t)}{\pi a^2} \left( 1 - \frac{\rho^2}{a^2} \right) \hat{\phi}, \quad \rho < a,
\qquad \rho^2 = (R - R_0)^2 + (z - z_\text{axis}(t))^2,

with $\hat\phi = (-y, x, 0) / R$ the toroidal direction. The magnetic axis stays on the midplane,
$z_\text{axis} = 0$, unless a vertical displacement is asked for; see
[#vertical-displacement-event].

All of this lives in one short input, `plasma_current_density.i`, which is the file to edit to change the
source. It defines the vector function `plasma_current_density` as a `ParsedVectorFunction` of $x$, $y$, $z$
and $t$, built from a `PiecewiseLinear` waveform for $I_p(t)$; that waveform can equally be read from a CSV
file of a measured discharge. Only the values of $\vec{J}$ in the source blocks (`source_blocks`, the
plasma by default) are used.

!listing canary/test/tests/electromagnetics/tokamak_wall/plasma_current_density.i
         start=major_radius

$\vec J$ should be divergence free, as any physical current density is: the magnetoquasistatic model has no
room for a current that piles up charge. The plasma current the solve actually receives is reported as the
postprocessor `plasma_current`, which shows when the two differ.

## Mesh id=mesh

The mesh is built by `tokamak_wall_mesh.py` with gmsh, by revolving a triangulation of the poloidal half-plane
about the $z$ axis in 64 toroidal layers and splitting the prisms into tetrahedra. Revolving gives elements that
are long in the toroidal direction and short across the 3 cm wall, which an isotropic 3-D mesh generator would
not, and keeps the cost down: the default mesh has 464096 tetrahedra, of which 49152 are in the wall, with two
through its thickness. The column of air around the axis, $R < 0.4$ m, is meshed with unstructured tetrahedra
instead, since revolved elements touching the axis would be degenerate. The result is a full 3-D mesh of the
whole torus, with no symmetry assumed, so non-axisymmetric sources need nothing changed.

The blocks are `plasma` ($\rho < a$), `vacuum` (between the plasma and the wall), `wall`, and `air`
(everything else); the outer surface of the air cylinder is the sideset `far_field`. Every dimension and element
size can be changed on the script's command line, and `--refine` scales all the element sizes and layer counts
together for convergence studies. The script writes the mesh as gmsh's `.msh` and converts it to Exodus with
MOOSE's `--mesh-only`, since MFEM picks up block and sideset names from Exodus files only.

## Formulation id=formulation

### The Vector Potential Solve id=vector-potential

The magnetic field is written as $\vec{B} = \nabla \times \vec{A}$ with the vector potential $\vec{A}$ in
$H(\text{curl})$, in the temporal gauge, where the electric field is $\vec{E} = -\partial\vec{A}/\partial t$.
Ampere's law with the eddy current $\sigma\vec{E}$ and the applied plasma current $\vec{J}$ gives

!equation
\left( \sigma \frac{\partial \vec{A}}{\partial t}, \vec{A}' \right)
+ \left( \nu \nabla \times \vec{A}, \nabla \times \vec{A}' \right)
= \left( \vec{J}, \vec{A}' \right),

with $\nu = 1/\mu_0$ everywhere, integrated in time by implicit Euler. Nothing about the form of $\vec{J}(t)$ is
assumed. The kernels are `MFEMTimeDerivativeVectorFEMassKernel` and `MFEMCurlCurlKernel`, with the source applied
through the mixed kernels described below; the linear system at each step is solved by conjugate gradients
preconditioned by Hypre's auxiliary-space Maxwell solver (AMS).

Physically the conductivity outside the wall is zero, which would leave every gradient supported there in the
null space of the operator. A conductivity of $10^{-2}$ S/m, eight orders of magnitude below the wall's, is given to the
plasma, vacuum and air instead to regularise it. Its effect is linear in it, and negligible: raising it a
hundredfold, to 1 S/m, changes the wall current by $5 \times 10^{-5}$ of its peak, and lowering it a hundredfold
by $5 \times 10^{-7}$, while the solver takes the same 27 iterations a step.

The far-field surface carries $\hat{n} \times \vec{A} = 0$, so that $\vec{B} \cdot \hat{n} = 0$ there: no flux
leaves the air cylinder. The far field of the device is that of a current loop, which falls off as
the cube of the distance, and at 6 m the boundary costs about 0.2% in the wall current; see
[#convergence].

### Making the Source Divergence Free at Every Step id=source

The right hand side has to be orthogonal to every discrete gradient outside the wall. Those are held only by
the tiny regularising conductivity, and a source that drives them produces potentials, and with them currents
in the wall, in inverse proportion to it. $\vec{J}$ itself is divergence free, but its
projection onto the mesh is not quite: the meshed edge of the source region is faceted, and a little current
leaks across it.

The sub-app `plasma_current_source.i` removes that leak. At every time step, at the time the step is being
taken to, it projects the user's $\vec J$ onto $H(\text{curl})$ as $\vec{J}_h$, and solves

!equation
\left( \nabla\phi, \nabla\psi \right) = \left( \vec{J}_h, \nabla\psi \right)_\text{source}
\quad \forall \psi,

for a correction potential $\phi$ over the whole mesh, pinned on `far_field`. The current the main app uses,
$\vec{J}_h$ restricted to the source blocks minus $\nabla\phi$, then satisfies
$(\vec{J}_h - \nabla\phi, \nabla\psi) = 0$ for every $\psi$ by construction. This is a Helmholtz projection: the
smallest correction that makes the source discretely compatible. The two pieces are copied to the main app
with `MultiAppMFEMCopyTransfer` and applied by `MFEMMixedVectorMassKernel` over the source blocks and
`MFEMMixedVectorGradientKernel` over the whole mesh. The correction is a gradient, so it carries no net toroidal
current.

Because $\vec{J}$ can be any function of time, this happens every step: the sub-app is a `TransientMultiApp`
stepped with the main app, executed at `TIMESTEP_BEGIN`. That needs both apps to evaluate time-dependent
coefficients at the end-of-step time $t_{n+1}$ before the solve, as implicit Euler requires.

!listing canary/test/tests/electromagnetics/tokamak_wall/plasma_current_source.i
         start=[Kernels]
         end=[Solvers]

### Reported Quantities id=reported-quantities

The toroidal current in the wall is the integral of $\sigma\vec{E} \cdot \hat\phi / (2\pi R)$ over the wall's
volume, which for an axisymmetric current is exactly the current through a poloidal cross-section, and in
general is that current averaged around the torus. Its vertical moment, $\int J_\phi z \, dA$, is the same
integral weighted by $z$; it vanishes for currents symmetric about the midplane. Both weights are projected
onto $H(\text{curl})$ once and integrated against $\vec{E} = -\partial\vec{A}/\partial t$ with
`MFEMVectorFEInnerProductIntegralPostprocessor`, which also gives the ohmic power $\int \sigma|\vec{E}|^2 dV$.
The probe fields are point samples of $\vec{B}$. The plasma current the solve received is reported the same way,
from $\vec{J}_h$.

!listing canary/test/tests/electromagnetics/tokamak_wall/tokamak_wall.i

## Verification id=verification

### The Filament Model id=filament-model

When the plasma and the wall are both axisymmetric, the induced current is purely toroidal, there is no
electrostatic potential, and the problem reduces to circuits. `filament_model.py` cuts the wall's poloidal
cross-section into small square cells, each of which is a closed toroidal loop $k$ carrying a current $I_k$, and
cuts the plasma into filaments carrying its current profile. Faraday's law around each wall loop is

!equation
\sum_j L_{kj} \frac{dI_j}{dt} + R_k I_k = -\frac{d\Psi_k}{dt},

with $R_k$ the loop's toroidal resistance, $L$ the inductance matrix of the wall loops, and $\Psi_k(t)$ the
plasma's flux through loop $k$. The mutual inductance of two coaxial circular loops is Maxwell's formula in
complete elliptic integrals; the self inductance of a cell is that of a thin ring through the geometric mean
distance of its rectangular cross-section, and neighbouring cells are coupled through Gauss points in both so
that their finite size is accounted for. The system is integrated exactly in time, through the modes of
$R\vec{v} = \lambda L\vec{v}$, for a plasma flux that is piecewise linear in time; it can also be integrated by
implicit Euler at the finite element time step, which isolates the finite element model's spatial error.

The model shares nothing with the finite element solution: no mesh, no vector potential and no far-field
boundary. It converges quickly as the cells are refined:

| Wall cells | Slowest wall time constant (ms) | Peak wall current (kA) |
| :- | -: | -: |
| 1 $\times$ 87 | 22.7535 | 804.91 |
| 2 $\times$ 174 | 22.7531 | 804.72 |
| 4 $\times$ 348 | 22.7530 | 804.68 |
| 8 $\times$ 695 | 22.7530 | 804.68 |

The 4 by 348 cell model is used as the reference. Its slowest mode, the wall's $L/R$ time, is 22.75 ms; the
next are 9.0 and 8.5 ms. The wall power during the current quench is the quantity most sensitive to the cells,
since it depends on how the current is distributed through the wall's thickness; going from 4 to 8 cells
through the thickness raises its peak by 0.2%.

### Results id=results

!media media/tokamak_wall_comparison.png
       id=comparison-figure
       caption=The default disruption scenario: `tokamak_wall.i` at its default resolution and time step against
               the filament model, integrated exactly in time.
       style=width:100%;margin-left:auto;margin-right:auto;

While the plasma current is ramped up over 10 ms, less than half the wall's $L/R$ time, the wall carries an
opposing current that reaches 800 kA, 80% of the plasma current, and shields most of its field from the centre
of the machine. Over the flat top the wall current decays, letting the field in. The 2 ms current quench then
swings the wall current by almost a megaamp, to $+663$ kA: the wall takes over two thirds of the vanished plasma
current, so that the flux through the centre barely changes during the quench and then decays with the wall,
over tens of milliseconds. The wall dissipates 37 MW at the end of the ramp, 25 MW at the end of the quench,
and 743 kJ in all.

The finite element solution follows the reference throughout. Its largest differences over the whole run,
relative to the peak of each quantity, are below. The time step is the same 0.25 ms in the finite element model
and the implicit Euler reference, so the second column is the spatial error alone.

| Quantity | Peak | Largest difference from exact | from implicit Euler |
| :- | -: | -: | -: |
| Wall current | $-799.5$ kA | 0.68% | 0.40% |
| $B_z$ at the centre | 434.9 mT | 0.87% | 0.27% |
| $B_z$ at the outboard probe | $-459.1$ mT | 1.43% | 1.34% |
| Ohmic power in the wall | 36.7 MW | 7.5% | 8.7% |
| Energy dissipated in the wall | 743.5 kJ | | +0.42% |

The plasma current the solve received is within 0.03% of the specified one. The wall current's vertical moment,
zero by symmetry, is at most 0.6 kA m, 0.2% of the wall current times the wall's minor radius: the
unstructured parts of the mesh are not quite symmetric about the midplane.

The ohmic power is the exception. Away from the quench it agrees as well as everything else, but during the
quench the finite element model dissipates 1 to 3 MW more than the reference, while carrying the same net
current. The excess is in the distribution of the current, not its total, and it does not shrink steadily with
the mesh (see below); lowest order Nedelec elements represent the rapidly changing, strongly rotational
electric field in the long toroidal elements poorly. Second order elements (`fe_order=SECOND flux_order=FIRST`)
halve it on the coarse mesh. The dissipated energy, which is the integral of the power, is not affected to any
degree that matters: it agrees to 0.4%.

### Convergence id=convergence

`tokamak_wall_mesh.py --refine` scales every element size and the number of toroidal and wall layers together.
The runs below use a 1 ms time step and are compared with the filament model integrated by implicit Euler at the
same step, so they show the spatial error alone; the last row repeats the default mesh with the air cylinder
doubled to 12 m.

| Mesh | Tetrahedra | Toroidal layers | Wall layers | Wall current | $B_z$ centre | $B_z$ outboard | Energy | Peak power |
| :- | -: | -: | -: | -: | -: | -: | -: | -: |
| `--refine 0.5` | 70035 | 32 | 1 | 0.51% | 2.07% | 5.02% | $-0.19$% | 5.8% |
| default | 464096 | 64 | 2 | 0.40% | 0.22% | 1.33% | +0.37% | 7.2% |
| `--refine 1.5` | 1522992 | 96 | 3 | 0.33% | 0.32% | 0.33% | +0.18% | 5.4% |
| default, 12 m air | 685380 | 64 | 2 | 0.18% | 0.29% | 1.41% | +0.61% | 6.9% |

The outboard probe sits 5 cm from the edge of the plasma, where $B_z$ changes by some 50 mT across an element,
and it is a point value of the lowest order $\vec B = \nabla \times \vec A$, which is constant in each element; its
error falls from 5% to 1.3% and 0.3% as the mesh is refined. The wall current is within half a percent already on
the coarsest mesh, and what remains of its error on the default one is mostly the far-field boundary: doubling the
air cylinder takes it from 0.40% to 0.18%.

With the default mesh, halving the time step halves the difference from the exact solution, as implicit Euler
should:

| Time step | Wall current | $B_z$ centre |
| -: | -: | -: |
| 1 ms | 2.75% | 3.44% |
| 0.5 ms | 1.39% | 1.75% |
| 0.25 ms | 0.68% | 0.87% |

The plasma current the solve receives is itself a discretisation of the specified one. Projecting the parabolic
profile onto the mesh misses some of its curvature in the poloidal plane, which loses current as the square of
the element size there, and samples it on the chords of the toroidal layers, which lie inside the circles they
approximate and gain current as the square of the layer angle. The ratio of the received to the specified
current, for a few combinations of the two:

| Plasma element size | 32 layers | 64 layers | 128 layers |
| -: | -: | -: | -: |
| 6 cm | 1.0057 | 0.9948 | 0.9919 |
| 4 cm | 1.0105 | 0.9997 | |
| 3 cm | 1.0122 | | |

The default, 64 layers and 4 cm, is within 0.03%; the two errors partly cancel there, but each is below 0.5%.

## Vertical Displacement Event id=vertical-displacement-event

Vertical displacement events, in which an unstable plasma drifts up or down into the wall before or while it is
quenched, are among the worst cases for the forces on a tokamak's vessel. `plasma_current_density.i` models one
crudely: with `vertical_displacement` set, the plasma rises by that distance at a steady rate over the flat top,
from 20 to 30 ms, and is quenched where it ends up. The current density is then not a fixed profile scaled by a
function of time. Nothing in the solve changes, since the sub-app re-projects and corrects $\vec J$ at every step;
the source blocks just have to take in the vacuum the plasma moves into:

```bash
<app> -i tokamak_wall.i source_blocks='plasma vacuum' plasma_current_source:vertical_displacement=0.08
```

The filament model handles the moving plasma as well, with $\Psi_k(t) = I_p(t) M_k(z_\text{axis}(t))$.

!media media/tokamak_wall_vde.png
       id=vde-figure
       caption=A vertical displacement event of 8 cm: the plasma rises from 20 to 30 ms and is then quenched.
       style=width:100%;margin-left:auto;margin-right:auto;

The net wall current hardly changes, but as the plasma rises the wall carries currents that oppose the motion,
running against the plasma current in the upper half of the wall and with it in the lower half, so that their
vertical moment grows to $-50$ kA m by the start of the quench. The quench reverses it, along with the net
current. The finite element solution agrees with the reference as well as for the stationary plasma:

| Quantity | Peak | Largest difference from exact | from implicit Euler |
| :- | -: | -: | -: |
| Wall current | $-800.0$ kA | 0.58% | 0.48% |
| Vertical moment of the wall current | $-49.6$ kA m | 2.11% | 1.18% |
| $B_z$ at the centre | 431.7 mT | 0.78% | 0.36% |
| $B_z$ at the outboard probe | $-459.3$ mT | 2.12% | 2.14% |
| Energy dissipated in the wall | 759.3 kJ | | +0.91% |

The outboard probe is less accurate than with the stationary plasma, because the displaced plasma brings its
edge, and the steepest field, closer to it. The received plasma current rises by 0.4% as the plasma moves
into the finer elements of the vacuum.

## Files id=files

| File | Role |
| :- | :- |
| `tokamak_wall.i` | The eddy current solve. This is the input to run. |
| `plasma_current_source.i` | Sub-app: projects $\vec J$ at each step and makes it divergence free. |
| `plasma_current_density.i` | The plasma current density, included by the sub-app. Edit this to change the source. |
| `tokamak_wall_mesh.py` | Builds `tokamak_wall.e` with gmsh. |
| `tokamak_wall.e` | The mesh. |
| `filament_model.py` | The axisymmetric filament model used as the reference. |
| `plot_tokamak_wall.py` | Compares a run with the filament model, and draws the figures on this page. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

All of them live in `test/tests/electromagnetics/tokamak_wall`.

## Running id=running

```bash
mpirun -n 16 <path to a MOOSE app built with MFEM support> -i tokamak_wall.i
```

or, from the repository root, through the test harness:

```bash
./run_tests --re Tokamak_Wall -p 16
```

The full run, 240 steps of 0.25 ms, takes about 10 minutes on 16 processes and peaks at around 1 GB per
process; the verification test that runs it is marked `heavy`. The two regression tests, of the default and the
vertical displacement event, take 20 steps of 2 ms each and under a minute on 8 processes.

The postprocessors are written to `OutputData/TokamakWall.csv` after every step, and the fields - $\vec A$,
$\vec B$, $\vec E$ and the toroidal current density - to `OutputData/TokamakWall` for ParaView every 2 ms. Every parameter
at the top of the inputs can be changed from the command line, and the sub-app's with the prefix
`plasma_current_source:`. For example, a coarser time step and a 2 MA plasma:

```bash
<app> -i tokamak_wall.i dt=1e-3 plasma_current_source:peak_plasma_current=2e6
```

To compare a run with the filament model:

```bash
python3 plot_tokamak_wall.py OutputData/TokamakWall.csv
```

## References id=references

- F. W. Grover, *Inductance Calculations*, Dover (1946) - mutual inductance of coaxial circles, and geometric
  mean distances.
- J. C. Maxwell, *A Treatise on Electricity and Magnetism*, vol. 2, art. 701 (1873) - the mutual inductance of
  two coaxial circles in elliptic integrals.
- R. Albanese and G. Rubinacci, *Finite element methods for the solution of 3D eddy current problems*,
  Advances in Imaging and Electron Physics, 102, 1 (1998) - the $\vec{A}$ formulations of the eddy current
  problem.
