# Cahn-Hilliard Spinodal Decomposition with MFEM

The Cahn-Hilliard equation describes a conserved order parameter $c$ - a concentration, say - in a mixture
that lowers its free energy by separating into two phases. Started from small noise about a uniform
state inside the spinodal region, the mixture decomposes into an interconnected pattern of the two phases,
which then coarsens as the interfaces between them shorten.

The MOOSE phase field module solves it with dedicated kernels, and its example
[Parsed_SplitCH.i](https://github.com/idaholab/moose/blob/master/modules/phase_field/examples/cahn-hilliard/Parsed_SplitCH.i)
is the problem set up here. This example solves the same problem with the MFEM backend of MOOSE, using only
objects in the base framework, and compares the result with the phase field module's.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

The free energy of the mixture is

!equation
F = \int_\Omega f(c) + \frac{\kappa}{2} \left| \nabla c \right|^2 \, d\Omega,
\qquad f(c) = W \left( c^2 - 1 \right)^2,

a double well with minima at $c = \pm 1$ and a gradient energy penalising interfaces between them. The
order parameter evolves down the gradient of $F$ while conserving its total,

!equation
\frac{\partial c}{\partial t} = \nabla \cdot \left( M \nabla \mu \right),
\qquad \mu = \frac{\delta F}{\delta c} = f'(c) - \kappa \nabla^2 c.

All quantities are dimensionless, with the values of `Parsed_SplitCH.i`.

| Quantity | Value |
| :- | :- |
| Barrier height $W$ | 1/4 |
| Mobility $M$ | 1 |
| Gradient energy coefficient $\kappa$ | 0.5 |
| Domain | $60 \times 60$ square |
| Mesh | $150 \times 150$ bilinear quadrilaterals |
| Initial condition | $c$ uniformly distributed in $[-0.1, 0.1]$ from node to node |
| Time step | $\Delta t = 2$, to $t = 20$ |

With these values the equilibrium interface is $c = \tanh(x / \delta)$ with
$\delta = \sqrt{\kappa / 2W} = 1$, which the mesh spacing of 0.4 resolves with a few elements, and the
fastest growing wavelength of the initial instability is $2 \pi \sqrt{2 \kappa / \lvert f''(0) \rvert} = 2 \pi$.

## Formulation id=formulation

The phase field module splits the fourth order equation into two second order ones, one for $c$ and one
for the whole chemical potential $\mu$, with $f'(c)$ in the equation for $\mu$:

!equation
\left( \frac{\partial c}{\partial t}, v \right) + \left( M \nabla \mu, \nabla v \right) = 0, \qquad
\left( \mu, q \right) - \left( f'(c), q \right) - \left( \kappa \nabla c, \nabla q \right) = 0.

That term cannot be written with the MFEM objects in the base framework. Kernels whose coefficients depend on
the solution are only re-evaluated as Newton iterates when they are non-linear form integrators, and the
only scalar one, `MFEMNLDiffusionKernel`, is a diffusion term that may act only on its own variable.
$(f'(c), q)$ is neither.

Only the gradient energy part of the chemical potential, $w = -\kappa \nabla^2 c$, is split off here instead.
The bulk part is kept in the equation for $c$ where, because $\nabla f'(c) = f''(c) \nabla c$, it becomes a
non-linear diffusion of $c$:

!equation
\left( \frac{\partial c}{\partial t}, v \right)
+ \left( M f''(c) \nabla c, \nabla v \right)
+ \left( M \nabla w, \nabla v \right) = 0, \qquad
\left( w, q \right) - \left( \kappa \nabla c, \nabla q \right) = 0.

The first is exactly what `MFEMNLDiffusionKernel` provides, with $k(c) = M f''(c) = 4 M W (3 c^2 - 1)$, and it
is given $\mathrm{d}k/\mathrm{d}c = 24 M W c$ as well so that Newton's method has the exact Jacobian. $k$ is
negative where $|c| < 1/\sqrt{3}$, inside the spinodal, and the uphill diffusion it describes there is what
drives the separation. The remaining terms are bilinear: `MFEMTimeDerivativeMassKernel`, `MFEMMassKernel`, and
`MFEMMixedGradGradKernel` coupling $c$ and $w$ each way. $c$ and $w$ are both first order $H^1$ fields, and
the coupled system is solved monolithically by `MFEMNewtonNonlinearSolver` with an LU factorisation from
MUMPS at each iteration, since its Jacobian is not symmetric. `MFEMTransient` integrates in time by implicit
Euler.

In the continuum the two splittings are the same equation. Discretely they differ: the phase field module's
$\mu$ is the projection of $f'(c) - \kappa \nabla^2 c$ onto the finite element space, where here
$\nabla f'(c)$ is taken exactly. They therefore agree to within the discretisation error, which the
[comparison](#comparison) below measures. The full chemical potential is recovered after each step as the
auxiliary variable `chemical_potential` $= f'(c) + w$, which is the variable the phase field module solves
for.

### Boundary and Initial Conditions id=boundary-conditions

All four sides are left natural, which is zero flux, $M \nabla \mu \cdot \vec{n} = 0$, together with
$\kappa \nabla c \cdot \vec{n} = 0$, a 90 degree contact angle. Testing the equation for $c$ with $v = 1$ then
shows that $\int_\Omega c \, d\Omega$ is conserved exactly. `Parsed_SplitCH.i` is periodic instead;
MFEM-MOOSE cannot yet build a periodic mesh.

`Parsed_SplitCH.i` starts from `RandomIC`, whose values depend on the random number generator and the
partitioning. Here the noise is a hash of the nodal coordinates,

!equation
c_0(x, y) = 0.1 \left( 2 \, \mathrm{fract}\!\left( 43758.5453 \sin(12.9898 x + 78.233 y) \right) - 1 \right),

which is just as uncorrelated from node to node, the same in serial and in parallel, and can be given to the
phase field module too, so that both start from identical nodal values.

### Free Energy and Execution Order id=execution-order

The total free energy is reported in its two parts. The gradient energy is the integral of
$\frac{\kappa}{2} \nabla c \cdot \nabla c$, with $\nabla c$ held exactly in the lowest order $H(\text{curl})$
space by `MFEMGradAux`. The bulk energy is the integral of $f(c)$, taken with the same
`MFEMVectorFEInnerProductIntegralPostprocessor` applied to a unit vector field.

These objects read $c$ through parsed coefficients, which MFEM-MOOSE does not count as dependencies when it
orders executed objects. At the end of a time step that is harmless, because $c$ comes from the solve, but at
the start the initial condition on $c$ is itself one of the objects being ordered, and they can run before it.
Everything is therefore computed at the end of each time step only, and the CSV output has no row for
$t = 0$.

!listing canary/test/tests/phase_field/cahn_hilliard/cahn_hilliard.i

## Files id=files

| File | Role |
| :- | :- |
| `cahn_hilliard.i` | The MFEM solve. This is the input to run. |
| `cahn_hilliard_phase_field.i` | The same problem with the phase field module, for the comparison. Not run by the test harness. |
| `plot_cahn_hilliard.py` | Draws the figures below and prints the comparison. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

All of them live in `test/tests/phase_field/cahn_hilliard`.

## Running id=running

```bash
<path to a MOOSE app built with MFEM support> -i cahn_hilliard.i
```

or, from the repository root, through the test harness:

```bash
./run_tests --re Cahn_Hilliard
```

It takes about 35 seconds on one process. The postprocessors - the mean of $c$ and the total free energy -
are written to `OutputData/CahnHilliard.csv` after every step, $c$ and the chemical potential along the line
$y = 30$ to `OutputData/CahnHilliard_*_line_0011.csv` at the end, and the fields to
`OutputData/CahnHilliard` for ParaView. Every parameter at the top of the input can be changed from the
command line, for example a longer run on a coarser mesh:

```bash
<app> -i cahn_hilliard.i end_time=200 elements_per_side=75
```

`cahn_hilliard_phase_field.i` needs an application that includes the phase field module, such as
`phase_field-opt`, and `plot_cahn_hilliard.py` lists the runs it reads.

## Results id=results

### Spinodal decomposition and coarsening id=coarsening

Within the first few steps the noise is amplified at the wavelengths the instability favours and $c$
separates into a labyrinth of the two phases, close to $\pm 1$, divided by interfaces of the equilibrium
width. The labyrinth then coarsens: narrow fingers retract, small domains dissolve into their neighbours, and
the pattern grows steadily larger. The interfaces meet the walls at right angles, as the natural boundary
condition requires.

!media media/cahn_hilliard_snapshots.png
       id=snapshots-figure
       caption=The order parameter $c$ at $t = 20$, the end of the run the test makes, and at $t = 60$ and
               $t = 200$ from the same input run on with `end_time=200`.
       style=width:100%;margin-left:auto;margin-right:auto;

The mean of $c$ is $8.0330449 \times 10^{-4}$ at every step, to round-off, out to $t = 200$, and the total free
energy falls monotonically, from 893.0 after the first step to 644.2 at $t = 20$ and 335.4 at $t = 200$. Newton's
method converges quadratically from the previous step's solution, taking up to nine iterations in the first
steps, when the separation is fastest, and four or five once the pattern is coarsening.

### Comparison with the phase field module id=comparison

`cahn_hilliard_phase_field.i` solves the same problem with the phase field module, through the same
`Modules/PhaseField/Conserved` action and `REVERSE_SPLIT` formulation as `Parsed_SplitCH.i`, on the same mesh,
from the same nodal initial values and with the same boundary conditions, time integrator and step. The two
differ only in how the chemical potential is split and discretised.

!media media/cahn_hilliard_comparison.png
       id=comparison-figure
       caption=Left: the total free energy from the two solvers. Right: $c$ along the line $y = 30$ at
               $t = 20$.
       style=width:100%;margin-left:auto;margin-right:auto;

The free energies agree to 0.1% throughout and the profiles of $c$ along the centre line lie on top of each
other. The differences that remain are the discretisation error of the two formulations, which the table below
shows by repeating both at half and twice the resolution. The noise is sampled at the nodes, so each mesh
starts from different noise and finds a different pattern; what is compared is the two solvers on each mesh.

| Elements per side | RMS $\Delta c$ | max $\lvert \Delta c \rvert$ | $F_\text{MFEM}$ | $F_\text{PF}$ | $\Delta F / F_\text{PF}$ |
| -: | -: | -: | -: | -: | -: |
| 75 | $1.13 \times 10^{-1}$ | 1.27 | 635.16 | 638.37 | $-5.0 \times 10^{-3}$ |
| 150 | $1.65 \times 10^{-2}$ | 0.258 | 644.20 | 644.80 | $-9.3 \times 10^{-4}$ |
| 300 | $3.82 \times 10^{-3}$ | 0.044 | 667.03 | 667.18 | $-2.1 \times 10^{-4}$ |

Here $\Delta c = c_\text{MFEM} - c_\text{PF}$ and $\Delta F = F_\text{MFEM} - F_\text{PF}$, at $t = 20$,
over every node. Both differences fall roughly fourfold with each halving of the element size, the second
order convergence expected of bilinear elements, so the two solvers converge to the same solution. The
largest nodal differences are at interfaces, where a small shift in position is a large change in $c$: on
the coarsest mesh, where an interface spans only a few elements, the two solutions place some interfaces
differently and the largest difference is most of the jump across one. On the mesh used here it is 0.26,
and the RMS difference is under 2%.

The phase field module took 26 seconds and MFEM 34 seconds for the 150 element run on one process, both
factorising their Jacobians with MUMPS.

## References id=references

- J. W. Cahn and J. E. Hilliard, *Free energy of a nonuniform system. I. Interfacial free energy*,
  J. Chem. Phys., 28, 258 (1958).
- [Phase field module: solving the Cahn-Hilliard equation](https://mooseframework.inl.gov/modules/phase_field/Solving.html)
  - the split formulation and the example this one follows.
