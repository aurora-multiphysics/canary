# HTS Tape AC Loss: The T-A Formulation

This is the first benchmark of the [HTS Modelling Workgroup](https://htsmodelling.com/benchmark-1), a single
high temperature superconducting tape carrying a sinusoidal transport current, solved with the T-A formulation
(Zhang et al.; Liang et al.). The benchmark, the power law the superconductor follows and the Norris formula
it is checked against are set out in [HTS Tape AC Loss](hts_tape_ac_loss.md), which solves it with the
$\vec{H}$ formulation on the same mesh. This page covers what the T-A formulation does differently, and
compares the two.

The T-A formulation is the usual one for coated conductors, whose superconducting layer is a hundred to ten
thousand times thinner than it is wide. It treats that layer as a thin strip, with the current density
uniform through its thickness, and splits the problem in two: a current vector potential $T$ in the strip,
from which the current density follows, and a magnetic vector potential $A$ everywhere, from which the
field follows. Free space then only has to satisfy magnetostatics. It needs no artificial resistivity,
carries no current, and has no time derivative, which is what makes the formulation cheap for the large
coils and stacks of tape it was devised for.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Formulation id=formulation

As in the $\vec{H}$ formulation, only the cross-section of the tape is modelled, with the current along $z$
and the field in the $xy$ plane. The current density in the tape and the flux density everywhere are

!equation
J_z = \nabla \times \vec{T}, \qquad \vec{B} = \nabla \times (A \hat{z}),

with $\vec{T}$ in $H(\text{curl})$ on the tape and $A$ in $H^1$ on the whole cross-section. $A$ satisfies
magnetostatics with the current in the tape as its source, and $\vec{T}$ satisfies Faraday's law in the tape
with $\vec{E} = \rho(\lvert \vec{J} \rvert) \vec{J}$ and the power law resistivity $\rho = (E_c/J_c)
(\lvert \vec{J} \rvert / J_c)^{n-1}$:

!equation
\left( \nu \nabla A, \nabla A' \right) - \left( \nabla \times \vec{T}, A' \right)_\text{tape} = 0,
\qquad
\left( \rho \, \nabla \times \vec{T}, \nabla \times \vec{T}' \right)_\text{tape}
+ \left( \frac{\partial A}{\partial t}, \nabla \times \vec{T}' \right)_\text{tape} = 0,

with $\nu = 1/\mu_0$. The second is $E_z + \partial A / \partial t = -\partial V / \partial z$ tested against
$\nabla \times \vec{T}'$: the gradient of the electric scalar potential along the tape is uniform across it,
and drops out of the weak form because no $\vec{T}'$ has a tangential component on the tape's boundary.
Neither the voltage nor a constraint on the current appears; the current is fixed through $\vec{T}$ instead
([see below](#bcs)).

The equations are assembled as one block system, with the power law term by
[MFEMNLCurlCurlKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMNLCurlCurlKernel.html), the
coupling into the $A$ equation by
[MFEMMixedScalarCurlKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMMixedScalarCurlKernel.html)
and the coupling into the $\vec{T}$ equation by
[MFEMMixedScalarWeakCurlKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMMixedScalarWeakCurlKernel.html),
whose trial variable is the time derivative of $A$.

### The Thin Strip id=thin-strip

The thin strip approximation takes the current density as uniform through the thickness of the tape, so
that $\vec{T} = T_y(x) \hat{y}$, normal to the tape, and $J_z = \partial T_y / \partial x$. It is usually made
by collapsing the tape onto a line, and solving for $T$ on that line coupled to $A$ on the surrounding
surface. MOOSE's MFEM interface has no coupling between a variable on a line and one on the surface it lies
in, so the approximation is made instead on the tape's own 1 µm cross-section, meshed one element thick.

On a lowest order $H(\text{curl})$ element, $\vec{T}$ has one degree of freedom on each edge, its tangential
component there. Holding the tangential $\vec{T}$ to zero on the faces of the tape leaves only the edges
that cross it free, which are the values of $T_y$ across the width, and each element of the tape then
carries the constant current density $(T_{y,\text{right}} - T_{y,\text{left}})/h$. That is exactly the thin
strip, with the tape's true thickness kept for the field around it. The input relies on the tape being one
element thick, which is the default of `hts_tape_mesh.py`: with more, $\vec{T}$ also has components along
the tape on the edges inside it, and a null space of gradients there that nothing fixes.

### Free Space id=free-space

$\vec{T}$ has no meaning outside the tape, but MFEM's finite element spaces span the whole mesh, so its edges
in free space need an equation of their own. They are given $(\vec{T}, \vec{T}')_\text{air} = 0$, and since
the tangential $\vec{T}$ is zero on every boundary that free space shares with the tape, $\vec{T}$ and the
current density are exactly zero there. The edges that carry $\vec{T}$ in the tape belong to no element of
free space, so this has no effect on the tape whatever it is.

## Symmetry and Boundary Conditions id=bcs

Only the quarter of the cross-section with $x \ge 0$ and $y \ge 0$ is modelled, on the mesh of the
$\vec{H}$ formulation ([see there](hts_tape_ac_loss.md#mesh)). $A$ is even in both $x$ and $y$, so its normal
derivative vanishes on both axes, which is the natural boundary condition. On the arc of radius $R = 60$ mm
that bounds free space, $A = 0$, which makes the field there tangential to the arc, as that of a line
current on the axis would be. This is the same approximation as the $\vec{H}$ formulation makes on the arc.

The current through the quarter of the tape is the circulation of $\vec{T}$ round its boundary, $I(t)/4$.
The tangential $\vec{T}$ is zero on the face and the end of the tape and on its centre plane $y = 0$, through
none of which current flows, and on its centre line $x = 0$ it is

!equation
T_y = -\frac{I(t)}{4 b},

with $b$ the half-thickness. The current is therefore imposed exactly, with none of it carried by free space,
and the `TapeCurrent` postprocessor reports $I(t)$ to rounding. `hts_tape_mesh.py` names the boundaries this
needs, `tape_boundary` for the face and the end of the tape and `tape_centre_line` for the edge where it meets
the $y$ axis.

Only the difference between the values of $T_y$ at the two ends of the tape matters, and the time-dependent
value is put at the centre rather than at the end for the sake of Newton's method. Each time step starts from
the previous solution with the new boundary value, so the whole of the step's increment in $T$ first appears
as current in the single element next to that boundary. At the end of the tape, that element is 9.4 µm wide
and already at $J_c$: the increment adds $5 J_c$ to it, and Newton then spends twenty or so iterations working
its way back down the power law, the residual falling by $(1 - 1/n)^n \approx 0.36$ in each. At the centre the
element is 187 µm wide and below $J_c$, and the increment adds a quarter of $J_c$ to it.

!listing canary/test/tests/superconductivity/hts_tape_ac_loss/hts_tape_ta.i

## Solving id=solver

### Second Order $A$ id=second-order

$A$ is second order. The $\vec{T}$ equation sees $A$ only through its average over each element of the
tape, and at first order that average barely registers a current density alternating from one element to
the next. In the middle of the tape, where the current has not penetrated, the resistivity is negligible
and nothing else holds the current density in place, and at first order it comes out ragged by 5% of $J_c$
from one element to the next. At second order its roughness, the sum of the second differences of the
element values, is a hundred and eighty times smaller and no larger than that of the $\vec{H}$ formulation,
and the profile agrees with the $\vec{H}$ formulation's to 0.1% of $J_c$ away from the current front. The
loss is little affected: 4.1624e-4 J/m at first order and 4.1645e-4 J/m at second, against 4.1663e-4 J/m for
the $\vec{H}$ formulation. Second order doubles the cost of each Newton iteration.

### Newton's Method and Its Residual id=newton

The linear systems are solved directly with MUMPS: they couple a magnetostatic problem to one whose
resistivity spans tens of orders of magnitude, and are small, about 10000 unknowns.

The two equations are in different units, and left as they are, a relative tolerance on Newton's method
sees only one of them. A row of the $\vec{T}$ equation is the difference in $E + \partial A/\partial t$
between neighbouring elements of the tape, times the half-thickness and the time step, so that a mismatch
of $E_c$ is $2.5 \times 10^{-15}$; a row of the $A$ equation is a current imbalance at a node times the time
step, and at the start of a step these are around $10^{-5}$. A relative tolerance of $10^{-10}$ on the
whole then leaves the $\vec{T}$ equation unconverged by more than half of $E_c$. Each equation is therefore
divided through by the tolerance it is to be solved to, $10^{-3} E_c$ for the $\vec{T}$ equation and
$10^{-6} I_c$ for the $A$ equation, and Newton stops when the residual falls below one.

The current tolerance is far below the current that any step adds, so that a step cannot pass before
Newton has iterated, which with a looser one some steps of the ramp did. The field tolerance is well clear
of the rounding floor of the $\vec{T}$ equation, near $10^{-5} E_c$, at which the line search starts to fail,
and the loss is the same to eleven figures for field tolerances of $10^{-2}$, $10^{-3}$ and $10^{-4} E_c$. At
$n = 25$ Newton then takes one to eight iterations a step, three or four in most, and the backtracking line
search is used in a handful of steps.

### The Ramp id=ramp

The first fortieth of a cycle is taken in steps ten times shorter. From the virgin state the current grows by
half or more in each of the first few steps of 400 to the cycle, the front up to which it has penetrated the
tape moves further in a step than the elements at the edge are wide, and Newton stalls in several of them.
In steps of 4000 to the cycle it does not. The state at the first peak, and so the loss, which is taken from
there on, does not depend on how the ramp is taken.

## Files id=files

| File | Role |
| :- | :- |
| `hts_tape_ta.i` | The solve. This is the input to run. |
| `hts_tape.i` | The same problem with the $\vec{H}$ formulation, documented [separately](hts_tape_ac_loss.md). |
| `hts_tape_mesh.py` | Builds `hts_tape.e`, shared by both, with gmsh. Run with `--help` for the meshing controls. |
| `run_ac_loss_sweep.py` | Runs the whole benchmark with either formulation, extrapolates in time, and writes the tables and figures. |

## Running id=running

```bash
cd test/tests/superconductivity/hts_tape_ac_loss
moose_test-opt -i hts_tape_ta.i                                    # I_m/I_c = 0.5, n = 25
moose_test-opt -i hts_tape_ta.i current_ratio=0.9 n_value=101
./run_ac_loss_sweep.py --formulation ta --exe moose_test-opt --plot  # the whole benchmark
```

Each run takes about seventy seconds on one core. `OutputData/HTSTapeTACSV.csv` holds the current in the
tape, the power dissipated in it and the loss per cycle at every step, the last only complete at the end of
the run, and `OutputData/HTSTapeTAPeakProfile_current_density_profile_0190.csv` the current density along
the middle of the tape at the first peak of the current, which is the 190th step.

## Results id=results

The loss per cycle, extrapolated in time from 400 and 800 steps per cycle as for the $\vec{H}$ formulation,
against the Norris formula and against the $\vec{H}$ formulation's extrapolated loss on the same mesh:

| $I_m/I_c$ | Norris (J/m) | $n = 25$ (J/m) | from Norris | from $\vec{H}$ | $n = 101$ (J/m) | from Norris | from $\vec{H}$ |
| :- | :- | :- | :- | :- | :- | :- | :- |
| 0.1 | 6.024e-7 | 7.933e-7 | +31.7% | +1.29% | 6.736e-7 | +11.8% | +0.76% |
| 0.2 | 9.757e-6 | 1.147e-5 | +17.6% | +0.50% | 1.042e-5 | +6.8% | +0.40% |
| 0.3 | 5.044e-5 | 5.543e-5 | +9.9% | +0.24% | 5.254e-5 | +4.2% | -0.42% |
| 0.4 | 1.644e-4 | 1.717e-4 | +4.5% | +0.17% | 1.700e-4 | +3.5% | +0.12% |
| 0.5 | 4.185e-4 | 4.183e-4 | -0.0% | -0.04% | 4.280e-4 | +2.3% | -0.19% |
| 0.6 | 9.176e-4 | 8.787e-4 | -4.2% | -0.22% | 9.287e-4 | +1.2% | -0.33% |
| 0.7 | 1.832e-3 | 1.674e-3 | -8.6% | -0.31% | 1.832e-3 | +0.0% | -0.43% |
| 0.8 | 3.461e-3 | 2.983e-3 | -13.8% | -0.39% | 3.409e-3 | -1.5% | -0.53% |
| 0.9 | 6.454e-3 | 5.094e-3 | -21.1% | -0.49% | 6.196e-3 | -4.0% | -0.74% |

!media media/hts_tape_ta_ac_loss.png
       id=ta-ac-loss
       caption=The loss per cycle with the T-A formulation against the amplitude of the transport current, for
               both exponents, with the Norris formula for a thin strip in the critical state.

The two formulations agree to within half a percent at $n = 25$ from $I_m/I_c = 0.2$ up, and the differences
from Norris are those of a finite exponent [described for the $\vec{H}$ formulation](hts_tape_ac_loss.md#results).
The T-A loss is slightly the lower at large amplitudes, where the current density in the penetrated bands is
above $J_c$, and slightly the higher at small ones. At $I_m/I_c = 0.1$ it is 1.3% the higher, where the current
penetrates only three of the edge elements: that is also where the $\vec{H}$ formulation's loss moves most
with the mesh, its refined meshes coming out 0.9 to 1.9% above its shipped one, so that the T-A value sits
among them.

### $n = 101$ id=n101

At $n = 101$ the agreement is to 0.8%, but with a caveat. In every run, once in the cycle, Newton's method runs
out of iterations in one to four consecutive steps, and `MFEMTransient` carries on regardless. It happens
when the reversed current, which re-enters the tape from its edge after each peak, crosses into the next of
the narrow elements at the edge, at a moment that depends on the amplitude, 6.45 ms after the start at
$I_m/I_c = 0.9$ and 7.4 ms at 0.3, and not on the time step: halving the step does not remove it. The element
the front enters has to swing from about $+0.9 J_c$ to $-J_c$ in one step, through the part of the power law
where it is flat, and the only other stiffness on its current is its own inductance, which for an element
9.4 µm wide is small. At $n = 25$ the power law is soft enough for Newton to manage, and no step of the sweep
fails to converge.

The steps that fail are few, Newton recovers within a few steps after them, and the losses stay within 0.8%
of the $\vec{H}$ formulation's, but the one that most visibly suffers, $I_m/I_c = 0.3$, comes out with its
loss at 800 steps below its loss at 400, against the trend of every other case. The $n = 101$ results should
be read with that in mind. `run_ac_loss_sweep.py` counts the failed steps of each run. The $\vec{H}$
formulation, in which the time derivative of the field stiffens the same elements, converges in every step
at $n = 101$ by its own test, in up to 22 iterations at $I_m/I_c = 0.9$.

!media media/hts_tape_ta_current_density.png
       id=ta-current-density
       caption=The current density along the middle of the tape at the first peak of the current, with the
               T-A formulation, against the critical state profile. The steps are the elements, on which the
               current density is constant.

The current density profiles at the first peak, which comes before any of the failed steps, have the shape of
the critical state and agree with those of the $\vec{H}$ formulation: at $n = 25$ and $I_m/I_c = 0.5$ to within
0.2% of $J_c$ away from the current front, both peaking at $1.07 J_c$ at the edge and carrying $0.33 J_c$ in the
middle.

## Convergence id=convergence

### Time Step id=timestep-convergence

Backward Euler is first order, and as with the $\vec{H}$ formulation, each halving of the step halves the
change in the loss, at $n = 25$ and $I_m/I_c = 0.5$:

| Steps per cycle | Loss (J/m) | Change |
| :- | :- | :- |
| 200 | 4.1460e-4 | |
| 400 | 4.1645e-4 | 1.85e-6 |
| 800 | 4.1738e-4 | 9.3e-7 |
| 1600 | 4.1784e-4 | 4.6e-7 |

Richardson extrapolation of the last two gives 4.1831e-4 J/m, and of 400 and 800 steps 4.1830e-4, against
4.1849e-4 for the $\vec{H}$ formulation. The shipped 400 steps, which the test runs, are 0.44% short of the
limit. No step of these runs fails to converge.

### Order of $A$ id=order-convergence

At 400 steps per cycle, $n = 25$ and $I_m/I_c = 0.5$:

| $A$ | Loss (J/m) | Roughness of the profile at the peak | Run time |
| :- | :- | :- | :- |
| first order | 4.1624e-4 | 3.9 | 34 s |
| second order (shipped) | 4.1645e-4 | 0.022 | 70 s |
| $\vec{H}$ formulation | 4.1663e-4 | 0.028 | 8 s |

The roughness is the sum of the absolute second differences of the element values of $J/J_c$ along the
profile, up to 5 mm from the centre ([see above](#second-order)).

## Notes id=notes

The test suite runs the $n = 25$, $I_m/I_c = 0.5$ case at 400 steps per cycle, and compares the history of
the tape current, the power and the loss, and the current density profile at the first peak. The rest of the
benchmark is `run_ac_loss_sweep.py --formulation ta`, which takes about an hour; run the $\vec{H}$ formulation's
sweep first for the comparison column.

The T-A formulation is the slower of the two here, by about a factor of nine. Its unknowns are not fewer:
$\vec{T}$ spans the whole mesh though it is used only in the tape, and second order $A$ has about four times
as many unknowns as the mesh has nodes. Its cost advantage is in problems with many tapes, where the
$\vec{H}$ formulation has to mesh every one of them through its thickness and the T-A formulation only along
its width, and free space is where most of the unknowns are.

## References id=references

- [HTS Modelling Workgroup, Benchmark #1](https://htsmodelling.com/benchmark-1) - the specification.
- H. Zhang, M. Zhang and W. Yuan, *An efficient 3D finite element method model based on the T-A
  formulation for superconducting coated conductors*, Supercond. Sci. Technol., 30, 024005 (2017).
- F. Liang, S. Venuturumilli, H. Zhang, M. Zhang, J. Kvitkovic, S. Pamidi, Y. Wang and W. Yuan, *A finite
  element model for simulating second generation high temperature superconducting coils/stacks with large
  number of turns*, J. Appl. Phys., 122, 043903 (2017).
- W. T. Norris, *Calculation of hysteresis losses in hard superconductors carrying ac: isolated conductors and
  edges of thin sheets*, J. Phys. D: Appl. Phys., 3, 489 (1970).
