# HTS Tape AC Loss: A Superconducting Tape Carrying an Alternating Current

This is the first benchmark of the [HTS Modelling Workgroup](https://htsmodelling.com/benchmark-1): a single
high temperature superconducting (HTS) tape carrying a sinusoidal transport current, for which the energy
dissipated per cycle, the AC loss, is to be computed. Coated conductor tapes of this kind are what the magnets
of compact fusion devices are wound from, and their AC loss is what sets the heat load those magnets put on
their cryogenics whenever their current changes.

A superconductor is not an ordinary conductor with a small resistance. Below its critical current density
$J_c$ it carries current with practically no electric field at all, and above it the field rises extremely
steeply. A transport current therefore does not spread across the tape as it would in copper: it enters at
the edges at $J_c$ and only fills the tape as it approaches the critical current $I_c$. The flux that has
entered the edges is pushed back out and reversed every half cycle, and that is where the loss comes from.
In the limit of an infinitely steep transition, the critical state, the loss of a thin strip has a closed
form due to Norris, which makes this benchmark a verification problem as well as a comparison between codes.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

The same benchmark is solved with the T-A formulation, which treats the tape as a thin strip and free space
as magnetostatics, in [HTS Tape AC Loss: The T-A Formulation](hts_tape_ac_loss_ta.md), on the same mesh.

## Problem Specification id=specification

| Quantity | Value |
| :- | :- |
| Tape cross-section | 12 mm wide, 1 µm thick, in free space |
| Critical current | $I_c = 300$ A, so $J_c = 2.5 \times 10^{10}$ A/m$^2$ |
| Electric field criterion | $E_c = 10^{-4}$ V/m |
| Power law exponent | $n = 25$; also run at $n = 101$ |
| Transport current | $I(t) = I_m \sin(2 \pi f t)$, $f = 50$ Hz |
| Amplitudes | $I_m / I_c = 0.1, 0.2, \ldots, 0.9$ |
| Reported | loss per cycle per unit length; current density along the middle of the tape at peak current |

The benchmark page states the geometry, $I_c$, the frequency and the amplitudes. The exponent and $E_c$ are
those of the published reference results by Grilli, which are given for $n = 25$ and
for $n = 101$, and $E_c = 10^{-4}$ V/m is the usual criterion.

### The Superconductor id=power-law

The tape follows the power law

!equation
\vec{E} = E_c \left( \frac{\lvert \vec{J} \rvert}{J_c} \right)^{n} \frac{\vec{J}}{\lvert \vec{J} \rvert},
\qquad
\rho(\lvert \vec{J} \rvert) = \frac{E_c}{J_c} \left( \frac{\lvert \vec{J} \rvert}{J_c} \right)^{n-1},

a resistivity that is negligible well below $J_c$ and rises as the 24th power of the current density above
it. As $n \rightarrow \infty$ it becomes the critical state: $\lvert \vec{J} \rvert = J_c$ wherever there is
an electric field, and $\lvert \vec{J} \rvert < J_c$ only where there is none.

### Analytic Reference id=norris

For a strip in the critical state, thin compared to its width, Norris found the loss per cycle per unit
length

!equation
Q = \frac{\mu_0 I_c^2}{\pi} \left[ (1 - F) \ln(1 - F) + (1 + F) \ln(1 + F) - F^2 \right],
\qquad F = \frac{I_m}{I_c}.

At the peak of the current, the current density across a strip of half-width $a$ is (Norris; Brandt and
Indenbom)

!equation
\frac{J(x)}{J_c} = \begin{cases}
  \dfrac{2}{\pi} \arctan \sqrt{\dfrac{a^2 - b^2}{b^2 - x^2}}, & \lvert x \rvert < b, \\[2ex]
  1, & b < \lvert x \rvert < a,
\end{cases}
\qquad b = a \sqrt{1 - F^2},

saturated at $J_c$ in the bands the current has penetrated and below it in the middle, which the current
has not yet reached.

Neither is the answer for a finite $n$. With the power law, the current density in the penetrated bands is
$J_c (E/E_c)^{1/n}$: less than $J_c$ where the electric field is weaker than $E_c$ and more than $J_c$ where
it is stronger. At small amplitudes the field is weak, the current spreads further in than the critical state
would put it, and the loss is higher than Norris's; at large amplitudes the reverse. Both differences close as
$n$ grows, which is the check made below.

## Formulation id=formulation

The tape is long, so only its cross-section in the $xy$ plane is modelled, with the current along $z$ and the
magnetic field in the plane. The problem is solved for the magnetic field $\vec{H}$ in $H(\text{curl})$, the
usual choice for superconductors (Brambilla et al.), with $J_z = \nabla \times \vec{H}$. Faraday's law
with $\vec{E} = \rho \vec{J}$ gives the weak form

!equation
\left( \mu_0 \frac{\partial \vec{H}}{\partial t}, \vec{H}' \right)
+ \left( \rho(\lvert \nabla \times \vec{H} \rvert) \nabla \times \vec{H}, \nabla \times \vec{H}' \right) = 0.

The curl-curl term is non-linear in the tape, where it is assembled by
[MFEMNLCurlCurlKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMNLCurlCurlKernel.html) from
$\rho$ and from $\lvert \vec{J} \rvert \, \partial \rho / \partial \lvert \vec{J} \rvert = (n-1) \rho$, which
is what the Newton Jacobian needs. The time derivative makes the operator definite everywhere, so there is no
null space to deal with.

### Free Space id=free-space

Free space has no resistivity, but the $\vec{H}$ formulation needs one there: it is the curl-curl term that
makes $\vec{H}$ curl free outside the tape. The resistivity only has to be large enough that free space carries
no current that matters, and making it much larger than that does harm.

The curl-curl term over an element of size $h$ is $\Delta t \, \rho / (\mu_0 h^2)$ times its mass term. Next
to the tip of the tape, $h$ is the tape's 0.5 µm half-thickness, and the commonly used 1 Ω m makes that ratio
$10^{16}$. The tape's own equations, which are almost pure mass term at that scale, then sit below the
rounding error of free space's, and Newton stops with the tape unconverged without noticing. The symptom is a
loss that moves by a percent or more with the Newton tolerance, or with any refinement of the mesh near the tip
([see below](#free-space-convergence)).

The resistivity of free space is therefore graded with distance $d$ from the tape,

!equation
\rho_\text{air} = \min \left( 1\ \Omega\,\text{m}, \ \max \left( 10^{-9}\ \Omega\,\text{m}, \ k d^2 \right) \right),
\qquad k = 10^4\ \Omega/\text{m},

which keeps the ratio at $\Delta t \, k d^2 / (\mu_0 h^2)$, around $10^3$ to $10^4$, wherever the elements
grow in proportion to their distance from the tape, and caps it far away. The current that free space carries
goes as the integral of $1/\rho_\text{air}$ over it, which stays small: where the transport current changes
fastest, free space carries about $3 \times 10^{-6}$ of its amplitude. With a uniform resistivity it would be
$2 \times 10^{-6}$ at 1 Ω m and $2 \times 10^{-4}$ at $10^{-2}$ Ω m, so a uniform value low enough to cure the
conditioning would carry too much current.

### Time Integration and the Loss id=loss

MFEMTransient integrates in time with backward Euler, 400 steps per cycle. The tape starts from the virgin
state, with no current, and is in the cyclic state from the first peak of the current on. The loss is taken
over the half cycle between the positive and negative peaks, $T/4 \le t \le 3T/4$, and doubled:

!equation
Q = 2 \int_{T/4}^{3T/4} \int_\text{tape} \vec{E} \cdot \vec{J} \, dA \, dt.

The power is integrated over the tape from an auxiliary variable holding $E_c J_c (\lvert J \rvert /
J_c)^{n+1}$, which is exact since the curl of the lowest order $H(\text{curl})$ space is constant on each
element, and over time as the sum of the power at the end of each step, which is what backward Euler
dissipates over it.

## Symmetry, Boundary Conditions and Mesh id=mesh

The tape and its current are symmetric about both axes, so only the quarter with $x \ge 0$ and $y \ge 0$ is
modelled. $J_z$ is even in $x$ and in $y$, so $H_x$ is odd in $y$ and $H_y$ odd in $x$, and the tangential field
vanishes on both axes: both carry $\vec{n} \times \vec{H} = 0$.

The free space is a quarter disc of radius $R = 60$ mm, and the transport current is imposed through the field
on its arc, which is taken as that of a line current on the axis,

!equation
\vec{H} = \frac{I(t)}{2 \pi r^2} (-y, x).

Its circulation round the boundary of the quarter is $I(t)/4$, so the tangential condition fixes the current
through the quarter whatever the field looks like inside it; only the shape of the field on the arc is
approximate, by the difference between a line current and a 12 mm strip seen from 60 mm. The lowest order
$H(\text{curl})$ space takes the tangential field at the midpoint of each boundary edge, which overestimates the
integral along an edge subtending an angle $\delta$ by a fraction $\delta^2/12$. The arc is cut into 1 mm
edges to put that at $2 \times 10^{-5}$, and the `TapeCurrent` postprocessor reports the current that actually
flows in the tape.

`hts_tape_mesh.py` builds the mesh with [gmsh](https://gmsh.info). The tape is a structured block of 100
quadrilaterals across its half-width and one through its half-thickness, graded so that the elements at its
edge, where the current enters and the field is singular, are 20 times narrower than those at its centre:
187 µm down to 9.4 µm. Free space is meshed with about 2700 triangles, which grow at a fixed rate with distance
from the tip of the tape, from its 0.5 µm half-thickness up to 4 mm. The blocks are `Tape` and `Air`, and the
boundaries `x_axis`, `y_axis` and `far_field`.

## Solver id=solver

The linear systems are solved directly with MUMPS. They are small, about 4400 unknowns, and a direct solve
keeps every Newton update exact. AMS preconditioned conjugate gradients converges on them too, now that free
space is graded, and gives the same loss to six figures, but takes around 300 iterations a solve and fourteen
times as long; with a uniform free space resistivity it stalls at a relative residual of $10^{-7}$, and the
inexact Newton updates that leaves converge only linearly.

Newton's method converges every step to a relative residual of $10^{-10}$, in two or three iterations. The
residual is dominated by free space, which is linear, so it falls by eight orders of magnitude well before the
tape has converged: at $10^{-8}$ the loss is still 0.2% from its converged value, and from $10^{-10}$ on it no
longer moves. The backtracking line search is never used at $n = 25$ and 400 steps per cycle, but at
$n = 101$, or with fewer steps, the full Newton step where the current front is moving can overshoot so far up
the power law that the residual overflows and the run fails. Backtracking brings it back.

!listing canary/test/tests/superconductivity/hts_tape_ac_loss/hts_tape.i

## Files id=files

| File | Role |
| :- | :- |
| `hts_tape.i` | The solve. This is the input to run. |
| `hts_tape_ta.i` | The same problem with the T-A formulation, documented [separately](hts_tape_ac_loss_ta.md). |
| `hts_tape_mesh.py` | Builds `hts_tape.e` with gmsh. Run with `--help` for the meshing controls. |
| `hts_tape.e` | The mesh, as built by `hts_tape_mesh.py` with its defaults. |
| `run_ac_loss_sweep.py` | Runs the whole benchmark, extrapolates in time, and writes the table and figures below. |

## Running id=running

```bash
cd test/tests/superconductivity/hts_tape_ac_loss
moose_test-opt -i hts_tape.i                         # I_m/I_c = 0.5, n = 25
moose_test-opt -i hts_tape.i current_ratio=0.9 n_value=101
./run_ac_loss_sweep.py --exe moose_test-opt --plot    # the whole benchmark
```

Each run takes about ten seconds on one core. `OutputData/HTSTapeCSV.csv` holds the current in the tape, the
power dissipated in it and the loss per cycle at every step, the last only complete at the end of the run, and
`OutputData/HTSTapePeakProfile_current_density_profile_0100.csv` the current density along the middle of the
tape at the first peak of the current.

To rebuild the mesh, which needs the gmsh Python module and a MOOSE executable to convert it to Exodus:

```bash
./hts_tape_mesh.py --exe /path/to/moose_test-opt
```

## Results id=results

The loss per cycle, extrapolated in time from 400 and 800 steps per cycle (see [below](#timestep-convergence)),
against the Norris formula:

| $I_m/I_c$ | Norris (J/m) | $n = 25$ (J/m) | difference | $n = 101$ (J/m) | difference |
| :- | :- | :- | :- | :- | :- |
| 0.1 | 6.024e-7 | 7.832e-7 | +30.0% | 6.685e-7 | +11.0% |
| 0.2 | 9.757e-6 | 1.142e-5 | +17.0% | 1.038e-5 | +6.4% |
| 0.3 | 5.044e-5 | 5.529e-5 | +9.6% | 5.276e-5 | +4.6% |
| 0.4 | 1.644e-4 | 1.714e-4 | +4.3% | 1.698e-4 | +3.3% |
| 0.5 | 4.185e-4 | 4.185e-4 | +0.0% | 4.288e-4 | +2.5% |
| 0.6 | 9.176e-4 | 8.807e-4 | -4.0% | 9.317e-4 | +1.5% |
| 0.7 | 1.832e-3 | 1.679e-3 | -8.3% | 1.840e-3 | +0.4% |
| 0.8 | 3.461e-3 | 2.995e-3 | -13.5% | 3.428e-3 | -1.0% |
| 0.9 | 6.454e-3 | 5.119e-3 | -20.7% | 6.243e-3 | -3.3% |

!media media/hts_tape_ac_loss.png
       id=ac-loss
       caption=The loss per cycle against the amplitude of the transport current, for both exponents, with the
               Norris formula for a thin strip in the critical state.

The differences from Norris are those [described above](#norris) of a finite exponent: higher than the critical
state at small amplitudes, where the electric field in the tape is weaker than $E_c$, lower at large ones,
where it is stronger, and crossing near $I_m/I_c = 0.5$. At $n = 101$ every one of them is smaller, by a factor
of three to six, as the critical state is approached. They are also what the reference results show. Those are
published only as figures, so the comparison below is with values read off them, to about 2% on the linear
plot at large amplitudes and to about 10% on the logarithmic one at small:

| $n$ | $I_m/I_c$ | reference, read from figure (J/m) | MOOSE (J/m) |
| :- | :- | :- | :- |
| 25 | 0.1 | 7.3e-7 | 7.83e-7 |
| 25 | 0.2 | 1.2e-5 | 1.14e-5 |
| 25 | 0.9 | 5.20e-3 | 5.12e-3 |
| 101 | 0.7 | 1.80e-3 | 1.84e-3 |
| 101 | 0.8 | 3.34e-3 | 3.43e-3 |
| 101 | 0.9 | 6.19e-3 | 6.24e-3 |

Below $I_m/I_c = 0.7$ the reference markers for the two exponents overlap each other and the Norris curve
too closely to be read separately.

!media media/hts_tape_current_density.png
       id=current-density
       caption=The current density along the middle of the tape at the first peak of the current, with the
               critical state profile. The steps are the elements, on which the current density is constant.

The current density profiles have the shape of the critical state, sharpening towards it as $n$ grows. In the
penetrated bands the current density sits above $J_c$, by 7% at $n = 25$ and $I_m/I_c = 0.5$ and by 13% at
0.9, and by 1% and 2% at $n = 101$, the field there being stronger than $E_c$ at the peak of the current. The
middle of the tape, which the current has not reached, carries $0.33 J_c$ at $I_m/I_c = 0.5$, the critical
state value of one third, and $0.66 J_c$ and $0.70 J_c$ at 0.9 for the two exponents against $0.71 J_c$. All of
these agree with the reference profiles.

## Convergence id=convergence

### Time Step id=timestep-convergence

Backward Euler is first order, and the loss converges as the step: each halving of the step halves the change,
at $n = 25$ and $I_m/I_c = 0.5$.

| Steps per cycle | Loss (J/m) | Change |
| :- | :- | :- |
| 100 | 4.1115e-4 | |
| 200 | 4.1480e-4 | 3.65e-6 |
| 400 | 4.1663e-4 | 1.83e-6 |
| 800 | 4.1756e-4 | 9.3e-7 |
| 1600 | 4.1803e-4 | 4.7e-7 |

Richardson extrapolation of the last two gives 4.1850e-4 J/m, and of 400 and 800 steps 4.1849e-4. The
shipped 400 steps, which the test runs, are 0.45% short of the limit; the results above are all extrapolated
from 400 and 800.

### Mesh id=mesh-convergence

At 400 steps per cycle and $n = 25$, varying the mesh controls of `hts_tape_mesh.py` one at a time from the
shipped mesh (100 elements across the half-width, one through the half-thickness, edge grading 20, free space
out to 60 mm):

| Mesh | $I_m/I_c = 0.1$ (J/m) | $I_m/I_c = 0.5$ (J/m) | $I_m/I_c = 0.9$ (J/m) |
| :- | :- | :- | :- |
| shipped | 7.808e-7 | 4.1663e-4 | 5.0928e-3 |
| 50 elements across | | 4.1688e-4 | |
| 200 elements across | 7.936e-7 | 4.1739e-4 | 5.1078e-3 |
| 400 elements across | 7.918e-7 | 4.1684e-4 | 5.1000e-3 |
| 2 elements through | 7.881e-7 | 4.1653e-4 | 5.0923e-3 |
| 4 elements through | | 4.1676e-4 | |
| edge grading 5 | | 4.1847e-4 | 5.1075e-3 |
| edge grading 50 | 7.938e-7 | 4.1750e-4 | |
| edge grading 200 | 7.958e-7 | | |
| free space to 120 mm | | 4.1654e-4 | 5.0923e-3 |

At half the critical current and above, the mesh moves the loss by no more than 0.4%. At $I_m/I_c = 0.1$ the
current penetrates only 30 µm from the edge, three of the shipped mesh's edge elements, and the shipped mesh
is 1.6% below the refined ones; that is small beside the 30% by which the loss there differs from Norris's. At
$n = 101$ the shipped mesh is within 0.7% of the refined ones at $I_m/I_c = 0.1$ and 0.2% at 0.5.

### Free Space Resistivity and Newton Tolerance id=free-space-convergence

The loss at $n = 25$ and $I_m/I_c = 0.5$, 400 steps per cycle, with the free space resistivity uniform and
graded, against the relative tolerance of Newton's method and the mesh at the tip of the tape:

| Free space | Newton tolerance | Mesh | Loss (J/m) | Steps not converged |
| :- | :- | :- | :- | :- |
| uniform 1 Ω m | $10^{-8}$ | shipped | 4.1736e-4 | 6 |
| uniform 1 Ω m | $10^{-8}$ | 400 elements across | 4.0857e-4 | 128 |
| uniform 1 Ω m | $10^{-8}$ | 2 elements through | 4.1292e-4 | 136 |
| uniform 0.1 Ω m | $10^{-6}$ | shipped | 4.1748e-4 | 0 |
| uniform 0.1 Ω m | $10^{-8}$ | shipped | 4.1713e-4 | 6 |
| uniform 0.1 Ω m | $10^{-10}$ | shipped | 4.1661e-4 | 280 |
| graded | $10^{-6}$ | shipped | 4.1749e-4 | 0 |
| graded | $10^{-8}$ | shipped | 4.1749e-4 | 0 |
| graded | $10^{-10}$ | shipped | 4.1663e-4 | 0 |
| graded | $10^{-12}$ | shipped | 4.1662e-4 | 6 |
| graded | $10^{-10}$ | 400 elements across | 4.1684e-4 | 0 |
| graded | $10^{-10}$ | 2 elements through | 4.1653e-4 | 0 |

With a uniform 1 Ω m, any refinement at the tip of the tape leaves Newton stalled at the rounding floor of free
space in a hundred or more steps, and the loss falls by 1 to 2%. With 0.1 Ω m the floor is lower, but the loss
only reaches its converged value at a tolerance that leaves Newton stalled on the floor in almost every step.
Graded, Newton converges everywhere, the loss is settled from $10^{-10}$ on, and it no longer depends on the mesh
at the tip. The steps that miss
$10^{-12}$ are those at the peaks of the current, where the step barely changes anything and the first
iteration already reaches the rounding floor.

## Notes id=notes

The test suite runs the $n = 25$, $I_m/I_c = 0.5$ case at 400 steps per cycle, and compares the history of
the tape current, the power and the loss, and the current density profile at the first peak. The rest of the
benchmark is `run_ac_loss_sweep.py`, which takes about ten minutes.

The benchmark also asks for the flux density along the middle of the tape. $\vec{H}$ is available for that, but
the field normal to a 1 µm tape is singular at its edges and is not reported here.

## References id=references

- [HTS Modelling Workgroup, Benchmark #1](https://htsmodelling.com/benchmark-1) - the specification.
- F. Grilli, [results for Benchmark #1](https://htsmodelling.com/hubfs/benchmarks/1/downloads/B1_results.pdf),
  Karlsruhe Institute of Technology - the reference losses and current density profiles for $n = 25$ and
  $n = 101$.
- W. T. Norris, *Calculation of hysteresis losses in hard superconductors carrying ac: isolated conductors and
  edges of thin sheets*, J. Phys. D: Appl. Phys., 3, 489 (1970).
- E. H. Brandt and M. Indenbom, *Type-II-superconductor strip with current in a perpendicular magnetic
  field*, Phys. Rev. B, 48, 12893 (1993).
- R. Brambilla, F. Grilli and L. Martini, *Development of an edge-element model for AC loss computation of
  high-temperature superconductors*, Supercond. Sci. Technol., 20, 16 (2007).
