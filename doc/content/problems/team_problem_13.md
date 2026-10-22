# TEAM Problem 13: The Saturating Steel Channel

[TEAM](https://www.compumag.org/wp/team/) Problem 13 is a non-linear magnetostatic benchmark from the
COMPUMAG workshop series. A racetrack coil drives flux up a steel plate, across a 0.5 mm air gap and back
down two steel channels whose open sides face the plate. The excitation is large enough to saturate the
steel, so the reluctivity depends on the local flux density and the problem is genuinely non-linear.

What makes it a good test is that the two things which set the answer pull in opposite directions. The gap
and the fringing around it carry most of the reluctance of the circuit, and they are three-dimensional; the
steel carries almost none of it until it saturates, at which point its reluctivity rises by a factor of
fifteen over a few tenths of a tesla. Getting the flux right needs both the geometry and the B-H curve
handled properly.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

All dimensions are from Fig. 1 of the problem statement. The steel is 3.2 mm thick throughout.

| Quantity | Value |
| :- | :- |
| Centre plate | 3.2 mm ($x$) by 50 mm ($y$) by 126.4 mm ($z$), centred on the origin |
| Channels | C-sections, 3.2 mm thick, 50 mm wide in $y$, 126.4 mm tall |
| Channel flanges | reach in to $x = \pm 2.1$ mm, at $\lvert z \rvert = 60$ to $63.2$ mm |
| Channel webs | at $\lvert x \rvert = 122.1$ to $125.3$ mm, full height |
| Channel offset | one at $15 \le y \le 65$ mm, the other at $-65 \le y \le -15$ mm |
| Plate-to-channel gap | 0.5 mm nominal |
| Coil | rounded rectangle, 150 mm inner opening with R25 corners, 25 mm thick, 100 mm tall |
| Excitation | 1000 AT (the benchmark also specifies 3000 AT) |

The two channels are not aligned with each other: one is rotated 180 degrees about the $z$ axis relative to
the other, which is why they sit on opposite sides in $y$. The coil surrounds the centre plate and passes
inside both channels, so the flux it drives up the plate divides between the two channels and returns down
their webs.

### The B-H Curve id=bh-curve

The steel is specified by 27 tabulated $(B, H)$ points up to $B = 1.8$ T and, above that, by Eq. (1) of the
problem statement,

!equation
B = \mu_0 H + \left( a H^2 + b H + c \right), \qquad 1.8 \le B \le 2.22\ \text{T},

with $a = -2.381 \times 10^{-10}$, $b = 2.327 \times 10^{-5}$ and $c = 1.590$, becoming
$B = \mu_0 H + M_s$ once the steel is fully saturated.

Table 4 of the problem statement asks participants which interpolation they used between the tabulated
points and lists splines first, so `bh_curve.py` fits a monotone cubic (Fritsch-Carlson) through all 27
points and joins it at $B = 1.8$ T to Eq. (1) inverted for $H(B)$. It reproduces every tabulated point
exactly and $H(B)$ is monotone throughout, which is what makes the Jacobian below positive definite.

Two small liberties are taken with the specification, both reported by `./bh_curve.py --check`:

- Eq. (1) gives $H = 9424.35$ A/m at $B = 1.8$ T where the table says 9423, a difference of 0.014%. The
  last tabulated point is moved onto Eq. (1) so that the two branches join continuously.
- The parabola $aH^2 + bH + c$ peaks at 2.158556 T rather than at the quoted $M_s = 2.16$ T. Taking its own
  peak as the saturation magnetisation keeps $B(H)$ continuously differentiable where the branches meet;
  using 2.16 exactly would put a 1.4 mT step in the curve and make $H(B)$ non-monotone.

The resulting reluctivity $\nu = H/B$ runs from 6400 m/H as $B \rightarrow 0$, down to a minimum of 341 m/H
at $B \approx 1.05$ T, and back up past 19000 m/H at $B = 2.2$ T.

## Symmetry and Mesh id=mesh

The geometry and the coil current are unchanged by reflection in $z = 0$, so only the half of the model
above that plane is meshed, inside a 500 mm by 500 mm by 200 mm box of free space. The problem statement
also allows a quarter model, using the 180 degree rotational symmetry about the $z$ axis, but that needs a
periodic boundary condition; the half model needs only the reflection.

`team13_mesh.py` builds the mesh with [gmsh](https://gmsh.info). The coil is a filleted box with a smaller
filleted box cut out of it, which reproduces the rounded rectangle of the problem statement exactly: the
corner arcs of a box of half-width $a$ filleted at radius $r$ are centred at $(\pm(a-r), \pm(a-r))$. Element
size is set by a small number of axis-aligned box fields, finest in the two 0.5 mm gaps, next finest through
the 3.2 mm steel, and coarsening out to the truncation boundary.

The blocks are named for what they are used for. `Coil` and `Air` are free space; `Channel` is what is left
of the steel of the two channels. The rest of the steel is cut into 25 measurement blocks, one centred on
each of the sections where Table 1 of the problem statement reports a flux density: `PlateZ000` to
`PlateZ060` up the centre plate, `FlangeX002` to `FlangeX122` along the top flange of the $+x$ channel, and
`WebZ000` to `WebZ060` down its web. The other channel is the mirror of that one and carries no blocks of its
own. Two sidesets are named: `symmetry`, the $z = 0$ plane, and `far_field`, the outer faces of the free
space box.

Every one of those blocks is steel and has to be listed in the input's `steel_blocks`, since that is what the
reluctivity and the non-linear curl-curl kernel are applied to. A block left out of it gets no curl-curl term
at all and the field in it diverges.

MFEM reads named blocks and sidesets from Exodus files only - its gmsh reader gives bare integer attributes
- so the script runs the mesh back through MOOSE with `--mesh-only` to convert it. That step needs a MOOSE
executable; the rest needs only gmsh.

## Formulation id=formulation

The problem is solved for the magnetic vector potential $\vec{A}$, with $\vec{B} = \nabla \times \vec{A}$, in
$H(\text{curl})$:

!equation
\left( \nu\!\left(\left|\nabla \times \vec{A}\right|\right) \nabla \times \vec{A}, \nabla \times \vec{A}' \right)
= \left( \vec{J}_s, \vec{A}' \right)

with $\nu$ taken from the [B-H curve](team_problem_13.md#bh-curve) inside the steel and $1/\mu_0$ elsewhere.
`MFEMNLCurlCurlKernel` asks for $\nu$ and for $\left|\vec{B}\right| \partial \nu / \partial \left|\vec{B}\right|$,
which the input forms from $H(B)$ and $H'(B)$ as

!equation
\nu = \frac{H(B)}{B}, \qquad
B \frac{\partial \nu}{\partial B} = H'(B) - \frac{H(B)}{B}.

Both are finite at $B = 0$, where they tend to $H'(0)$ and to zero, but neither expression can be evaluated
there, so the input switches to those limits below $B = 10^{-6}$ T. This matters more than it looks: floor
the denominator instead and $\nu$ comes out as zero rather than 6400 m/H wherever the field vanishes, which
leaves the steel with no reluctivity at all on the first Newton iteration and the solve never recovers.

Because $H(B)$ is monotone, $\nu + B\,\partial\nu/\partial B = H'(B) > 0$, so the Jacobian stays positive
definite even where $\nu$ is falling, and the solve converges quadratically - ten Newton iterations from a
zero initial guess.

### The Coil Current id=coil-current

The benchmark's coil is a winding carrying a uniform ampere-turn density. With $d(x,y)$ the distance from
the square core the corner arcs are centred on, that current is

!equation
\vec{J}_\text{applied} = \frac{NI}{A}\, \hat{z} \times \nabla d,

of magnitude $NI/A$ everywhere in the winding, $A$ being the 25 mm by 100 mm cross-section. The half model
holds half the coil and its mirror image supplies the other half, so the current density is the same as in
the full problem.

Writing that down is not the end of it. The curl-curl operator has every gradient in its null space, so its
right hand side has to be orthogonal to every discrete gradient, and the projection of $\vec{J}_\text{applied}$
onto the mesh is not: the meshed coil boundary is faceted, so a little current leaks across it. The
sub-app `team13_coil_source.i` takes that leak back out. It solves, over the whole mesh,

!equation
\left( \nabla\phi, \nabla\psi \right) = \left( \vec J_\text{applied}, \nabla\psi \right)_\text{coil}
\qquad \text{for all } \psi,

so that the current the main app uses, $\vec{J}_\text{applied}$ over the coil minus $\nabla\phi$ over the
whole domain, satisfies $(\vec{J}, \nabla\psi) = 0$ for every $\psi$ by construction. This is a Helmholtz
projection - the smallest correction that makes the source discretely divergence free - so the current stays
as close to the specified uniform density as it can. Solving a conduction problem in the coil instead, the
way an open coil source does, would drive the current to the inside of the bends and stop representing a
winding.

$\phi$ is pinned on `far_field`, which is not only to make that solve non-singular: it matches the main
problem, whose tangential condition on the same boundary leaves exactly the gradients of potentials constant
there in its null space. The two pieces reach the main app separately, as `applied_current` and
`correction_potential`, because the first is integrated over the coil alone and the second over the whole
domain; a single combined field would have to carry the coil's sharp edge, which an $H(\text{curl})$ field
cannot.

The current does not depend on the solution, so the sub-app runs once, at `INITIAL`.

!listing canary/test/tests/electromagnetics/team/problem_13/team13_coil_source.i

### Boundary Conditions id=boundary-conditions

Reflecting in $z = 0$ leaves the coil current unchanged, so $B_z$ is even there and $B_x$ and $B_y$ are odd:
the field is normal to the symmetry plane. That is the natural condition of this weak form,
$\vec{n} \times \vec{H} = 0$, so the plane is left free and carries no constraint at all. The far-field faces
carry $\vec{n} \times \vec{A} = 0$, which is $\vec{B} \cdot \vec{n} = 0$, holding the flux inside the box.

### Removing the null space id=null-space

A compatible right hand side makes the problem solvable, but $\vec{A}$ is still only determined up to a
gradient, and the Newton iterates wander along that null space until $\nabla \times \vec{A}$ is lost to
cancellation. A small mass term $\varepsilon(\vec{A}, \vec{A}')$ removes it. At $\varepsilon = 10^2$ it is
around six orders of magnitude below the curl-curl term at the scale of this model and does not measurably
perturb $\vec{B}$: the reported flux densities are unchanged to four decimal places over
$1 \le \varepsilon \le 10^4$.

With the null space gone the system is no longer singular, and `MFEMHypreAMS` must be told so. Leaving
`singular = true` set costs AMS most of its effectiveness here - the inner solves stop converging, the Newton
steps become inaccurate, and the non-linear iteration stalls at a few percent residual instead of converging.

That the source really is divergence free shows up in $\vec{A}$. A mass term of this kind absorbs whatever
part of the right hand side lies in the null space, and turns it into a gradient of size proportional to
$1/\varepsilon$. With the corrected source there is nothing to absorb, and $\int|\vec{A}|^2$ comes out at
$2.9096 \times 10^{-9}$, $2.9095 \times 10^{-9}$ and $2.9083 \times 10^{-9}$ at $\varepsilon = 10^2$, $10^3$
and $10^4$ - the same number. With the analytic current projected straight onto the mesh, the same three runs
gave $3.2674 \times 10^{-1}$, $3.2674 \times 10^{-3}$ and $3.2677 \times 10^{-5}$, exactly the
$1/\varepsilon^2$ of a right hand side that is not compatible. Correcting the source moved the reported flux
densities by about 0.15%.

#### Other ways to remove the null space id=null-space-alternatives

Two others were tried, both on the shipped mesh and both before the source was corrected, so the numbers
below go with the analytic source.

A tree-cotree gauge, `MFEMTreeCotreeGaugeEssentialConstraint`, fixes the lowest order edge degrees of freedom
on a spanning tree of the mesh to zero, removing the null space outright rather than perturbing the operator.
It agreed with the mass term on every reported flux density to $3 \times 10^{-4}$ and on the field energy to
$1 \times 10^{-3}$, which was a useful check on both. It also brought $\int|\vec{A}|^2$ down from
$3.27 \times 10^{-1}$ to $1.45 \times 10^{-7}$, which was the reason to want it at the time - with an
incompatible source the mass term leaves a potential some 2000 times the physical field, so `a_field` is only
meaningful through $\vec{B} = \nabla \times \vec{A}$. Correcting the source removes that motivation, since
the mass term now leaves no gradient at all. The gauge costs about 1.8 times the CG work, AMS being built for
the ungauged operator.

`MFEMHypreAMS`'s `block` parameter, which names the subdomains where the mass term is zero and hands hypre an
interior-node marker for them, does not work here. It is a hint to the preconditioner, not a gauge: it tells
AMS where the operator is singular but does not remove the null space. That is enough in the eddy current
problem it is meant for, where a conducting region carries a $\sigma\,\partial\vec{A}/\partial t$ term that
removes the null space there and the source lives inside that conductor. Neither holds here - nothing is
conducting - and with the whole domain marked the interior-node correction annihilates the Newton update
itself, leaving the residual at exactly its starting value before it blows up. Confining the mass term to the
steel and marking only the free space fails the same way, and adding the marker on top of a working mass term
changes nothing measurable.

### Newton damping id=damping

The reluctivity rises fifteenfold between $B = 1.5$ T and $1.8$ T, and an undamped Newton update aimed at an
element sitting in that knee overshoots to the far side of it. The next update brings it back, and the
iteration settles into a limit cycle rather than converging. It shows up whenever the discretisation resolves
the knee more sharply, or the steel is driven further into it:

| Case | Undamped | With backtracking |
| :- | -: | -: |
| 1000 AT, first order | 10 iterations, $1.1 \times 10^{-8}$ | unchanged: the full step is always accepted |
| 3000 AT, first order | 50 iterations, stalls at $3 \times 10^{-3}$ | 11 iterations, $7.6 \times 10^{-9}$ |
| 1000 AT, second order | 50 iterations, stalls at $8.8 \times 10^{-3}$ | 10 iterations, $3.2 \times 10^{-10}$ |

`line_search = backtracking` on `MFEMNewtonNonlinearSolver` shortens the update until the residual norm meets
an Armijo sufficient-decrease condition. The shipped 1000 AT first order case never needs it - the full step
always passes, and the results are unchanged to the last digit - but the two harder cases go from not
converging at all to converging in about the same number of iterations as the easy one. The second order case
went from 50 iterations and 56 minutes without converging to 10 iterations and 14 minutes with.

Worth knowing if the same thing shows up elsewhere: the stalled solutions were not wrong, just unconverged by
the residual test. Every reported value in the 50-iteration second order run agreed with the converged one to
$3.7 \times 10^{-5}$, and the 3000 AT run agreed to $2 \times 10^{-6}$ between iterations 20 and 50. The
field settles long before the residual does, because what is chattering is a handful of elements whose
reluctivity flips either side of the knee without moving the flux.

!listing canary/test/tests/electromagnetics/team/problem_13/team13.i

## Files id=files

| File | Role |
| :- | :- |
| `team13.i` | The magnetostatic solve. This is the input to run. |
| `team13_coil_source.i` | Sub-app: the divergence-free coil current. Launched automatically. |
| `team13.e` | The mesh. |
| `team13_mesh.py` | Builds `team13.e` with gmsh. |
| `bh_curve.py` | Turns the benchmark B-H curve into the expressions in `team13.i`. |
| `plot_team13_flux.py` | Draws both comparison figures below from the postprocessor CSV. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

All of them live in `test/tests/electromagnetics/team/problem_13`.

## Running id=running

The main app launches the coil source sub-app itself, so only one input needs to be given:

```bash
<path to a MOOSE app built with MFEM support> -i team13.i
```

or, from the repository root, through the test harness:

```bash
./run_tests --re TEAM_Problem_13
```

Results are written to `OutputData/`. `TEAM13CSV.csv` holds the 25 section averages of Table 1,
`TEAM13CSV_air_line_0001.csv` and `TEAM13CSV_recommended_points_0001.csv` hold Tables 2 and 3, and `TEAM13/`
holds a ParaView collection containing `a_field`, `b_field` and `b_magnitude`.

To rebuild the mesh, which needs gmsh and a MOOSE executable:

```bash
cd test/tests/electromagnetics/team/problem_13
./team13_mesh.py --exe /path/to/app-opt        # or set MOOSE_DIR
```

`--help` lists the element size controls.

The test runs the 1000 AT case. The benchmark's other excitation is a command line argument. The figures show
1000 AT only, that being the case there are measurements for, but the 3000 AT run can be drawn alongside:

```bash
<app> -i team13.i ampere_turns=3000 \
      Outputs/ReportedPostprocessors/file_base=OutputData/TEAM13CSV_3000AT
./plot_team13_flux.py --csv-3000 OutputData/TEAM13CSV_3000AT.csv
```

It converges in 11 Newton iterations and takes about the same time as the 1000 AT case. That is down to the
line search; see [Newton damping](team_problem_13.md#damping) below.

## Results id=results

Table 1 of the problem statement reports the average flux density in the steel over 25 sections, and Tables 2
and 3 report the flux density at points in the air and at four points recommended for comparison between
methods. All three are reported by the input.

Two choices go into reporting Table 1, and both are worth stating.

The first is what to average. Table 1 is headed "average flux density $\left|\vec{B}\right|$", but the
measurements behind it were taken with search coils wound round the steel, and a search coil reads the flux
through the section, which is the component of $\vec{B}$ along it and not the magnitude. What is reported
here is the component, $\vec{B} \cdot \hat{n}$ with $\hat{n}$ the direction the flux runs in that part of
the circuit: down the plate, along the flange towards the plate, and up the web. That makes every value
positive and directly comparable with the table. The distinction only matters where the flux is turning, but
there it matters a great deal - see [below](team_problem_13.md#turning).

The second is that a section is a plane, and an average over one is not something the postprocessors can take
on an interior face, so each is reported as the average over a block of steel centred on it. The flange and
web blocks are 5 mm bands. The plate blocks tile the plate in 10 mm slabs, the stations there being 10 mm
apart. Four bands sit at the end of their part and are clipped, so they are one-sided rather than centred:
positions 8 and 18 at the two ends of the flange, position 19 at the top of the web, and position 7 at the top
of the plate. Position 25, at the foot of the web, is clipped at the symmetry plane but is still centred, the
rest of its band lying below it.

The benchmark asks for both of its excitations, and both are tabulated below, but it only published
measurements at 1000 AT, so that is the only column the computed values can be checked against. The 3000 AT
column is what the same model gives with the steel driven much further into saturation. The 1000 AT and
3000 AT columns of Tables 1, 2 and 3 are blank in the problem statement - they are the tables participants
were asked to fill in - and the one set of readings it does print is headed "1000 AT, measured".

Of the 25 sections, Fig. 5 of the problem statement singles out three as the ones to compare against
experiment rather than against other methods: positions 1, 12 and 25, at the foot of the plate, half way
along the flange and at the foot of the web. It adds position 27 of Table 2, in the air, as a point where the
flux density is high and does not change sharply. Those four are the benchmark's own choice of where a
calculation should be judged, and the model is 8.4%, 9.3% and 8.0% low at the first three.

The measurements come in two sets, taken on specimens whose gaps came out at 0.52 mm and 0.47 mm against the
nominal 0.5 mm. The gap carries most of the reluctance of the circuit, so the spread between them is a fair
measure of how tightly the experiment pins the answer.

| No. | Section | 1000 AT (T) | 3000 AT (T) | $G = 0.52$ mm | $G = 0.47$ mm | Error at 1000 AT |
| :- | :- | -: | -: | -: | -: | -: |
| 1 | plate, $z = 0$ | 1.231 | 1.842 | 1.333 | 1.354 | -8.4% |
| 2 | plate, $z = 10$ mm | 1.220 | 1.834 | 1.329 | 1.339 | -8.6% |
| 3 | plate, $z = 20$ mm | 1.187 | 1.809 | 1.286 | 1.304 | -8.3% |
| 4 | plate, $z = 30$ mm | 1.130 | 1.759 | 1.225 | 1.245 | -8.5% |
| 5 | plate, $z = 40$ mm | 1.044 | 1.663 | 1.129 | 1.138 | -7.9% |
| 6 | plate, $z = 50$ mm | 0.911 | 1.466 | 0.985 | 0.982 | -7.4% |
| 7 | plate, $z = 60$ mm | 0.578 | 0.914 | 0.655 | 0.674 | -13.1% |
| 8 | flange, $x = 2.1$ mm | 0.329 | 0.526 | 0.259 | 0.263 | +25.9% |
| 9 | flange, $x = 10$ mm | 0.433 | 0.727 | 0.453 | 0.451 | -4.3% |
| 10 | flange, $x = 20$ mm | 0.515 | 0.901 | 0.554 | 0.563 | -7.8% |
| 11 | flange, $x = 30$ mm | 0.580 | 1.041 | 0.637 | 0.641 | -9.3% |
| 12 | flange, $x = 40$ mm | 0.637 | 1.167 | 0.698 | 0.706 | -9.3% |
| 13 | flange, $x = 50$ mm | 0.690 | 1.283 | 0.755 | 0.763 | -9.1% |
| 14 | flange, $x = 60$ mm | 0.740 | 1.391 | 0.809 | 0.819 | -9.1% |
| 15 | flange, $x = 80$ mm | 0.824 | 1.559 | 0.901 | 0.907 | -8.8% |
| 16 | flange, $x = 100$ mm | 0.869 | 1.620 | 0.945 | 0.958 | -8.7% |
| 17 | flange, $x = 110$ mm | 0.878 | 1.619 | 0.954 | 0.968 | -8.6% |
| 18 | flange, $x = 122.1$ mm | 0.883 | 1.604 | 0.956 | 0.968 | -8.2% |
| 19 | web, $z = 60$ mm | 0.886 | 1.600 | 0.960 | 0.971 | -8.3% |
| 20 | web, $z = 50$ mm | 0.890 | 1.610 | 0.965 | 0.973 | -8.2% |
| 21 | web, $z = 40$ mm | 0.896 | 1.619 | 0.970 | 0.982 | -8.2% |
| 22 | web, $z = 30$ mm | 0.902 | 1.629 | 0.974 | 0.985 | -7.9% |
| 23 | web, $z = 20$ mm | 0.906 | 1.636 | 0.981 | 0.991 | -8.1% |
| 24 | web, $z = 10$ mm | 0.909 | 1.640 | 0.984 | 0.995 | -8.1% |
| 25 | web, $z = 0$ | 0.911 | 1.641 | 0.985 | 0.995 | -8.0% |

Twenty-three of the twenty-five sit between 4% and 9% below the measurements, and flat across the whole
circuit: the plate, the flange and the web are all short by about the same amount, which is what a single
under-resolved reluctance looks like rather than an error local to one part. Refining the mesh moves the whole
set up together, as the [convergence study](team_problem_13.md#convergence) below shows.

!media media/team13_positions.png
       id=all-positions
       caption=All 25 sections of Table 1 at 1000 AT, against the position number the problem statement
               gives them. The three runs are the three parts of the magnetic circuit and are not
               joined to each other: position 7 is the top of the plate and position 8 the tip of the flange,
               with the 0.5 mm
               gap between them. The computed values track the measurements around the whole circuit and sit
               a little below them everywhere, which is the mesh resolution in the gap rather than anything
               local.
       style=width:100%;margin-left:auto;margin-right:auto;

### Where the flux turns id=turning

The two that stand out, positions 7 and 8, are the two sections where the flux is turning rather than running
along the part, and they are worth a note because they are also where the choice of reported quantity bites.

Averaging $\left|\vec{B}\right|$ rather than the component through the section makes no difference over most
of the circuit and a large one at those two:

| Section | Average component | Average $\lvert \vec{B} \rvert$ | Ratio |
| :- | -: | -: | -: |
| plate, $z = 0$ | 1.231 | 1.231 | 1.00 |
| plate, $z = 30$ mm | 1.130 | 1.131 | 1.00 |
| plate, $z = 50$ mm | 0.911 | 0.939 | 1.03 |
| plate, $z = 60$ mm | 0.578 | 0.747 | 1.29 |
| flange, $x = 2.1$ mm | 0.329 | 0.440 | 1.34 |
| flange, $x = 10$ mm | 0.433 | 0.484 | 1.12 |
| flange, $x = 60$ mm | 0.740 | 0.743 | 1.00 |
| web, $z = 0$ | 0.911 | 0.911 | 1.00 |

Position 8 is the extreme case, and the reason is the geometry. The plate is 50 mm wide, centred on $y = 0$,
and the flange is 50 mm wide starting at $y = 15$ mm, so the two face each other across the gap over only
$15 \le y \le 25$ mm - a fifth of the flange's width. All the flux crossing the gap at the tip goes through
that strip. Sampling across the width at $x = 2.2$ mm shows it plainly: 1.24 to 1.51 T over
$16 \le y \le 25$ mm, saturated, falling to 0.39 T by $y = 31$ mm and to between 0.03 and 0.10 T over the
outer half of the flange, which carries only transverse leakage. That leakage adds to an average of the
magnitude while carrying no flux through the section at all, which is why reporting the magnitude there
overstated the section by two thirds.

What is left at position 8 is the band average. Over the section at $x = 2.2$ mm the component comes to
0.279 T against a measured 0.261 T, so the field itself is right to about 7%, in line with the rest of the
circuit; but the section value climbs to 0.370 T by $x = 4.5$ mm, and averaging across the band takes the
reported figure to 0.329 T. Position 7 is the same story at the top of the plate, where the flux is on its way
out into the flange: its band runs from 55 mm to the top of the plate at 63.2 mm and averages over the whole
of that turn.

### Element order id=element-order

The offset above is spatial discretisation error, and raising the polynomial order removes most of it. Running
the same mesh with second order Nedelec elements for $\vec{A}$ - `fe_order=SECOND flux_order=FIRST`, which
also takes the H(div) and L2 spaces carrying the derived fields up an order - gives:

| | First order | Second order |
| :- | -: | -: |
| Mean error over the 23 well-posed sections | -8.2% | **-1.6%** |
| Range | -9.3% to -4.3% | -2.5% to +1.6% |
| Position 7, at the turn out of the plate | -13.1% | -3.8% |
| Degrees of freedom | 0.25M | 1.4M |
| Runtime | 1.3 min | 14 min |

That is a better answer than first order reached on any mesh in the
[study below](team_problem_13.md#convergence), including the 865k element one at -3.8%, and it comes from 210k
elements. It confirms what the mesh study implied: the flat few percent the first order results sit below the
measurements is discretisation error in the gap, not something wrong with the model.

Position 8 is the exception and moves the other way, from +25.9% to +31.7%. That is the expected sign. Its
error is not discretisation but the band average standing in for a plane section, and a more accurate solution
resolves the climb across that band better, so the average moves further from the value at the section itself.

Second order needs the line search to converge at all; see [Newton damping](team_problem_13.md#damping).

!media media/team13_flux.png
       id=flux-profile
       caption=Flux density up the centre plate at 1000 AT, positions 1 to 7 of Table 1. The computed profile
               has the shape of the measurement and sits a few percent below it on the shipped mesh; refining
               the mesh in the gap closes most of that. The last point, at z = 60 mm, averages over the turn
               out into the flange and falls further short.
       style=width:80%;margin-left:auto;margin-right:auto;

### Tables 2 and 3 id=point-values

The problem statement asks for two further sets, neither of which it gives measurements for - they are there
to be compared between methods. Both are point values rather than section averages, so unlike Table 1 the
magnitude is the right quantity for them. Both are sampled straight from
$\left|\nabla \times \vec{A}\right|$ rather than from the projected `b_magnitude`, so they carry no extra
projection error, and both are written to their own CSV files beside the postprocessor output.

Table 2, positions 26 to 36, is the flux density in the air under the flange, along the line $y = 20$ mm,
$z = 55$ mm. It falls from 34.4 mT at $x = 10$ mm to 1.6 mT at $x = 110$ mm: this is the leakage across the
channel window, strongest where the window is narrowest, next to the plate.

Table 3 is four points chosen because a method that is getting the field right should agree there.

| Point | Coordinates (mm) | Where | Computed (T) |
| :- | :- | :- | -: |
| 1 | (2.2, 15.1, 60.1) | in the flange, at the corner nearest the gap | 2.268 |
| 2 | (2.0, 14.9, 50.9) | in the gap, beside the plate | 0.030 |
| 3 | (1.5, 0, 55) | on the face of the plate, near the top | 0.323 |
| 4 | (1.5, 0, 25) | on the face of the plate, half way up | 1.151 |

Point 1 is deep into saturation, which is what the corner of the flange facing the gap should be. Point 4
sits between the plate averages either side of it, 1.187 T at $z = 20$ mm and 1.130 T at $z = 30$ mm, as it
should. Point 3 is much lower than the plate average at that height because it is on the face of the plate
rather than inside it, where the flux is already turning out towards the flange.

## Convergence id=convergence

Almost all of the discretisation error is in the 0.5 mm gap. The gap and the fringing around it carry most
of the reluctance of the circuit, so how finely they are resolved sets the flux; the steel carries little of
it and tolerates a coarse mesh. Element count on its own is a poor guide:

| Gap element size (mm) | Steel element size (mm) | Elements | Mean error, positions 1 to 6 |
| -: | -: | -: | -: |
| 1.2 | 4.0 | 144k | -12.2% |
| 0.9 | 3.0 | 320k | -8.9% |
| 0.7 | 2.3 | 679k | -6.9% |
| 0.25 | 4.0 | 204k | -7.4% |
| 0.25 | 3.0 | 383k | -5.6% |
| 0.15 | 2.5 | 865k | -3.8% |

The error is measured against the mean of the two specimens. Reading down the gap column rather than the
element column: 204k elements with a 0.25 mm gap mesh beat 679k elements with a 0.7 mm one, at less than a
third of the size. The shipped mesh is the fourth row, which runs in about a minute; halving the gap size
again roughly halves the remaining error for four times the elements.

What is left at the finest mesh tried is a few percent. Some of it is still discretisation, some is the free
space box being truncated 250 mm out, and some is the comparison itself: the reported values are averages
over a block of finite height rather than over the plane section a search coil encloses.

The number of Newton iterations is not uniform across those meshes. The shipped mesh converges in ten, and
the finest in a similar number, but the 679k mesh needed more than thirty. The solver has no line search, so
a mesh that resolves the corner of the flange end more sharply than it resolves the gap can put an iterate
deep into saturation and take a while to walk back out.

## Notes id=notes

The test suite runs the 1000 AT first order case. The second order case is a test too, marked heavy so that it
is only picked up by `./run_tests --heavy`, since it takes a quarter of an hour. The 3000 AT case is tabulated
above but is not a test and is left off the figures, there being no measurements to read it against.

## References id=references

<!-- textlint-disable terminology -->

- T. Nakata, N. Takahashi and K. Fujiwara, *Summary of results for benchmark problem 13 (3-D non-linear
  magnetostatic model)*, COMPEL, 11, 3, 345 (1992).
- [Problem statement](https://www.compumag.org/wp/wp-content/uploads/2018/06/problem13.pdf) - the geometry,
  the B-H curve, the measured flux densities and the reporting tables used above.

<!-- textlint-enable terminology -->
