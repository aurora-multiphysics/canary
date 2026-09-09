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
| Channel flanges | reach in to $x = \pm 2.1$ mm, at $\left|z\right| = 60$ to $63.2$ mm |
| Channel webs | at $\left|x\right| = 122.1$ to $125.3$ mm, full height |
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

The blocks are named for what they are used for. `Coil` and `Air` are free space; `Channel` is the steel of
the two channels. The rest of the steel is cut into measurement blocks, each centred on one of the sections
where the benchmark reports a flux density: `PlateZ000` to `PlateZ060` up the centre plate, `WebZ000` at the
foot of a web and `FlangeX060` halfway along a flange. Two sidesets are named: `symmetry`, the $z = 0$ plane,
and `far_field`, the outer faces of the free space box.

MFEM reads named blocks and sidesets from Exodus files only - its gmsh reader gives bare integer attributes
- so the script runs the mesh back through MOOSE with `--mesh-only` to convert it. That step needs a MOOSE
executable; the rest needs only gmsh.

## Formulation id=formulation

The problem is solved for the magnetic vector potential $\vec{A}$, with $\vec{B} = \nabla \times \vec{A}$, in
$H(\text{curl})$:

!equation
\left( \nu\!\left(\left|\nabla \times \vec{A}\right|\right) \nabla \times \vec{A}, \nabla \times \vec{A}' \right)
= \left( \vec{J}_s, \vec{A}' \right)

with $\nu$ taken from the [B-H curve](#bh-curve) inside the steel and $1/\mu_0$ elsewhere.
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
\left( \nabla\phi, \nabla\psi \right) = \left( \vec{J}_\text{applied}, \nabla\psi \right)_\text{coil}
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
gradient, and the Newton iterates wander along that null space until the curl is lost to cancellation. A
small mass term $\varepsilon(\vec{A}, \vec{A}')$ removes it. At $\varepsilon = 10^2$ it is around six orders
of magnitude below the curl-curl term at the scale of this model and does not measurably perturb $\vec{B}$:
the reported flux densities are unchanged to four decimal places over $1 \le \varepsilon \le 10^4$.

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
meaningful through its curl. Correcting the source removes that motivation, since the mass term now leaves no
gradient at all. The gauge costs about 1.8 times the CG work, AMS being built for the ungauged operator.

`MFEMHypreAMS`'s `block` parameter, which names the subdomains where the mass term is zero and hands hypre an
interior-node marker for them, does not work here. It is a hint to the preconditioner, not a gauge: it tells
AMS where the operator is singular but does not remove the null space. That is enough in the eddy current
problem it is meant for, where a conducting region carries a $\sigma\,\partial\vec{A}/\partial t$ term that
removes the null space there and the source lives inside that conductor. Neither holds here - nothing is
conducting - and with the whole domain marked the interior-node correction annihilates the Newton update
itself, leaving the residual at exactly its starting value before it blows up. Confining the mass term to the
steel and marking only the free space fails the same way, and adding the marker on top of a working mass term
changes nothing measurable.

!listing canary/test/tests/electromagnetics/team/problem_13/team13.i

## Files id=files

| File | Role |
| :- | :- |
| `team13.i` | The magnetostatic solve. This is the input to run. |
| `team13_coil_source.i` | Sub-app: the divergence-free coil current. Launched automatically. |
| `team13.e` | The mesh. |
| `team13_mesh.py` | Builds `team13.e` with gmsh. |
| `bh_curve.py` | Turns the benchmark B-H curve into the expressions in `team13.i`. |
| `plot_team13_flux.py` | Draws the comparison below from the postprocessor CSV. |
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

Results are written to `OutputData/`. `TEAM13CSV.csv` holds the reported flux densities and `TEAM13/` holds a
ParaView collection containing `a_field`, `b_field` and `b_magnitude`.

To rebuild the mesh, which needs gmsh and a MOOSE executable:

```bash
cd test/tests/electromagnetics/team/problem_13
./team13_mesh.py --exe /path/to/app-opt        # or set MOOSE_DIR
```

`--help` lists the element size controls. To change the excitation to the benchmark's other case, add
`ampere_turns=3000` on the command line.

## Results id=results

The benchmark reports the average flux density in the steel over a set of sections. Each section has a block
of its own in the mesh, centred on it, and the reported values are the average of $\left|\vec{B}\right|$ over
that block. The block for $z = 0$ is the top half of its band, the rest of it lying below the symmetry plane,
and the block for $z = 60$ mm runs from 55 mm to the top of the plate rather than being centred.

The measurements come in two sets, taken on specimens whose gaps came out at 0.52 mm and 0.47 mm against the
nominal 0.5 mm. The gap carries most of the reluctance of the circuit, so the spread between them is a fair
measure of how tightly the experiment pins the answer.

| Position | Section | Computed (T) | Measured, $G = 0.52$ mm | Measured, $G = 0.47$ mm |
| :- | :- | -: | -: | -: |
| 1 | plate, $z = 0$ | 1.236 | 1.333 | 1.354 |
| 2 | plate, $z = 10$ mm | 1.225 | 1.329 | 1.339 |
| 3 | plate, $z = 20$ mm | 1.192 | 1.286 | 1.304 |
| 4 | plate, $z = 30$ mm | 1.136 | 1.225 | 1.245 |
| 5 | plate, $z = 40$ mm | 1.053 | 1.129 | 1.138 |
| 6 | plate, $z = 50$ mm | 0.941 | 0.985 | 0.982 |
| 7 | plate, $z = 60$ mm | 0.741 | 0.655 | 0.674 |
| 14 | flange, $x = 60$ mm | 0.765 | 0.809 | 0.819 |
| 25 | web, $z = 0$ | 0.929 | 0.985 | 0.995 |

!media media/team13_flux.png
       id=flux-profile
       caption=Flux density up the centre plate at 1000 AT. The computed profile has the shape of the
               measurement and sits a few percent below it on the shipped mesh; refining the mesh in the
               gap closes most of that, as the table below shows.
       style=width:80%;margin-left:auto;margin-right:auto;

On the shipped mesh the computed flux runs 6 to 9% below the measurements over positions 1 to 6, and the
whole profile moves up towards them under refinement. Position 7 is the exception, and it is not really
comparable: its block runs from 55 mm to the top of the plate rather than being centred on 60 mm, and it
sits exactly where the flux is turning out of the plate and into the flange, so it averages over the
steepest part of the profile.

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

The benchmark also specifies a 3000 AT case, which saturates the steel much harder; `ampere_turns=3000` on
the command line switches to it. Only the 1000 AT case is exercised by the test suite, and only the 1000 AT
measurements are quoted above, because those are the ones printed in the problem statement.

## References id=references

- T. Nakata, N. Takahashi and K. Fujiwara, *Summary of results for benchmark problem 13 (3-D non-linear
  magnetostatic model)*, COMPEL, 11, 3, 345 (1992).
- [Problem statement](https://www.compumag.org/wp/wp-content/uploads/2018/06/problem13.pdf) - the geometry,
  the B-H curve, the measured flux densities and the reporting tables used above.
