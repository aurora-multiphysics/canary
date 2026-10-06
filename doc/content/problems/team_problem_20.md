# TEAM Problem 20: The Static Force Problem

[TEAM](https://www.compumag.org/wp/team/) Problem 20 is a non-linear magnetostatic benchmark from the
COMPUMAG workshop series, set to test the calculation of magnetic force. A coil wound round a steel centre
pole drives flux up the pole, across into the top plates of a steel yoke, down the yoke's two legs and back
across a 1.5 mm air gap into the foot of the pole. The quantity asked for is the force on the pole, at four
excitations from 1000 AT, where the steel is well below the knee of its B-H curve, to 5000 AT, where the pole
carries over 2 T.

Force is a harder quantity to get right than flux density. It depends on the square of the field, it is
concentrated in the air gaps where the field changes fastest, and the obvious ways of computing it from a
finite element solution - integrating the Maxwell stress over a surface, or differentiating the energy with
respect to a displacement - are sensitive to exactly where the surface cuts the mesh or to remeshing noise.
This example takes the force as a volume integral over a thin shell of air around the pole, which is
insensitive to both, and checks it two ways.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

All dimensions are from Fig. 1 of the problem statement. The origin is at the centre of the underside of
the yoke, with $z$ up the pole.

| Quantity | Value |
| :- | :- |
| Centre pole | 25 mm ($x$) by 10 mm ($y$) by 98.5 mm ($z$), from $z = 26.5$ to $125$ mm |
| Yoke | a rectangular frame 127 mm ($x$) by 25 mm ($y$) by 150 mm ($z$), 25 mm thick |
| Yoke window | $\lvert x \rvert \le 38.5$ mm, $25 \le z \le 125$ mm |
| Opening in the top plates | 27 mm across, $\lvert x \rvert \le 13.5$ mm |
| Gaps | 1.5 mm under the foot of the pole; 1 mm from the top corner of the pole to the top plates |
| Coil | rounded square, 39 mm inner and 75 mm outer, R5 and R23 corners, 18 mm thick |
| Coil height | 96.6 mm, 1.7 mm clear of the yoke top and bottom |
| Excitation | 1000, 3000, 4500 and 5000 AT |

The top of the pole is level with the underside of the top plates, so it does not overlap them: the flux
leaving the pole has to turn and cross 1 mm sideways to reach them, and much of it fringes. The foot of the
pole faces the bottom of the yoke squarely across its 1.5 mm gap. The yoke is 25 mm deep in $y$ and the
coil 75 mm across, so the coil passes through the window on either side of the pole and closes round
outside the yoke.

### The B-H Curve id=bh-curve

The steel is specified by the 38 tabulated $(B, H)$ points of Table 1 of the problem statement, up to
$B = 2.3$ T and $H = 135000$ A/m, and above that by Eq. (1),

!equation
B = \mu_0 H + M_s,

with the saturation magnetisation quoted as $M_s = 2.16$ T. Table 6 asks which interpolation was used
between the tabulated points and lists splines first, so `bh_curve.py` fits a monotone cubic
(Fritsch-Carlson) through all 38 points and joins it at 2.3 T to Eq. (1). It reproduces every tabulated
point exactly and $H(B)$ is monotone throughout.

One liberty is taken with the specification, and `./bh_curve.py --check` reports it. With $M_s = 2.16$ T,
Eq. (1) does not pass through the last tabulated point: at $H = 135000$ A/m it gives 2.3296 T where the table
says 2.3, a 30 mT step. $M_s$ is instead taken as $2.3 - \mu_0 \times 135000 = 2.1304$ T, which continues the
line from the last point, and the spline's end slope is clamped to $1/\mu_0$ so that the curve is
continuously differentiable across the join.

The resulting reluctivity $\nu = H/B$ runs from 2700 m/H as $B \rightarrow 0$, down to 504 m/H at
$B = 0.75$ T, and up to 58700 m/H at $B = 2.3$ T.

!media media/team20_bh_curve.png
       id=bh-figure
       caption=The B-H curve as the input uses it: the monotone spline through the 38 tabulated points,
               continued above 2.3 T by Eq. (1) of the problem statement with the saturation magnetisation
               adjusted to meet the last point.
       style=width:80%;margin-left:auto;margin-right:auto;

## Symmetry and Mesh id=mesh

The model is symmetric under reflection in $x = 0$ and in $y = 0$, so only the quarter with $x \ge 0$ and
$y \ge 0$ is meshed, which is the region Fig. 2 of the problem statement allows. It is not symmetric in $z$ -
the top plates have an opening and the bottom plate does not - so no further reduction is possible. The free
space box runs 250 mm out in $x$ and $y$ and 200 mm above and below the yoke.

`team20_mesh.py` builds the mesh with [gmsh](https://gmsh.info). The coil is a filleted box with a smaller
filleted box cut out of it, which reproduces the rounded square of the problem statement exactly, cut down to
the quarter. Element size is set by a small number of axis-aligned box fields: finest in the gap under the
pole and at its top corner (0.3 mm), then through the pole and the air round it (1 mm), the yoke (3 mm), the
coil (6 mm), and coarsening out to the truncation boundary. The shipped mesh has 141k tetrahedra.

The blocks are `Pole`, `Yoke`, `Coil` and `Air`, together with two measurement blocks, `PoleZ075` and
`YokeZ075`, which are 5 mm bands of the pole and of the yoke leg centred on $z = 75$ mm, where Table 3 of the
problem statement asks for the average flux density. Both are steel and have to be listed in the input's
`steel_blocks`. Three sidesets are named: `symmetry_x` and `symmetry_y`, the two symmetry planes, and
`far_field`, the outer faces of the box.

MFEM reads named blocks and sidesets from Exodus files only, so the script runs the mesh back through MOOSE
with `--mesh-only` to convert it. That step needs a MOOSE executable; the rest needs only gmsh.

## Formulation id=formulation

The field is solved for exactly as in the TEAM Problem 13 example: the magnetic vector potential
$\vec{A}$, with $\vec{B} = \nabla \times \vec{A}$, in $H(\text{curl})$,

!equation
\left( \nu\!\left(\left|\nabla \times \vec{A}\right|\right) \nabla \times \vec{A}, \nabla \times \vec{A}' \right)
= \left( \vec{J}_s, \vec{A}' \right),

with $\nu$ from the [B-H curve](#bh-curve) in the steel and $1/\mu_0$ elsewhere, solved by Newton's method
with `MFEMNLCurlCurlKernel` supplying the Jacobian from $\nu$ and
$B\,\partial\nu/\partial B = H'(B) - H(B)/B$. A small mass term, $10^2 (\vec{A}, \vec{A}')$, removes the
null space of the curl-curl operator, and the preconditioner is `MFEMHypreAMS` with `singular = false`.

The coil current is a uniform winding density, $NI/A$ over the 18 mm by 96.6 mm cross-section, running
anticlockwise seen from above so that it drives flux up the pole. As in TEAM Problem 13, the sub-app
`team20_coil_source.i` projects it onto the mesh and removes the part of it that leaks across the faceted
coil boundary, so that the source is discretely divergence free.

!listing canary/test/tests/electromagnetics/team/problem_20/team20_coil_source.i

### Boundary Conditions id=boundary-conditions

The planes $x = 0$ and $y = 0$ both contain the axis of the coil, and the field lies in them: $B_x = 0$ on
the first and $B_y = 0$ on the second. That is $\vec{B} \cdot \vec{n} = 0$, which is
$\vec{n} \times \vec{A} = 0$ - the same condition as on the far-field faces, where it holds the flux inside
the box. All three sidesets carry it.

Unlike the symmetry plane of TEAM Problem 13, these planes are crossed by the coil current, which flows
normally through them. The coil source sub-app pins its correction potential on all three sidesets for
that reason: the gradients left in the main problem's null space are those of potentials vanishing on the
whole tangential boundary, so the current only has to be divergence free in the interior and may cross the
symmetry planes freely.

### The Force id=force

The force on a body in air is the Maxwell stress integrated over any closed surface in the air around it,

!equation
\vec{F} = \oint_S \mathbf{T} \cdot \vec{n} \, dS, \qquad
\mathbf{T} = \frac{1}{\mu_0} \left( \vec{B} \vec{B} - \tfrac{1}{2} \left|\vec{B}\right|^2 \mathbf{I} \right).

Taken over one surface, that is sensitive to exactly where the surface cuts the mesh: the field from a first
order solution is constant on each element and jumps between them, and the stress is quadratic in it. So it is
averaged over a family of surfaces instead. With a weight $g$ that is 1 on the pole and falls to 0 across a
thin shell of air round it, the divergence theorem turns the surface integral into

!equation
\vec{F} = -\int_\text{shell} \mathbf{T} \cdot \nabla g \, dV,

which is the surface integral averaged over the level sets of $g$. When the shell is one element thick, this
is Coulomb's virtual work method: $-\nabla g$ is the derivative of the mesh coordinates with respect to a
rigid displacement of the pole. Its $z$ component, the one the benchmark asks for, is

!equation
F_z = \int \frac{\left|\vec{B}\right|^2}{2\mu_0} \frac{\partial g}{\partial z} \, dV
    - \int \frac{B_z}{\mu_0} \, \vec{B} \cdot \nabla g \, dV,

which the input evaluates as two `MFEMVectorFEInnerProductIntegralPostprocessor`s - the inner product of
$\nabla g$ with $\hat{z}$ and with $\vec{B}$, scaled by the coefficients above - combined by a
`ParsedPostprocessor` and multiplied by four for the four quarters of the pole. A positive force is upwards.

$g$ is the linear ramp $\max(0, 1 - d/\delta)$ in the distance $d$ from the pole, interpolated onto the
lowest order $H^1$ space so that it is exactly 1 at every node of the pole and exactly 0 at every node of the
yoke and the coil. The shell thickness $\delta$ has to stay under the 1 mm between the top corner of the pole
and the top plates, and the input uses 0.9 mm. Within that, the answer barely depends on it - on the coarsest
mesh of the [convergence study](#convergence):

| $\delta$ | 1000 AT | 5000 AT |
| :- | -: | -: |
| 0.01 mm (one element: Coulomb's method) | -8.262 N | -84.64 N |
| 0.2 mm | -8.262 N | -84.64 N |
| 0.5 mm | -8.256 N | -84.56 N |
| 0.9 mm | -8.239 N | -84.30 N |

#### Checking the force id=force-check

The same machinery was pointed at the coil, which also sits in air, with a shell round it and the coil's
Lorentz force, $\int (\vec{J} \times \vec{B}) \cdot \hat{z} \, dV$, computed directly alongside. The two
agreed to 1.7% at 1000 AT (-0.353 N against -0.359 N) and 1.3% at 5000 AT (-5.24 N against -5.31 N) on the
coarsest mesh. That checks the signs, the factors and the symmetry factor of four independently of the
benchmark measurements. The coil force is small and set by the leakage field, which is why it is the more
demanding of the two.

#### Execution order id=execution-order

The force postprocessors are executed on `FINAL` rather than at the end of the time step, and that matters.
MFEM aux kernels and postprocessors are ordered by the variables they name as variables. A coefficient that
reads an auxiliary variable - here `stress_bz`, which reads `b_z` - is not seen as a dependency, so a
postprocessor using it can run before the aux kernel that fills the variable. In development the tension term
did exactly that and read $B_z = 0$, which gave +8.38 N instead of -8.24 N. Every aux kernel runs at the end
of the time step and before `FINAL`, so the force is always taken from current fields.

The point samplers for Tables 2 and 4 avoid the same trap a different way: they are
`MFEMVariablePointValueSampler`s, which name `b_z` and `b_x` as variables, rather than coefficient samplers.

### Newton Damping id=damping

`line_search = backtracking` is set, as it is in TEAM Problem 13, but this problem does not need it. Every
excitation converges in four to ten Newton iterations with or without it, at first and second order, and
the results agree to $10^{-10}$. Without damping the first update overshoots - the residual rises by a factor
of three - and Newton then converges quadratically; with it, that step is shortened and the iteration count
moves by one or two either way. It is kept because it is cheap insurance for finer meshes or further
saturation, where TEAM Problem 13 showed an undamped iteration can settle into a limit cycle.

!listing canary/test/tests/electromagnetics/team/problem_20/team20.i

## Files id=files

| File | Role |
| :- | :- |
| `team20.i` | The magnetostatic solve. This is the input to run. |
| `team20_coil_source.i` | Sub-app: the divergence-free coil current. Launched automatically. |
| `team20.e` | The mesh. |
| `team20_mesh.py` | Builds `team20.e` with gmsh. |
| `bh_curve.py` | Turns the benchmark B-H curve into the expressions in `team20.i`, and draws it. |
| `plot_team20_force.py` | Draws the force figure below from the postprocessor CSVs. |
| `tests` | Test harness specification. |
| `gold/` | Reference postprocessor output. |

All of them live in `test/tests/electromagnetics/team/problem_20`.

## Running id=running

The main app launches the coil source sub-app itself, so only one input needs to be given:

```bash
<path to a MOOSE app built with MFEM support> -i team20.i
```

or, from the repository root, through the test harness, which runs all four excitations:

```bash
./run_tests --re TEAM_Problem_20
```

The other excitations are a command line argument:

```bash
<app> -i team20.i ampere_turns=3000 \
      Outputs/ReportedPostprocessors/file_base=OutputData/TEAM20CSV_3000AT
```

Results are written to `OutputData/`. `TEAM20CSV.csv` holds the force and the Table 3 averages,
`TEAM20CSV_gap_points_0001.csv` the Table 2 points, `TEAM20CSV_line_ab_0001.csv` and
`TEAM20CSV_line_cd_0001.csv` the two lines of Table 4, and `TEAM20/` a ParaView collection. Each case takes
about 50 s on one process, and about 10 s on eight. The figure below is drawn by `./plot_team20_force.py` once
all four have run.

To rebuild the mesh, which needs gmsh and a MOOSE executable:

```bash
cd test/tests/electromagnetics/team/problem_20
./team20_mesh.py --exe /path/to/app-opt        # or set MOOSE_DIR
```

`--help` lists the element size controls, and `mesh_file=` on the command line runs the input on another
mesh.

## Results id=results

### The force id=force-results

The problem statement publishes no measurements; they are in the summary of results listed under
[References](#references), which is not openly available. The values here are read off Fig. 6 of Hafla et
al., which plots them, and are good to about half a newton.

| Excitation | Measured | Computed | Error |
| :- | -: | -: | -: |
| 1000 AT | 8.2 N | 8.24 N | +0.5% |
| 3000 AT | 54.5 N | 56.58 N | +3.8% |
| 4500 AT | 75.1 N | 77.34 N | +3.0% |
| 5000 AT | 80.2 N | 82.69 N | +3.1% |

The force is downwards in every case, pulling the pole towards the yoke beneath it; the table gives its size.
It grows faster than the current up to 3000 AT - seven times the force for three times the current - and then
much more slowly, as the pole saturates:
the average flux density in it is 0.70 T at 1000 AT and 1.77 T at 3000 AT, but only 2.04 T at 5000 AT. The
computed force is a few percent high at the three saturated excitations, and the [convergence
study](#convergence) shows most of that is discretisation error, which second order elements remove.

!media media/team20_force.png
       id=force-figure
       caption=The downward force on the centre pole at the four excitations of the benchmark, computed on the
               shipped mesh with first order elements, against the measurements.
       style=width:80%;margin-left:auto;margin-right:auto;

### Flux densities id=flux-results

The problem statement asks for three further sets, none with measurements, for comparison between methods.
All are reported with the sign the field has: the coil drives flux up the pole and it returns down the yoke.

Table 2 is $B_z$ in the gap under the pole, at $P_1$ on its axis and at $P_2$ under its corner edge. Table 3
is the average of $B_z$ over the section of the pole ($\alpha$-$\beta$) and of the yoke leg ($\gamma$-$\delta$)
at $z = 75$ mm, reported as the average over the 5 mm bands centred there.

| Quantity | Position (mm) | 1000 AT | 3000 AT | 4500 AT | 5000 AT |
| :- | :- | -: | -: | -: | -: |
| $B_z$ at $P_1$ | (0, 0, 25.75) | 0.321 | 0.841 | 0.996 | 1.034 |
| $B_z$ at $P_2$ | (12.5, 5, 25.75) | 0.213 | 0.539 | 0.619 | 0.638 |
| average $B_z$, pole | $z = 75$ | 0.700 | 1.774 | 1.995 | 2.040 |
| average $B_z$, yoke leg | $z = 75$ | -0.133 | -0.359 | -0.439 | -0.460 |

$P_2$ sits directly under the corner edge of the pole, where the field is singular, so it is the least well
converged of these: it moves by several percent between meshes where $P_1$ moves by a fraction of one. Each
yoke leg returns half the pole's flux through a section two and a half times larger than the pole's, so it
runs at about a fifth of the pole's flux density, less what leaks across the window.

Table 4 is $B_x$ in the air at the top corner of the pole, where the flux crosses to the top plates: along
a-b, which runs up past the corner at $x = 13$ mm, $y = 4$ mm, and along c-d, which runs out across it at
$x = 13$ mm, $z = 125$ mm. Point 6 of a-b and point 16 of c-d are the same point.

| No. | $x$ (mm) | $y$ (mm) | $z$ (mm) | 1000 AT | 3000 AT | 4500 AT | 5000 AT |
| :- | -: | -: | -: | -: | -: | -: | -: |
| 1 (a) | 13 | 4 | 120 | 0.108 | 0.289 | 0.351 | 0.367 |
| 2 | 13 | 4 | 121 | 0.123 | 0.330 | 0.403 | 0.422 |
| 3 | 13 | 4 | 122 | 0.160 | 0.431 | 0.530 | 0.557 |
| 4 | 13 | 4 | 123 | 0.217 | 0.582 | 0.704 | 0.734 |
| 5 | 13 | 4 | 124 | 0.352 | 0.921 | 1.066 | 1.097 |
| 6 | 13 | 4 | 125 | 0.519 | 1.203 | 1.338 | 1.369 |
| 7 | 13 | 4 | 126 | 0.269 | 0.688 | 0.815 | 0.847 |
| 8 | 13 | 4 | 127 | 0.184 | 0.483 | 0.582 | 0.607 |
| 9 | 13 | 4 | 128 | 0.115 | 0.306 | 0.372 | 0.389 |
| 10 | 13 | 4 | 129 | 0.078 | 0.207 | 0.254 | 0.266 |
| 11 (b) | 13 | 4 | 130 | 0.064 | 0.171 | 0.210 | 0.220 |
| 12 (c) | 13 | 0 | 125 | 0.515 | 1.254 | 1.417 | 1.453 |
| 13 | 13 | 1 | 125 | 0.526 | 1.278 | 1.441 | 1.477 |
| 14 | 13 | 2 | 125 | 0.522 | 1.268 | 1.427 | 1.461 |
| 15 | 13 | 3 | 125 | 0.512 | 1.244 | 1.398 | 1.432 |
| 16 | 13 | 4 | 125 | 0.519 | 1.203 | 1.338 | 1.369 |
| 17 | 13 | 5 | 125 | 0.448 | 1.001 | 1.123 | 1.151 |
| 18 | 13 | 6 | 125 | 0.277 | 0.697 | 0.818 | 0.847 |
| 19 | 13 | 7 | 125 | 0.136 | 0.362 | 0.442 | 0.463 |
| 20 | 13 | 8 | 125 | 0.130 | 0.345 | 0.421 | 0.441 |
| 21 | 13 | 9 | 125 | 0.075 | 0.202 | 0.249 | 0.261 |
| 22 | 13 | 10 | 125 | 0.065 | 0.175 | 0.216 | 0.227 |
| 23 | 13 | 11 | 125 | 0.049 | 0.132 | 0.165 | 0.174 |
| 24 | 13 | 12 | 125 | 0.049 | 0.132 | 0.165 | 0.174 |
| 25 (d) | 13 | 13 | 125 | 0.015 | 0.042 | 0.052 | 0.055 |

$B_x$ peaks at the corner, $z = 125$ mm, and falls away either side; across the pole it holds up to its edge
at $y = 5$ mm and then drops sharply. These are values of a first order field, which is constant on each
element, so neighbouring points that fall in the same element report the same value - points 23 and 24 do.
Second order elements, or a finer mesh round the corner, resolve the profile more smoothly.

## Convergence id=convergence

The four meshes below were built by `team20_mesh.py` with every element size scaled together, the shipped
mesh lying between the first two. The force is reported as its size.

| Mesh | Elements | Order | 1000 AT | 3000 AT | 4500 AT | 5000 AT |
| :- | -: | :- | -: | -: | -: | -: |
| coarse | 59k | first | 8.24 | 57.38 | 78.79 | 84.30 |
| shipped | 141k | first | 8.24 | 56.58 | 77.34 | 82.69 |
| medium | 161k | first | 8.26 | 56.58 | 77.30 | 82.65 |
| fine | 423k | first | 8.31 | 56.24 | 76.59 | 81.86 |
| finest | 1.19M | first | 8.29 | 55.89 | 76.02 | 81.24 |
| coarse | 59k | second | 8.31 | 55.73 | 75.66 | 80.82 |
| medium | 161k | second | 8.35 | 55.69 | 75.52 | 80.67 |
| Measured | | | 8.2 | 54.5 | 75.1 | 80.2 |

At the three saturated excitations the first order force falls steadily towards the measurement as the mesh
is refined, and second order elements on the coarse mesh already do better than first order on the finest.
The two second order results agree with each other to 0.3%, which puts the converged answer within about 2% of
the measurement at 3000 AT and within 1% at 4500 and 5000 AT - inside the uncertainty of reading the
measurements off a plot. At 1000 AT every result lies within 2% of the measurement and of each other.

The test runs the shipped mesh at first order, which takes under a minute a case. Second order on the shipped
mesh takes two to three minutes a case on eight processes:

```bash
<app> -i team20.i fe_order=SECOND flux_order=FIRST
```

## References id=references

- N. Takahashi, T. Nakata and H. Morishige, *Summary of results for problem 20 (3-D static force problem)*,
  COMPEL, 14, 2/3, 57 (1995) - the measurements.
- [Problem statement](https://www.compumag.org/jsite/images/stories/TEAM/problem20.pdf) - the geometry, the
  B-H curve and the reporting tables used above.
- W. Hafla, A. Buchau and W. M. Rucker, *Solution of magnetostatic field problems with the integral equation
  method*, [ICS Newsletter, 13, 2 (2006)](https://www.compumag.org/jsite/images/stories/newsletter/ICS-06-13-2-Hafla.pdf)
  - Fig. 6 plots the measured force used above.
