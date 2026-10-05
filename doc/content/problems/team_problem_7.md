# TEAM Problem 7: Asymmetrical Conductor with a Hole

[TEAM](https://www.compumag.org/wp/team/) Problem 7 is a time-harmonic eddy current benchmark from the
COMPUMAG workshop series. A racetrack coil carrying a sinusoidal current sits 30 mm above a thick aluminium plate
with an off-centre square hole. The benchmark compares measurements with computed values of the flux density
above the plate and of the eddy current density on its two faces, at 50 Hz and 200 Hz.

The problem has no symmetry, so the whole of it has to be modelled. The hole makes the conductor multiply
connected. The surrounding free space has no conductivity, so the curl-curl operator is singular there and has
to be gauged. This example uses a tree-cotree gauge restricted to the non-conducting region to do that.

This example uses the MFEM backend of MOOSE, so it requires an application built with MFEM support. It is
skipped by `./run_tests` when that capability is absent.

## Problem Specification id=specification

| Quantity | Value |
| :- | :- |
| Plate | 294 mm by 294 mm by 19 mm, occupying $0 \le x, y \le 294$ mm, $0 \le z \le 19$ mm |
| Hole | 108 mm by 108 mm, through the plate, occupying $18 \le x, y \le 126$ mm |
| Plate conductivity | $3.526 \times 10^{7}$ S/m |
| Coil | a winding 25 mm wide and 100 mm tall, occupying $49 \le z \le 149$ mm, 30 mm above the plate |
| Coil outline | 200 mm square outside with R50 corners, 150 mm inside with R25 corners, centred on $(x, y) = (194, 100)$ mm |
| Excitation | 2742 ampere-turns, sinusoidal, peaking at $\omega t = 0$ |
| Frequencies | 50 Hz and 200 Hz; the skin depth in the plate is 12 mm and 6 mm |
| Relative permeability | 1.0, everywhere |

The benchmark measures two quantities, at $\omega t = 0°$ and $\omega t = 90°$:

- $B_z$ along the lines A1-B1 ($y = 72$ mm) and A2-B2 ($y = 144$ mm), both at $z = 34$ mm, 15 mm above the
  plate. There are 17 points at $x = 0$ to 288 mm, spaced 18 mm apart.
- $J_y$ along $y = 72$ mm on the top (A3-B3, $z = 19$ mm) and bottom (A4-B4, $z = 0$) faces of the plate, at
  $x = 0$, 18 and 126 to 288 mm. The line crosses the hole in between.

## Mesh id=mesh

Every part of the model is a prism standing on the $x$-$y$ plane: the plate, the hole through it, and the coil.
So `team7_mesh.py` triangulates the plan once and extrudes it through a stack of slabs, giving wedge elements
throughout. That puts layer boundaries exactly on the faces of the plate and the coil, and on the $z = 34$ mm
plane where the flux density is measured. It also lets the plate be resolved through its thickness, where the
skin depth sets the scale, without paying for the same resolution in plan.

The free space box is the outer boundary of the mesh the problem statement recommends. It runs from
$-1353$ to 1647 mm in $x$ and $y$, and from $-300$ to 449 mm in $z$. The supplied mesh, `team7.e`, has 41976
wedges:

- **In plan:** 12 mm elements over the plate and the coil, grading out to 500 mm at the outer boundary.
- **In $z$:** 6 layers through the plate, 2 in each of the two 15 mm gaps between the plate and the coil, 4
  through the coil, and 4 graded layers each above and below.

There are three blocks, `Plate`, `Coil` and `Air`; the hole is part of `Air`. There is one sideset, `far_field`,
the outside of the free space box.

The mesh is built with the [gmsh](https://gmsh.info) Python API and converted to Exodus by a `--mesh-only` pass
of a MOOSE executable, so that the block and sideset names reach MFEM:

```bash
cd test/tests/electromagnetics/team/problem_7
./team7_mesh.py --exe /path/to/moose_test-opt    # or set MOOSE_DIR
```

Run it with `--help` for the sizing controls. The defaults reproduce the supplied mesh.

## Formulation id=formulation

The fields vary as $e^{i\omega t}$, and are solved for as a complex magnetic vector potential
$\vec{A}$, with $\vec{B} = \nabla \times \vec{A}$, in $H(\text{curl})$. Inside the plate the electric scalar
potential is absorbed into $\vec{A}$ (the $\vec{A}^*$, or modified vector potential, formulation). The electric
field there is then $\vec{E} = -i\omega\vec{A}$ and the eddy current density is
$\vec{J} = -i\omega\sigma\vec{A}$. The weak form is

!equation
\left( \nu \nabla \times \vec{A}, \nabla \times \vec{A}' \right)
+ \left( i\omega\sigma \vec{A}, \vec{A}' \right)_\text{Plate}
= \left( \vec{J}_s, \vec{A}' \right)

where $\nu$ is the reluctivity of free space everywhere. The coil is a stranded winding: it carries the source
current $\vec{J}_s$ but no eddy current, so like the air it has no conductivity. The hole needs no special
treatment. It is free space, and the $\vec{A}^*$ formulation handles the multiply connected plate around it as it
stands.

The outer faces carry $\vec{n} \times \vec{A} = 0$, which is $\vec{B} \cdot \vec{n} = 0$, as the benchmark
asks.

### Gauging the Non-conducting Region id=gauge

Inside the plate the conductivity term removes the gradients from the null space of the curl-curl operator. In
the air and the coil nothing does, and the system is singular. The
`MFEMComplexTreeCotreeGaugeEssentialConstraint` is restricted to the `Air` and `Coil` blocks. It fixes the real
and imaginary parts of $\vec{A}$ to zero on the edges of a spanning tree of those blocks, which removes exactly
those gradients without touching the plate. The free space region threads the hole and so is not simply
connected, but the cotree handles that. The plate's edges seed the spanning forest and are never gauged.

### A Divergence-free Source id=source

A gauged system only reproduces the physical solution if its right hand side is orthogonal to the gradients it
removes, that is, if $\vec{J}_s$ is discretely divergence free. The benchmark's winding current has a uniform
magnitude of $2742 / (0.025 \times 0.1) = 1.097 \times 10^{6}$ A/m$^2$ and follows the coil around. It can be
written down analytically, but projecting it onto the mesh is not enough. The meshed coil is faceted where its
corners are curved, and the projected current leaks out across the facets.

So, as in TEAM Problem 13, the source is built in two parts:

!equation
\vec{J}_s = \vec{J}_\text{applied}\big|_\text{Coil} - \nabla\phi,
\qquad
\left( \nabla\phi, \nabla\psi \right) = \left( \vec{J}_\text{applied}, \nabla\psi \right)_\text{Coil}
\quad \forall \psi

This makes $(\vec{J}_s, \nabla\psi) = 0$ for every $\psi$. It is the smallest correction that makes the source
discretely divergence free. The two equations are solved in one system rather than with a sub-app, because a
complex problem cannot couple to real variables. $\phi$ is a second complex variable, `correction_potential`,
and the applied current is a complex auxiliary variable set by an initial condition. The $\phi$ equation does not
depend on $\vec{A}$, so the coupled system is block triangular.

The signs of the coupling kernels in the input reflect a convention of the MFEM backend. A coupled variable that
is not solved for is moved to the right hand side, which negates the kernel's contribution.

### Solving id=solving

The complex system is assembled as a real system of twice the size and factorised directly with MUMPS. For the
supplied mesh that is 219234 unknowns, including $\phi$.

!listing canary/test/tests/electromagnetics/team/problem_7/team7.i

## Reading the Benchmark Tables id=tables

The coil current is real in the solve, so it peaks at $\omega t = 0$, as the benchmark specifies. The field at an
instant is the real part of the phasor times $e^{i\omega t}$, so the two measured columns correspond to

!equation
F(\omega t = 0°) = \text{Re}\,\hat{F}, \qquad F(\omega t = 90°) = \text{Re}\,(i\hat{F}) = -\text{Im}\,\hat{F}.

The compumag.org copy of the problem statement labels its two $J_y$ tables A3-B3 ($z = 19$ mm) and A4-B4
($z = 0$). Those labels are the wrong way round. The table it labels A4-B4 has much the larger currents: three
times larger at $x = 288$ mm at 200 Hz, and twenty times larger mid-plate. A plate three skin depths thick
carries its larger currents on the face nearer the coil, which is the top. Published comparisons use that table
for $z = 19$ mm too; for example the measurements plotted in Fig. 8 of Kapidani et al. are its values. Here each
table is compared with the face its values physically belong to.

## Files id=files

| File | Role |
| :- | :- |
| `team7.i` | The time-harmonic solve. This is the input to run. |
| `team7.e` | The mesh, 41976 wedges. |
| `team7_mesh.py` | Builds `team7.e` with gmsh. |
| `plot_team7.py` | Compares the sampled fields with the measurements and draws the figures below. |
| `tests` | Test harness specification. |
| `gold/` | Reference output for both frequencies. |

All of them live in `test/tests/electromagnetics/team/problem_7`.

## Running id=running

```bash
<path to a MOOSE app built with MFEM support> -i team7.i                   # 50 Hz
<path to a MOOSE app built with MFEM support> -i team7.i frequency=200 \
    Outputs/ReportedPostprocessors/file_base=OutputData/TEAM7CSV_200Hz
```

or, from the repository root, through the test harness. The 200 Hz case is marked heavy:

```bash
./run_tests --re TEAM_Problem_7 --all-tests
```

Each run takes about 3 minutes on one process or 1 minute on four. The direct factorisation needs 6 to 7 GB of
memory in total, whatever the number of processes. Results do not depend on the number of processes: the sampled
values agree to 1 part in $10^{12}$ between runs on 1, 2 and 4.

Results are written to `OutputData/`:

- `TEAM7CSV_bz_a1_b1_0001.csv`, `TEAM7CSV_bz_a2_b2_0001.csv`, `TEAM7CSV_jy_a3_b3_0001.csv` and
  `TEAM7CSV_jy_a4_b4_0001.csv` hold the real and imaginary parts of every component of $\vec{B}$ or $\vec{J}$
  along the four lines, in SI units.
- `TEAM7CSV.csv` holds the time-averaged Joule loss in the plate.
- `TEAM7/` is a ParaView collection with `a_field`, `b_field`, `j_field` and the source current.

The $J_y$ points all lie on a face of the plate, and three of them on an edge of it or of the hole, where the
current density is discontinuous. Each is moved 1 µm into the plate so that it is sampled from the plate's side.

## Results id=results

`plot_team7.py` reads the sampled values from both runs, converts them to the benchmark's units and phases,
prints them against the measurements, and draws the figures:

```bash
cd test/tests/electromagnetics/team/problem_7
./plot_team7.py            # --no-plot to print the comparison only
```

It needs only `matplotlib`, and writes to `doc/content/media/`.

The benchmark tabulates each quantity at two instants, which together give its amplitude and phase. Published
comparisons, for example those of COMSOL and of Kapidani et al., plot a single curve per line instead:
$\text{sign}(\text{Re}\,F)\,|F|$, the magnitude carrying the sign of the $\omega t = 0°$ value. That form hides
almost all of the phase. The first figure uses it, so that these results can be read against those comparisons.
The other two separate the instants, and so test the phase as well.

!media media/team7_signed.png
       id=team7-signed
       caption=Every line as $\text{sign}(\text{Re}\,F)\,|F|$, the form published comparisons use, at 50 Hz and
               200 Hz. Hollow markers are the 50 Hz bottom-face measurements, which are not physically consistent;
               see below.
       style=width:100%;margin-left:auto;margin-right:auto;

!media media/team7_bz.png
       id=team7-bz
       caption=$B_z$ 15 mm above the plate along A1-B1 and A2-B2, at 50 Hz (left) and 200 Hz (right), at the two
               instants the benchmark reports.
       style=width:100%;margin-left:auto;margin-right:auto;

!media media/team7_jy.png
       id=team7-jy
       caption=$J_y$ along $y = 72$ mm on the top (A3-B3) and bottom (A4-B4) faces of the plate, at the two
               instants the benchmark reports. The curve is broken where the line crosses the hole.
       style=width:100%;margin-left:auto;margin-right:auto;

Most of these profiles pass through zero, so the table below gives the RMS difference over each line as a
percentage of the largest measured magnitude on that line, over both instants. "Signed" is the signed magnitude
of the first figure, and 0° and 90° are the two instants.

| Line | 50 Hz, signed | 50 Hz, 0° | 50 Hz, 90° | 200 Hz, signed | 200 Hz, 0° | 200 Hz, 90° |
| :- | -: | -: | -: | -: | -: | -: |
| A1-B1, $B_z$ | 3.1% | 2.8% | 0.6% | 3.1% | 3.0% | 2.9% |
| A2-B2, $B_z$ | 2.1% | 2.3% | 0.8% | 2.5% | 2.3% | 4.3% |
| A3-B3, $J_y$ (top) | 2.5% | 12.7% | 11.4% | 2.9% | 12.2% | 10.2% |
| A4-B4, $J_y$ (bottom) | -- | -- | -- | 12.4% | 15.4% | 15.1% |

The 50 Hz bottom-face row is left out: those measurements are not physically consistent, as explained below. By
the measure published comparisons use, the flux density and the top-face current agree to 2% to 3%.

### Flux density id=results-b

The flux density agrees closely at both frequencies and both instants.

The two lines differ from each other for physical reasons, and the model reproduces the difference. The coil
alone, solved at 0.001 Hz where there are no eddy currents, gives nearly the same profile on both, about 100 G
inside the coil. The plate's eddy currents cut that roughly in half at 50 Hz. A2-B2, at $y = 144$ mm, runs over
solid plate. A1-B1, at $y = 72$ mm, runs over the hole, and the current circulating around the hole changes the
field above it. Over the hole's left half the return flux on A1-B1 is $-18$ G to $-22$ G, against $-9$ G to
$-15$ G on A2-B2, and above its right edge, at $x = 126$ mm, the field is 78 G against 60 G.

At 200 Hz the $\omega t = 90°$ component is about 6 G in the coil interior against about 4 G measured. This comes
from a nearly uniform phase difference of 2.5° to 4° along both lines, against less than 1° at 50 Hz. The
$\omega t = 90°$ column is small, and a 3° shift of a 63 G field is the whole of the difference. The amplitude
agrees to about 1%. Changing the plate's conductivity by 10% either way moves the phase by only 0.3°, so an
uncertain conductivity does not explain it. The cause has not been found.

### Eddy current density on the top face id=results-top

On the top face the computed magnitude of $J_y$ is within 8% of the measured one at every point, at both
frequencies, and the signed magnitude agrees to under 3% RMS. Away from the hole the phase agrees to within 8°.

Near the hole it does not, and that drives the larger differences at the two instants. At 50 Hz the phase differs
by 18° to 27° at $x = 0$ and between $x = 126$ and 180 mm. At 200 Hz the difference is at the hole's edges: 18°
at $x = 18$ mm and 27° at $x = 126$ mm. At those two points the model matches the measured bottom-face current to
within 5° in phase and 4% in magnitude. The measurements therefore have the top face lagging the bottom face at
the hole's edges. Current diffuses into the plate from the top, the face nearer the coil, so the top face should
lead. That points to an error in the top-face measurements near the hole rather than in the model, but the
problem statement gives no detail of how the currents were measured with which to confirm it.

### Eddy current density on the bottom face id=results-bottom

At 200 Hz the bottom face agrees at the edges of the hole and mid-plate, where the current has almost died away.
The disagreement is at $x = 270$ and 288 mm, within a few millimetres of the plate's outer edge, where the current
turns to run down the side face.

At 50 Hz the bottom-face measurements are not physically consistent, and are not used in the comparison. The
ratio of the bottom-face current to the top-face current at the same point shows why. Through a plate
$d = 19$ mm thick, current driven from one face falls off and lags with depth. For a slab it would be reduced by
$|1/\cosh(kd)|$, with $k = (1 + i)/\delta$ and skin depth $\delta$, which at 50 Hz is a factor of 0.43 with a
lag of 91°. Mid-plate, between $x = 216$ and 252 mm, the model gives 0.42 to 0.50 with a lag of 68° to 76°. Over
the same stretch the measured 50 Hz bottom current is 0.8 to 0.9 of the top face's and leads it by 54° to 73°: it
is barely attenuated, and ahead in phase rather than behind. No sign convention for the table's columns
reconciles it. At 200 Hz, by contrast, the model and the measurements both reduce the current mid-plate to under
a tenth, as the slab does.

Oddly, its magnitudes match the model's almost value for value if moved one position along the line (18 mm),
for example 1.11 against 1.10 and 0.88 against 0.88, while its phases still do not. That suggests the table was
corrupted, perhaps in transcription, but without the original measurements it cannot be corrected. Published
comparisons do not use the bottom face at all.

### Joule loss id=results-loss

`JouleLoss`, the time-averaged loss in the plate, $\frac{1}{2}\int |\vec{J}|^2/\sigma$, is 4.90 W at 50 Hz and
10.45 W at 200 Hz. The benchmark does not ask for it, but it summarises the whole eddy current distribution in
one number. Khebir et al. report 9.60 W at 200 Hz with second order elements on 424000 tetrahedra, and 9.66 W with
first order elements on 755000.

## Convergence and Robustness id=convergence

None of the differences from the measurements described above depends on how the problem is discretised or
truncated. Each check below was run against the supplied setup:

- **Layers through the plate, 6 to 12, at 50 Hz.** The sampled values move by at most 0.5% of the peak on any
  line, so the plate is resolved through its thickness.
- **Plan refined to 6 mm along the edges of the plate and the hole, at both frequencies**
  (`--edge-size 6 --edge-band 20`, 73612 wedges). The overall agreement with the measurements is unchanged, as the
  RMS differences show. Individual values sampled directly above or on those edges move by up to 8% of the peak.
  For $B_z$ that is up to 5 G, where the field changes by up to 70 G over a few centimetres. Lowest-order
  Raviart-Thomas elements represent $B$ only to first order pointwise, so those points are the most mesh-sensitive
  in the comparison. The Joule loss falls 3% at 50 Hz and 5% at 200 Hz, to 4.75 W and 9.94 W.
- **Each wedge split into three tetrahedra, at 50 Hz.** The RMS differences change by at most 1.4% of the peak,
  so the results do not depend on the element type.
- **The free space box enlarged six-fold in $z$ and doubled in plan, at both frequencies.** The RMS differences
  change by at most 1.2% of the peak, so the truncation boundary the benchmark specifies is far enough away.
- **The gauge replaced by a small conductivity in the air and the coil, at 50 Hz.** A conductivity of 1 S/m,
  $3 \times 10^{-8}$ of the plate's, also removes the null space. The results agree with the gauged ones to
  round-off, so the gauge fixes the potential without affecting the physics.
- **The plate's conductivity changed by 10% either way, at 200 Hz.** The phase of $B_z$ moves by 0.3° and its
  amplitude by under 1%.

The edge-refined mesh needs 11.7 GB for the factorisation, against 7 GB for the supplied one, which is why it is
not the default. Second order elements would be the more economical route to a converged answer, but the
tree-cotree gauge requires a first order space.

## References id=references

- K. Fujiwara and T. Nakata, *Results for benchmark problem 7 (asymmetrical conductor with a hole)*, COMPEL
  9(3), 137-154, 1990.
- [Problem statement](https://www.compumag.org/wp/wp-content/uploads/2018/06/problem7.pdf), with the measured
  values.
- B. Kapidani, M. Merkel, S. Schöps and R. Vázquez, *Arbitrary order spline representation of cohomology
  generators for isogeometric analysis of eddy current problems*,
  [arXiv:2310.14313](https://arxiv.org/abs/2310.14313).
- A. Khebir, P. Dłotko, A. Kouki, B. Kapidani and R. Specogna, *$T$-$\Omega$ formulation with higher order
  hierarchical basis functions for non simply connected conductors*,
  [arXiv:1704.03694](https://arxiv.org/abs/1704.03694).
