# MFEM Example 27

A MOOSE reproduction of MFEM's parallel
[Example 27](https://github.com/mfem/mfem/blob/master/examples/ex27p.cpp) (`ex27p`). It solves the
Laplace problem

!equation
-\nabla^2 u = 0

on a rectangle with two circular holes, with a different boundary condition on each part of the
boundary:

| Boundary | Condition |
| :- | :- |
| Lower edge (attribute 1) | Neumann, $\mathbf{n} \cdot \nabla u = 1$ |
| Upper edge (attribute 2) | Robin, $\mathbf{n} \cdot \nabla u + u = 1$ |
| Left hole (attribute 3) | Dirichlet, $u = 0$ |
| Right hole (attribute 4) | Natural, $\mathbf{n} \cdot \nabla u = 0$ |

The two ends of the rectangle are identified with each other, so the solution is periodic in $x$.
Integrating the Laplacian by parts leaves the boundary term
$-\int_\Gamma (\mathbf{n} \cdot \nabla u) v$, so the Neumann and Robin conditions enter the weak form
as boundary integrals: $\int_\Gamma v$ on the right-hand side for the Neumann edge, and
$\int_\Gamma u v$ on the left-hand side plus $\int_\Gamma v$ on the right-hand side for the Robin
edge. The natural condition needs no term, and the Dirichlet condition is imposed by eliminating the
boundary values from the system.

MOOSE has no object named for the Robin condition, but [MFEMConvectiveHeatFluxBC] adds the same two
terms for the convective condition $\mathbf{n} \cdot \nabla u = h (T_\infty - u)$, which is a Robin
condition with $h = 1$ and $T_\infty = 1$ here.

`ex27p` builds its mesh in code: it joins the two ends of the rectangle, refines the mesh twice and
then moves the nodes so that the hole boundaries become circles. MOOSE cannot build this mesh
itself, so the input reads `periodic-rectangle-two-holes.mesh`, written by `ex27p` after these
steps. The periodicity is stored in the mesh, so no periodic boundary condition is needed. Each step
of `ex27p` maps onto one object in the input:

| `ex27p.cpp` | MOOSE input |
| :- | :- |
| `GenerateSerialMesh()` | [MFEMFileMesh] reading the mesh written by `GenerateSerialMesh()` |
| One parallel `UniformRefinement()` | `parallel_refine = 1` in [MFEMFileMesh] |
| `H1_FECollection` of order 1 | [MFEMScalarFESpace] with `fec_type = H1` and `fec_order = FIRST` |
| `DiffusionIntegrator` | [MFEMDiffusionKernel] |
| `BoundaryLFIntegrator` on attribute 1 | [MFEMBoundaryIntegratedBC] with `boundary = 1` |
| `MassIntegrator` and `BoundaryLFIntegrator` on attribute 2 | [MFEMConvectiveHeatFluxBC] with `boundary = 2` |
| `GetEssentialTrueDofs()` and `ProjectBdrCoefficient()` on attribute 3 | [MFEMScalarDirichletBC] with `boundary = 3` |
| `HyprePCG` preconditioned by `HypreBoomerAMG` | [MFEMHyprePCG] preconditioned by [MFEMHypreBoomerAMG] |

The test covers the default continuous discretization of `ex27p`. Its discontinuous Galerkin option
(`-dg`) needs the Neumann and Robin conditions on boundary faces, which these MOOSE objects do not
support. The test compares the $L^2$ norm of the solution, from [MFEML2Error] with `function = 0`,
against the value computed by `ex27p`.

!listing canary/test/tests/mfem/ex27p/ex27p.i

## Boundary Condition Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /BCs/MFEMConvectiveHeatFluxBC

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMBoundaryIntegratedBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMBoundaryIntegratedBC.html

[MFEMConvectiveHeatFluxBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMConvectiveHeatFluxBC.html

[MFEMScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMScalarDirichletBC.html

[MFEMHyprePCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHyprePCG.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
