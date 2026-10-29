# MFEM Example 26

A MOOSE reproduction of MFEM's parallel
[Example 26](https://github.com/mfem/mfem/blob/master/examples/ex26p.cpp) (`ex26p`). It solves the
Poisson problem of MFEM's Examples 0 and 1,

!equation
-\nabla^2 u = 1

with $u = 0$ on all boundaries of MFEM's `star.mesh`, now with fourth-order Lagrange elements and
a p-multigrid preconditioner. The multigrid works on a hierarchy of first, second and fourth-order
spaces on the same mesh. Each level above the coarsest is smoothed by a second-order Chebyshev
smoother, which needs only the action of the operator and its diagonal, so the fourth-order
problem can be partially assembled instead of stored as a matrix. The first-order level is solved
approximately by a few iterations of CG preconditioned by BoomerAMG. Each step of `ex26p` maps onto
one object in the input:

| `ex26p.cpp` | MOOSE input |
| :- | :- |
| Two serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 2` |
| `ParFiniteElementSpaceHierarchy` with orders 1, 2 and 4 | [MFEMFESpaceHierarchy] with `refinements = '2 4'` |
| `DiffusionIntegrator` and `DomainLFIntegrator` | [MFEMDiffusionKernel] and [MFEMDomainLFKernel] |
| `GeometricMultigrid` with a V-cycle | [MFEMGeometricMultigridSolver] |
| `OperatorChebyshevSmoother` of order 2 | [MFEMOperatorChebyshevSmoother] |
| `CGSolver` with `HypreBoomerAMG` on the coarsest level | [MFEMCGSolver] with [MFEMHypreBoomerAMG] |
| `AssemblyLevel::PARTIAL` on the finest level | `assembly_level = partial` in the `[Executioner]` |

`ex26p` partially assembles the second-order level as well, but
[MFEMGeometricMultigridSolver] currently crashes when a level other than the finest is partially
assembled, so the input assembles that level fully. This gives the same operator. The test
compares the $L^2$ norm of the solution, from [MFEML2Error] with `function = 0`, against the value
computed by `ex26p`.

!listing canary/test/tests/mfem/ex26p/ex26p.i

## Solver Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Solvers/MFEMGeometricMultigridSolver

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMFESpaceHierarchy]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMFESpaceHierarchy.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDomainLFKernel.html

[MFEMGeometricMultigridSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMGeometricMultigridSolver.html

[MFEMOperatorChebyshevSmoother]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMOperatorChebyshevSmoother.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
