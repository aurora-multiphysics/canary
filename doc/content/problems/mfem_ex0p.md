# MFEM Example 0

A MOOSE reproduction of MFEM's parallel
[Example 0](https://github.com/mfem/mfem/blob/master/examples/ex0p.cpp) (`ex0p`), the most basic
MFEM example. It solves the Poisson problem

!equation
-\nabla^2 u = 1

with $u = 0$ on the boundary of MFEM's two-dimensional `star.mesh`, refined once in parallel, using
first-order $H^1$ elements. Each step of `ex0p` maps onto one object in the input:

| `ex0p.cpp` | MOOSE input |
| :- | :- |
| `ParMesh` with one `UniformRefinement()` | [MFEMFileMesh] with `parallel_refine = 1` |
| `H1_FECollection` of order 1 | [MFEMScalarFESpace] with `fec_type = H1` |
| `GetBoundaryTrueDofs()` with `x = 0` | [MFEMScalarDirichletBC], on all boundaries by default |
| `DiffusionIntegrator` | [MFEMDiffusionKernel] |
| `DomainLFIntegrator` with coefficient 1 | [MFEMDomainLFKernel] |
| `CGSolver` preconditioned by `HypreBoomerAMG` | [MFEMCGSolver] preconditioned by [MFEMHypreBoomerAMG] |

`ex0p` applies BoomerAMG as a single V-cycle per CG iteration, whereas [MFEMHypreBoomerAMG] iterates
to a tolerance of `1e-5` by default, so the input sets `l_max_its = 1` and `l_tol = 0` to match.

The test compares the $L^2$ norm of the solution, from [MFEML2Error] with `function = 0`, against
the value $\|u\|_{L^2} = 0.37438237056131$ printed by `ex0p` after adding `x.ComputeL2Error(zero)`
with a zero `ConstantCoefficient` at the end of the example. Run with `OMP_NUM_THREADS=1`, `ex0p`
also takes the same CG iterations as this input, with identical residuals.

!listing canary/test/tests/mfem/ex0p/ex0p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMDiffusionKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMScalarDirichletBC.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDomainLFKernel.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
