# MFEM Example 1

A MOOSE reproduction of MFEM's parallel
[Example 1](https://github.com/mfem/mfem/blob/master/examples/ex1p.cpp) (`ex1p`). It solves the
Poisson problem

!equation
-\nabla^2 u = 1

with $u = 0$ on the external boundary of MFEM's two-dimensional `star.mesh`, using first-order $H^1$
elements. The mesh is refined uniformly in serial as many times as keeps it at or under 10,000
elements, which is four times for the 20 elements of `star.mesh`, and then twice more in parallel,
giving 82,561 unknowns. Each step of `ex1p` maps onto one object in the input:

| `ex1p.cpp` | MOOSE input |
| :- | :- |
| Four serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 4` and `parallel_refine = 2` |
| `H1_FECollection` of order 1 | [MFEMScalarFESpace] with `fec_type = H1` |
| `MarkExternalBoundaries()` with `x = 0` | [MFEMScalarDirichletBC], on all boundaries by default |
| `DiffusionIntegrator` | [MFEMDiffusionKernel] |
| `DomainLFIntegrator` with coefficient 1 | [MFEMDomainLFKernel] |
| `CGSolver` preconditioned by `HypreBoomerAMG` | [MFEMCGSolver] preconditioned by [MFEMHypreBoomerAMG] |

`ex1p` computes its number of serial refinements from the mesh and constrains only external
boundaries, whereas `serial_refine` is a fixed count and [MFEMScalarDirichletBC] applies to every
boundary attribute. The two agree on `star.mesh`, whose boundary is entirely external, but the input
would need adjusting for other meshes.

`ex1p` applies BoomerAMG as a single V-cycle per CG iteration, whereas [MFEMHypreBoomerAMG] iterates
to a tolerance of `1e-5` by default, so the input sets `l_max_its = 1` and `l_tol = 0` to match.

The test compares the $L^2$ norm of the solution, from [MFEML2Error] with `function = 0`, against
the value $\|u\|_{L^2} = 0.37972504807552$ printed by `ex1p` after adding `x.ComputeL2Error(zero)`
with a zero `ConstantCoefficient` at the end of the example. Run with `OMP_NUM_THREADS=1`, `ex1p`
also takes the same CG iterations as this input, with identical residuals.

!listing canary/test/tests/mfem/ex1p/ex1p.i

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
