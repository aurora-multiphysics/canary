# MFEM Example 14

A MOOSE reproduction of MFEM's parallel
[Example 14](https://github.com/mfem/mfem/blob/master/examples/ex14p.cpp) (`ex14p`). It solves the
Poisson problem

!equation
-\nabla^2 u = 1

with $u = 0$ on all boundaries of MFEM's two-dimensional `star.mesh`, using a discontinuous Galerkin
(DG) discretization. The solution may jump between elements, and continuity is enforced weakly by
terms on the element faces that penalize the jumps. The Dirichlet condition is imposed the same way,
on the boundary faces, instead of being eliminated from the system. The face terms are those of the
symmetric interior penalty method, with symmetry parameter $\sigma = -1$ and penalty parameter
$\kappa = (p + 1)^2 = 4$ for elements of order $p = 1$. Each step of `ex14p` maps onto one object
in the input:

| `ex14p.cpp` | MOOSE input |
| :- | :- |
| Four serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 4` and `parallel_refine = 2` |
| `DG_FECollection` of order 1, Gauss-Legendre basis | [MFEMScalarFESpace] with `fec_type = L2` and `basis = GaussLegendre` |
| `DiffusionIntegrator` | [MFEMDiffusionKernel] |
| `DGDiffusionIntegrator` on interior faces | [MFEMDGDiffusionKernel] |
| `DGDiffusionIntegrator` and `DGDirichletLFIntegrator` on boundary faces | [MFEMDGDiffusionBC] |
| `DomainLFIntegrator` | [MFEMDomainLFKernel] |
| `CGSolver` preconditioned by `HypreBoomerAMG` | [MFEMCGSolver] preconditioned by [MFEMHypreBoomerAMG] |

[MFEMDGDiffusionKernel] and [MFEMDGDiffusionBC] default to the same $\sigma$ and $\kappa$ as
`ex14p`. The test compares the $L^2$ norm of the solution, from [MFEML2Error] with `function = 0`,
against the value computed by `ex14p`.

!listing canary/test/tests/mfem/ex14p/ex14p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMDGDiffusionKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMDGDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDGDiffusionKernel.html

[MFEMDGDiffusionBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMDGDiffusionBC.html

[MFEMDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDomainLFKernel.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
