# MFEM Example 3

A MOOSE reproduction of MFEM's parallel
[Example 3](https://github.com/mfem/mfem/blob/master/examples/ex3p.cpp) (`ex3p`). It solves the
definite Maxwell problem

!equation
\nabla \times \nabla \times \mathbf{E} + \mathbf{E} = \mathbf{f}

on MFEM's three-dimensional `beam-tet.mesh`, with the right-hand side
$\mathbf{f} = (1 + \kappa^2) \mathbf{E}$ chosen so that the exact solution is
$\mathbf{E} = (\sin \kappa y, \sin \kappa z, \sin \kappa x)$ with $\kappa = \pi$. The tangential
component of $\mathbf{E}$ is set to that of the exact solution on all boundaries. The mesh is
refined uniformly in serial as many times as keeps it at or under 1,000 elements, which is once for
the 48 elements of `beam-tet.mesh`, and then twice more in parallel, giving 32,016 unknowns with
first-order Nédélec elements. Each step of `ex3p` maps onto one object in the input:

| `ex3p.cpp` | MOOSE input |
| :- | :- |
| One serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 1` and `parallel_refine = 2` |
| `ND_FECollection` of order 1 | [MFEMVectorFESpace] with `fec_type = ND` |
| `E_exact()` and `f_exact()` | [ParsedVectorFunction] objects `E_exact` and `f_exact` |
| All boundary attributes essential, with `x` projected from `E_exact()` | [MFEMVectorTangentialDirichletBC] with `E_exact` |
| `CurlCurlIntegrator` | [MFEMCurlCurlKernel] |
| `VectorFEMassIntegrator` | [MFEMVectorFEMassKernel] |
| `VectorFEDomainLFIntegrator` with `f_exact()` | [MFEMVectorFEDomainLFKernel] with `f_exact` |
| `HyprePCG` preconditioned by `HypreAMS` | [MFEMHyprePCG] preconditioned by [MFEMHypreAMS] |

The parsed functions take $\kappa$ as `3.141592653589793`, which is the double nearest to $\pi$, the
value of `M_PI` used by `ex3p`. In `ex3p`, `HyprePCG` starts from the vector it is given, which
holds the exact solution's values on the boundary and zero elsewhere, whereas [MFEMHyprePCG] starts
from zero unless `use_initial_guess = true`. The input sets it so that both start from the same
guess.

The test compares the $L^2$ error $\|\mathbf{E}_h - \mathbf{E}\|_{L^2}$, from [MFEMVectorL2Error]
with `function = E_exact`, against the value $0.39163092340866$ that `ex3p` prints at the end of the
example, here printed at full precision. This is the discretization error of first-order Nédélec
elements on this mesh, so it is not zero: with second-order elements (`-o 2`), `ex3p` gives
$0.0158$. Run with `OMP_NUM_THREADS=1`, `ex3p` also takes the same PCG iterations as this input,
with identical residuals.

!listing canary/test/tests/mfem/ex3p/ex3p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMCurlCurlKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[ParsedVectorFunction]: https://mooseframework.inl.gov/source/functions/MooseParsedVectorFunction.html

[MFEMVectorTangentialDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorTangentialDirichletBC.html

[MFEMCurlCurlKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMCurlCurlKernel.html

[MFEMVectorFEMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEMassKernel.html

[MFEMVectorFEDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEDomainLFKernel.html

[MFEMHyprePCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHyprePCG.html

[MFEMHypreAMS]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAMS.html

[MFEMVectorL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorL2Error.html
