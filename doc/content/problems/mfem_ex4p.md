# MFEM Example 4

A MOOSE reproduction of MFEM's parallel
[Example 4](https://github.com/mfem/mfem/blob/master/examples/ex4p.cpp) (`ex4p`). It solves the
$H(\mathrm{div})$ diffusion problem

!equation
-\nabla (\alpha \nabla \cdot \mathbf{F}) + \beta \mathbf{F} = \mathbf{f}

with $\alpha = \beta = 1$ on MFEM's two-dimensional `star.mesh`, with the right-hand side
$\mathbf{f} = (1 + 2\kappa^2) \mathbf{F}$ chosen so that the exact solution is
$\mathbf{F} = (\cos \kappa x \sin \kappa y, \cos \kappa y \sin \kappa x)$ with $\kappa = \pi$. The
normal component of $\mathbf{F}$ is set to that of the exact solution on all boundaries. The mesh is
refined uniformly in serial as many times as keeps it at or under 1,000 elements, which is twice for
the 20 elements of `star.mesh`, and then twice more in parallel, giving 10,400 unknowns with
lowest-order Raviart-Thomas elements. Each step of `ex4p` maps onto one object in the input:

| `ex4p.cpp` | MOOSE input |
| :- | :- |
| Two serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 2` |
| `RT_FECollection` of order 0 | [MFEMVectorFESpace] with `fec_type = RT` and `fec_order = CONSTANT` |
| `F_exact()` and `f_exact()` | [ParsedVectorFunction] objects `F_exact` and `f_exact` |
| All boundary attributes essential, with `x` projected from `F_exact()` | [MFEMVectorNormalDirichletBC] with `F_exact` |
| `DivDivIntegrator` | [MFEMDivDivKernel] |
| `VectorFEMassIntegrator` | [MFEMVectorFEMassKernel] |
| `VectorFEDomainLFIntegrator` with `f_exact()` | [MFEMVectorFEDomainLFKernel] with `f_exact` |
| `CGSolver` preconditioned by `HypreAMS` | [MFEMCGSolver] preconditioned by [MFEMHypreAMS] |

`ex4p` preconditions with AMS on two-dimensional meshes, as here, and with ADS on three-dimensional
ones. The parsed functions take $\kappa$ as `3.141592653589793`, which is the double nearest to
$\pi$, the value of `M_PI` used by `ex4p`. Both `CGSolver` and [MFEMCGSolver] start from the
boundary values of the exact solution by default, so no setting is needed to match the initial
guess.

The test compares the $L^2$ error $\|\mathbf{F}_h - \mathbf{F}\|_{L^2}$, from [MFEMVectorL2Error]
with `function = F_exact`, against the value $0.032335414550771$ that `ex4p` prints at the end of
the example, here printed at full precision. Run with `OMP_NUM_THREADS=1`, `ex4p` also takes the
same CG iterations as this input, with identical residuals.

!listing canary/test/tests/mfem/ex4p/ex4p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMDivDivKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[ParsedVectorFunction]: https://mooseframework.inl.gov/source/functions/MooseParsedVectorFunction.html

[MFEMVectorNormalDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorNormalDirichletBC.html

[MFEMDivDivKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDivDivKernel.html

[MFEMVectorFEMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEMassKernel.html

[MFEMVectorFEDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEDomainLFKernel.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreAMS]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAMS.html

[MFEMVectorL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorL2Error.html
