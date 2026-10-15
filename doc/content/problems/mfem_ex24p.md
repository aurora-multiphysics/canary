# MFEM Example 24

A MOOSE reproduction of MFEM's parallel
[Example 24](https://github.com/mfem/mfem/blob/master/examples/ex24p.cpp) (`ex24p`), in its default
gradient form. It takes $p = \sin x \sin y \sin z$, interpolated in a first-order $H^1$ space on
MFEM's `beam-hex.mesh`, and represents its gradient in the first-order Nédélec space
$H(\mathrm{curl})$ in three ways:

1. by solving the mixed problem

   !equation
   (\mathbf{E}, \mathbf{v}) = (\nabla p, \mathbf{v}) \quad \forall \mathbf{v} \in H(\mathrm{curl})

2. by applying the discrete gradient operator to $p$, and
3. by projecting the exact gradient $\nabla p$ onto the Nédélec space.

The gradient of an $H^1$ function already lies in $H(\mathrm{curl})$, so the first two agree. Each
step of `ex24p` maps onto one object in the input:

| `ex24p.cpp` | MOOSE input |
| :- | :- |
| Two serial and one parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 1` |
| `H1_FECollection` and `ND_FECollection` of order 1 | [MFEMScalarFESpace] and [MFEMVectorFESpace] |
| `ProjectCoefficient()` of $p$ in $H^1$ | Auxiliary variable `p` with an [MFEMScalarIC] |
| `VectorFEMassIntegrator` | [MFEMVectorFEMassKernel] |
| `MixedVectorGradientIntegrator` applied to $p$ | [MFEMMixedVectorGradientKernel] with `trial_variable = p` |
| `ParDiscreteLinearOperator` with `GradientInterpolator` | [MFEMGradAux] |
| `ProjectCoefficient()` of $\nabla p$ in $H(\mathrm{curl})$ | Auxiliary variable with an [MFEMVectorIC] |
| `ComputeL2Error()` against $\nabla p$ | [MFEMVectorL2Error] |

The mixed kernel has `coefficient = -1` because MOOSE moves the terms of variables it does not
solve for, here `p`, to the right-hand side. `ex24p` solves with hypre's PCG preconditioned by
diagonal scaling, which MOOSE does not wrap, so the input uses [MFEMCGSolver] preconditioned by
[MFEMOperatorJacobiSmoother], the same Jacobi preconditioning. Both converge to a relative tolerance
of $10^{-12}$. The test compares the three $L^2$ errors against the values computed by `ex24p`.

!listing canary/test/tests/mfem/ex24p/ex24p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMMixedVectorGradientKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[MFEMScalarIC]: https://mooseframework.inl.gov/source/mfem/ics/MFEMScalarIC.html

[MFEMVectorIC]: https://mooseframework.inl.gov/source/mfem/ics/MFEMVectorIC.html

[MFEMVectorFEMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEMassKernel.html

[MFEMMixedVectorGradientKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMMixedVectorGradientKernel.html

[MFEMGradAux]: https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMGradAux.html

[MFEMVectorL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorL2Error.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMOperatorJacobiSmoother]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMOperatorJacobiSmoother.html
