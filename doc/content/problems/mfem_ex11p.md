# MFEM Example 11

A MOOSE reproduction of MFEM's parallel
[Example 11](https://github.com/mfem/mfem/blob/master/examples/ex11p.cpp) (`ex11p`). It computes the
five lowest eigenmodes of the Laplacian

!equation
-\nabla^2 u = \lambda u

with $u = 0$ on all boundaries of MFEM's two-dimensional `star.mesh`. These are the fundamental
modes of vibration of a membrane of that shape clamped at its edge, with $\lambda$ the square of
each mode's frequency. Discretization with first-order Lagrange elements gives the generalized
eigenproblem $A \mathbf{x} = \lambda M \mathbf{x}$, where $A$ is the stiffness matrix and $M$ the
mass matrix. Each step of `ex11p` maps onto one object in the input:

| `ex11p.cpp` | MOOSE input |
| :- | :- |
| Two serial and one parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 1` |
| `H1_FECollection` of order 1 | [MFEMScalarFESpace] with `fec_type = H1` and `fec_order = FIRST` |
| `DiffusionIntegrator`, forming $A$ | [MFEMDiffusionKernel] |
| `MassIntegrator`, forming $M$ | Formed by [MFEMEigenproblem] |
| `MarkExternalBoundaries()` and `EliminateEssentialBCDiag()` | [MFEMScalarDirichletBC] with `coefficient = 0` |
| `nev = 5` | [MFEMEigenproblem] with `num_modes = 5` |
| `HypreLOBPCG` preconditioned by `HypreBoomerAMG` | [MFEMHypreLOBPCG] preconditioned by [MFEMHypreBoomerAMG] |
| `GetEigenvalues()` | [MFEMEigenvaluesPostprocessor] |

The test compares the five eigenvalues against those computed by `ex11p`.

!listing canary/test/tests/mfem/ex11p/ex11p.i

## Eigensolver Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Solvers/MFEMHypreLOBPCG

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMScalarDirichletBC.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMEigenproblem]: https://mooseframework.inl.gov/source/mfem/problem/MFEMEigenproblem.html

[MFEMHypreLOBPCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreLOBPCG.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEMEigenvaluesPostprocessor]: https://mooseframework.inl.gov/source/mfem/vectorpostprocessors/MFEMEigenvaluesPostprocessor.html
