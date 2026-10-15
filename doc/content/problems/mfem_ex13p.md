# MFEM Example 13

A MOOSE reproduction of MFEM's parallel
[Example 13](https://github.com/mfem/mfem/blob/master/examples/ex13p.cpp) (`ex13p`). It computes the
five lowest nonzero eigenmodes of the Maxwell eigenproblem

!equation
\nabla \times \nabla \times \mathbf{E} = \lambda \mathbf{E}

with $\mathbf{E} \times \mathbf{n} = 0$ on all boundaries of MFEM's `beam-tet.mesh`, an
$8 \times 1 \times 1$ box of tetrahedra. These are the resonant modes of a closed cavity with
perfectly conducting walls, with $\lambda$ the square of each mode's wavenumber. Discretization with
first-order Nédélec elements gives the generalized eigenproblem $A \mathbf{x} = \lambda M \mathbf{x}$,
where $A$ is the curl-curl matrix and $M$ the mass matrix. Every gradient field is in the null space
of $A$, so the AME eigensolver projects these zero modes out to find the nonzero ones. Each step of
`ex13p` maps onto one object in the input:

| `ex13p.cpp` | MOOSE input |
| :- | :- |
| Two serial and one parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 1` |
| `ND_FECollection` of order 1 | [MFEMVectorFESpace] with `fec_type = ND` and `fec_order = FIRST` |
| `CurlCurlIntegrator`, forming $A$ | [MFEMCurlCurlKernel] |
| `VectorFEMassIntegrator`, forming $M$ | Formed by [MFEMEigenproblem] |
| All boundary attributes essential | [MFEMVectorTangentialDirichletBC] with `vector_coefficient = '0 0 0'` |
| `nev = 5` | [MFEMEigenproblem] with `num_modes = 5` |
| `HypreAME` preconditioned by `HypreAMS` | [MFEMHypreAME] preconditioned by [MFEMHypreAMS] |
| `SetSingularProblem()` on `HypreAMS` | [MFEMHypreAMS] with `singular = true` |
| `GetEigenvalues()` | [MFEMEigenvaluesPostprocessor] |

The test compares the five eigenvalues against those computed by `ex13p`.

!listing canary/test/tests/mfem/ex13p/ex13p.i

## Eigensolver Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Solvers/MFEMHypreAME

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[MFEMCurlCurlKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMCurlCurlKernel.html

[MFEMEigenproblem]: https://mooseframework.inl.gov/source/mfem/problem/MFEMEigenproblem.html

[MFEMVectorTangentialDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorTangentialDirichletBC.html

[MFEMHypreAME]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAME.html

[MFEMHypreAMS]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAMS.html

[MFEMEigenvaluesPostprocessor]: https://mooseframework.inl.gov/source/mfem/vectorpostprocessors/MFEMEigenvaluesPostprocessor.html
