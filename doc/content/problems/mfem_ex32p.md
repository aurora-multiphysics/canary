# MFEM Example 32

A MOOSE reproduction of MFEM's parallel
[Example 32](https://github.com/mfem/mfem/blob/master/examples/ex32p.cpp) (`ex32p`). It computes the
five lowest nonzero eigenmodes of the Maxwell eigenproblem

!equation
\nabla \times \nabla \times \mathbf{E} = \lambda \varepsilon \mathbf{E}

with $\mathbf{E} \times \mathbf{n} = 0$ on all boundaries of MFEM's `inline-quad.mesh`, a unit
square. The mesh is two-dimensional, but the field $\mathbf{E}$ has all three components, each
depending on $x$ and $y$ only. The dielectric tensor is constant but anisotropic,

!equation
\varepsilon = \begin{pmatrix} 2 & 1/\sqrt{2} & 0 \\ 1/\sqrt{2} & 2 & 1/\sqrt{2} \\ 0 & 1/\sqrt{2} & 2 \end{pmatrix},

so it couples the in-plane components of $\mathbf{E}$ to the out-of-plane one. The field is
discretized with first-order Nédélec elements for three-dimensional vectors in two dimensions,
which keep the tangential in-plane components continuous between elements, as well as the
out-of-plane one. This gives the generalized eigenproblem $A \mathbf{x} = \lambda M \mathbf{x}$, where
$A$ is the curl-curl matrix and $M$ the mass matrix weighted by $\varepsilon$. As in MFEM's
[Example 13](https://github.com/mfem/mfem/blob/master/examples/ex13p.cpp), every gradient field is in
the null space of $A$, so the AME eigensolver projects these zero modes out to find the nonzero ones.
Each step of `ex32p` maps onto one object in the input:

| `ex32p.cpp` | MOOSE input |
| :- | :- |
| Two serial and one parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 1` |
| `ND_R2D_FECollection` of order 1 | [MFEMVectorFESpace] with `fec_type = ND` and `range_dim = 3` |
| `CurlCurlIntegrator`, forming $A$ | [MFEMCurlCurlKernel] |
| `VectorFEMassIntegrator` with the matrix $\varepsilon$, forming $M$ | [MFEMEigenproblem] with `rhs_matrix_coefficient` |
| All boundary attributes essential | [MFEMVectorTangentialDirichletBC] with `vector_coefficient = '0 0'` |
| `nev = 5` | [MFEMEigenproblem] with `num_modes = 5` |
| `HypreAME` preconditioned by `HypreAMS` | [MFEMHypreAME] preconditioned by [MFEMHypreAMS] |
| `SetSingularProblem()` on `HypreAMS` | [MFEMHypreAMS] with `singular = true` |
| `GetEigenvalues()` | [MFEMEigenvaluesPostprocessor] |

The zero tangential condition is given with two components: MFEM checks the boundary coefficient of
this space against two-component vectors, and rejects a three-component one. Zero values project to
zero either way, so the condition and the eliminated degrees of freedom are the same as in `ex32p`.
The test compares the five eigenvalues against those computed by `ex32p`.

!listing canary/test/tests/mfem/ex32p/ex32p.i

## Eigenproblem Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Problem/MFEMEigenproblem

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[MFEMCurlCurlKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMCurlCurlKernel.html

[MFEMEigenproblem]: https://mooseframework.inl.gov/source/mfem/problem/MFEMEigenproblem.html

[MFEMVectorTangentialDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorTangentialDirichletBC.html

[MFEMHypreAME]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAME.html

[MFEMHypreAMS]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAMS.html

[MFEMEigenvaluesPostprocessor]: https://mooseframework.inl.gov/source/mfem/vectorpostprocessors/MFEMEigenvaluesPostprocessor.html
