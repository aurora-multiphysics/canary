# MFEM Example 2

A MOOSE reproduction of MFEM's parallel
[Example 2](https://github.com/mfem/mfem/blob/master/examples/ex2p.cpp) (`ex2p`). It solves the
linear elasticity problem

!equation
-\nabla \cdot \sigma(\mathbf{u}) = 0, \qquad
\sigma(\mathbf{u}) = \lambda (\nabla \cdot \mathbf{u}) I + \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T)

for a cantilever beam made of two materials, with $\lambda = \mu = 50$ in material 1 and
$\lambda = \mu = 1$ in material 2. The beam, MFEM's two-dimensional `beam-tri.mesh`, is fixed with
$\mathbf{u} = 0$ on boundary 1 and pulled down by the traction
$\sigma(\mathbf{u}) \cdot \mathbf{n} = (0, -0.01)$ on boundary 2, with no traction elsewhere. The
mesh is refined uniformly in serial as many times as keeps it at or under 1,000 elements, which is
twice for the 16 elements of `beam-tri.mesh`, and then once more in parallel, giving 1,170 unknowns
with first-order vector $H^1$ elements. Each step of `ex2p` maps onto one object in the input:

| `ex2p.cpp` | MOOSE input |
| :- | :- |
| Two serial and one parallel `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 2` and `parallel_refine = 1` |
| Vector `H1_FECollection` of order 1, ordered `byVDIM` | [MFEMVectorFESpace] with `fec_type = H1` and `ordering = VDIM` |
| Essential boundary attribute 1 with `x = 0` | [MFEMVectorDirichletBC] on boundary `1` |
| `VectorBoundaryLFIntegrator` with the pull-down force | [MFEMVectorBoundaryIntegratedBC] on boundary `2` |
| `PWConstCoefficient` values of $\lambda$ and $\mu$ | [MFEMGenericFunctorMaterial] on blocks `1` and `2` |
| `ElasticityIntegrator` | [MFEMLinearElasticityKernel] |
| `HyprePCG` preconditioned by `HypreBoomerAMG` | [MFEMHyprePCG] preconditioned by [MFEMHypreBoomerAMG] |

The input passes the vector space to [MFEMHypreBoomerAMG] through `fespace`, which configures it with
`SetElasticityOptions()`, as `ex2p` does when run with `-elast`; by default `ex2p` uses
`SetSystemsOptions()` instead. The preconditioner only changes the path to the solution: the two runs
of `ex2p` give displacement $L^2$ norms that agree to about $10^{-11}$, relative.

The test compares the $L^2$ norm of the displacement, from [MFEMVectorL2Error] with
`function = '0 0'`, against the value $\|\mathbf{u}\|_{L^2} = 1.1067699081498$ printed by the
default `ex2p`. The norm is computed with a zero `VectorConstantCoefficient` straight after
`RecoverFEMSolution()`, because `ex2p` goes on to displace the mesh and invert the solution for
output. Run with `-elast` and `OMP_NUM_THREADS=1`, `ex2p` also takes the same PCG iterations as this
input, with identical residuals, and gives the same norm to every printed digit.

!listing canary/test/tests/mfem/ex2p/ex2p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMLinearElasticityKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMVectorFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMVectorFESpace.html

[MFEMVectorDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorDirichletBC.html

[MFEMVectorBoundaryIntegratedBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorBoundaryIntegratedBC.html

[MFEMGenericFunctorMaterial]: https://mooseframework.inl.gov/source/mfem/functormaterials/MFEMGenericFunctorMaterial.html

[MFEMLinearElasticityKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMLinearElasticityKernel.html

[MFEMHyprePCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHyprePCG.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEMVectorL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorL2Error.html
