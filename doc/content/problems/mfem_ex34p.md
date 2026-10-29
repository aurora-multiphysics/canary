# MFEM Example 34

A MOOSE reproduction of MFEM's parallel
[Example 34](https://github.com/mfem/mfem/blob/master/examples/ex34p.cpp) (`ex34p`). It solves the
magnetostatic problem

!equation
\nabla \times \nabla \times \mathbf{A} + \delta \mathbf{A} = \mathbf{J}

for the vector potential $\mathbf{A}$ on MFEM's `fichera-mixed.mesh`, a cube with one octant
removed, built from tetrahedra, hexahedra and prisms. The current density $\mathbf{J}$ flows only
through a conductor made of five of the mesh's elements, so it is computed on a submesh of the
conductor, in two steps. First the electric potential $\phi$ solves $-\nabla^2 \phi = 0$ in the
conductor, with $\phi = 0$ on one port (boundary 2) and $\phi = 1$ on the other (boundary 23).
Then $\mathbf{J} = -\nabla \phi$ is found in $H(\mathrm{div})$, with $\mathbf{J} \cdot \mathbf{n} = 0$
on the conductor's walls: the symmetry planes (boundaries 9 to 16) and its interface with the rest
of the domain, a boundary (25) that the submesh creates. $\mathbf{J}$ is zero outside the
conductor. $\mathbf{A} \times \mathbf{n} = 0$ on the outer boundaries, and the natural condition
holds on the symmetry planes. As `ex34p` notes, this $\mathbf{J}$ is not exactly divergence free, so
the small mass term with $\delta = 10^{-6}$ keeps the system solvable.

Each step uses a different linear solver, and a MOOSE application has one, so the input is split
into three applications that pass their results on with MultiApp transfers:

- `electric_potential.i` solves for $\phi$ on the conductor with CG preconditioned by BoomerAMG.
- `current_density.i` runs `electric_potential.i`, receives $\phi$, solves for $\mathbf{J}$ on the
  conductor with CG preconditioned by Jacobi, and transfers it to the full mesh.
- `ex34p.i` runs `current_density.i`, receives $\mathbf{J}$, and solves for $\mathbf{A}$ with PCG
  preconditioned by AMS.

Each step of `ex34p` maps onto objects in these inputs:

| `ex34p.cpp` | MOOSE input |
| :- | :- |
| `SetAttribute()` on five elements of `fichera-mixed.mesh` | `fichera-mixed-conductor.mesh`, with those elements in block 2 |
| One serial and one parallel `UniformRefinement()` call | [MFEMFileMesh] with `serial_refine = 1` and `parallel_refine = 1` |
| `ParSubMesh::CreateFromDomain()` | [MFEMDomainSubMesh] with `block = 2` |
| `DiffusionIntegrator` for $\phi$ | [MFEMDiffusionKernel] |
| Mass and gradient integrators for $\mathbf{J}$ | [MFEMVectorFEMassKernel] and [MFEMMixedVectorGradientKernel] |
| $\mathbf{J} \cdot \mathbf{n} = 0$ on the conductor walls | [MFEMVectorNormalDirichletBC] |
| `ParSubMesh::Transfer()` of $\mathbf{J}$ to the full mesh | [MFEMSubMeshTransfer] |
| Passing $\phi$ and $\mathbf{J}$ between the steps | [FullSolveMultiApp] and [MultiAppMFEMCopyTransfer] |
| `CurlCurlIntegrator` and `VectorFEMassIntegrator` with $\delta$ | [MFEMCurlCurlKernel] and [MFEMVectorFEMassKernel] |
| `VectorFEDomainLFIntegrator` with $\mathbf{J}$ | [MFEMVectorFEDomainLFKernel] with `vector_coefficient = J` |
| $\mathbf{A} \times \mathbf{n} = 0$ except on the symmetry planes | [MFEMVectorTangentialDirichletBC] |
| `HyprePCG` preconditioned by `HypreAMS` | [MFEMHyprePCG] preconditioned by [MFEMHypreAMS] |
| `CurlInterpolator` for $\mathbf{B} = \nabla \times \mathbf{A}$ | [MFEMCurlAux] |

`ex34p` marks the conductor by changing the attribute of five elements after reading the mesh.
MOOSE selects submeshes by block, and cannot change the block of an element by its number, so the
input reads a copy of the mesh with those attributes already changed. `ex34p` preconditions the
current density solve with hypre's diagonal scaling, which MOOSE does not wrap, so the input uses
[MFEMOperatorJacobiSmoother], the same Jacobi preconditioning; both converge to a relative tolerance
of $10^{-12}$. The test compares the $L^2$ norms of $\mathbf{A}$ and $\mathbf{B}$, from
[MFEMVectorL2Error] with `function = '0 0 0'`, against the values computed by `ex34p`.

!listing canary/test/tests/mfem/ex34p/ex34p.i

!listing canary/test/tests/mfem/ex34p/current_density.i

!listing canary/test/tests/mfem/ex34p/electric_potential.i

## Transfer Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Transfers/MFEMSubMeshTransfer

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMDomainSubMesh]: https://mooseframework.inl.gov/source/mfem/submeshes/MFEMDomainSubMesh.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMVectorFEMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEMassKernel.html

[MFEMMixedVectorGradientKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMMixedVectorGradientKernel.html

[MFEMVectorNormalDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorNormalDirichletBC.html

[MFEMSubMeshTransfer]: https://mooseframework.inl.gov/source/mfem/transfers/MFEMSubMeshTransfer.html

[FullSolveMultiApp]: https://mooseframework.inl.gov/source/multiapps/FullSolveMultiApp.html

[MultiAppMFEMCopyTransfer]: https://mooseframework.inl.gov/source/mfem/transfers/MultiAppMFEMCopyTransfer.html

[MFEMCurlCurlKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMCurlCurlKernel.html

[MFEMVectorFEDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMVectorFEDomainLFKernel.html

[MFEMVectorTangentialDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMVectorTangentialDirichletBC.html

[MFEMHyprePCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHyprePCG.html

[MFEMHypreAMS]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreAMS.html

[MFEMCurlAux]: https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMCurlAux.html

[MFEMOperatorJacobiSmoother]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMOperatorJacobiSmoother.html

[MFEMVectorL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMVectorL2Error.html
