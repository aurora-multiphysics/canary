# MFEM Example 35

A MOOSE reproduction of MFEM's parallel
[Example 35](https://github.com/mfem/mfem/blob/master/examples/ex35p.cpp) (`ex35p`), in its default
scalar form. It solves the damped harmonic oscillator

!equation
-\nabla^2 u - \omega^2 u + i \omega \sigma u = 0

for a complex field $u$ on MFEM's `fichera-mixed.mesh`, with angular frequency $\omega = 2\pi$ and
damping $\sigma = 2$. The oscillation is driven through the boundary: on a port made of boundaries 7,
8, 11 and 12, $u$ is the second lowest eigenmode of the Laplacian on the port, with zero values on
the port's edges; on the rest of the boundary $u = 0$. This is how a waveguide is fed with one of
its modes. The boundary value is real, but the damping term couples the real and imaginary parts of
$u$, so the solution is complex inside the domain.

The port mode is an eigenproblem and the oscillator a complex problem, two different problem types
in MOOSE, so the input is split into two applications:

- `port_mode.i` computes five eigenmodes on a submesh made of the port boundaries, with LOBPCG
  preconditioned by BoomerAMG and the random seed of `ex35p`, and transfers the second one to the
  full mesh.
- `ex35p.i` runs `port_mode.i`, receives the mode, and solves the complex problem with it as the
  boundary value.

Each step of `ex35p` maps onto objects in these inputs:

| `ex35p.cpp` | MOOSE input |
| :- | :- |
| One serial and one parallel `UniformRefinement()` call | [MFEMFileMesh] with `serial_refine = 1` and `parallel_refine = 1` |
| `ParSubMesh::CreateFromBoundary()` | [MFEMBoundarySubMesh] |
| `H1_FECollection` of order 1 on the two-dimensional port | [MFEMGenericFESpace] with `fec_name = H1_2D_P1` |
| `DiffusionIntegrator` and `MassIntegrator` for the port modes | [MFEMDiffusionKernel]; [MFEMEigenproblem] forms the mass |
| `HypreLOBPCG` preconditioned by `HypreBoomerAMG` | [MFEMHypreLOBPCG] preconditioned by [MFEMHypreBoomerAMG] |
| `ParSubMesh::Transfer()` of the mode to the full mesh | [MFEMSubMeshTransfer] from `port_mode_1` |
| Passing the mode to the oscillator | [FullSolveMultiApp] and [MultiAppMFEMShapeEvaluationTransfer] |
| `ParSesquilinearForm` with diffusion and mass terms | [MFEMComplexKernel] of [MFEMDiffusionKernel] and [MFEMMassKernel] |
| All boundary attributes essential, with the mode on the port | [MFEMComplexScalarDirichletBC] |
| SuperLU, with `-slu` | [MFEMSuperLU] |

A complex problem in MOOSE holds only complex variables, so the oscillator receives the real port
mode as the real part of a complex field, with a transfer that evaluates the mode at the nodes of the
full mesh. By default `ex35p` solves the complex system with FGMRES to a relative tolerance of
$10^{-6}$, preconditioned block by block with BoomerAMG on a separate real operator. MOOSE cannot
form that preconditioner, so the input solves directly with [MFEMSuperLU], as `ex35p -slu` does. The
test compares the $L^2$ norm of $u$, $\sqrt{\|u_r\|^2 + \|u_i\|^2}$ from [MFEMComplexL2Error] with
zero functions, against the value computed by `ex35p -slu`, which differs from the default
iterative solution by about $2 \times 10^{-8}$ relative.

!listing canary/test/tests/mfem/ex35p/ex35p.i

!listing canary/test/tests/mfem/ex35p/port_mode.i

## Transfer Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Transfers/MultiAppMFEMShapeEvaluationTransfer

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMBoundarySubMesh]: https://mooseframework.inl.gov/source/mfem/submeshes/MFEMBoundarySubMesh.html

[MFEMGenericFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMGenericFESpace.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMEigenproblem]: https://mooseframework.inl.gov/source/mfem/problem/MFEMEigenproblem.html

[MFEMHypreLOBPCG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreLOBPCG.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEMSubMeshTransfer]: https://mooseframework.inl.gov/source/mfem/transfers/MFEMSubMeshTransfer.html

[FullSolveMultiApp]: https://mooseframework.inl.gov/source/multiapps/FullSolveMultiApp.html

[MultiAppMFEMShapeEvaluationTransfer]: https://mooseframework.inl.gov/source/mfem/transfers/MultiAppMFEMShapeEvaluationTransfer.html

[MFEMComplexKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMComplexKernel.html

[MFEMMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMMassKernel.html

[MFEMComplexScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMComplexScalarDirichletBC.html

[MFEMSuperLU]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMSuperLU.html

[MFEMComplexL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMComplexL2Error.html
