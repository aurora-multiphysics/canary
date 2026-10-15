# MFEM Example 22

A MOOSE reproduction of MFEM's parallel
[Example 22](https://github.com/mfem/mfem/blob/master/examples/ex22p.cpp) (`ex22p`), in its default
scalar H1 form. It solves the complex-valued damped harmonic oscillator

!equation
-\nabla \cdot \left( \frac{1}{\mu} \nabla u \right) - \omega^2 \epsilon u + i \omega \sigma u = 0

on MFEM's `inline-quad.mesh`, the unit square, with $\mu = \epsilon = 1$, $\sigma = 20$ and
$\omega = 10$. The oscillation is driven through the boundary, where $u$ is set to the exact
solution $u = e^{-i \kappa y}$ with $\kappa = \sqrt{\mu \omega (\epsilon \omega - i \sigma)}$: a
wave travelling in $y$ that decays as it goes. Each step of `ex22p` maps onto one object in the
input:

| `ex22p.cpp` | MOOSE input |
| :- | :- |
| One serial and one parallel `UniformRefinement()` call | [MFEMFileMesh] with `serial_refine = 1` and `parallel_refine = 1` |
| `ParComplexGridFunction` on an order 1 `H1_FECollection` | [MFEMComplexVariable] on an [MFEMScalarFESpace] |
| `ParSesquilinearForm` with `DiffusionIntegrator` (real part only) | [MFEMComplexKernel] with [MFEMDiffusionKernel] components |
| `MassIntegrator`, $-\omega^2 \epsilon$ (real) and $\omega \sigma$ (imaginary) | [MFEMComplexKernel] of [MFEMMassKernel] |
| `ProjectBdrCoefficient()` with the exact solution | [MFEMComplexScalarDirichletBC] |
| `ComputeL2Error()` on the real and imaginary parts | [MFEMComplexL2Error] |

`ex22p` solves with FGMRES, preconditioned block by block with BoomerAMG built from a separate real
operator. MOOSE cannot form that preconditioner, so the input solves the system directly with
[MFEMSuperLU]. `ex22p` converges to a relative tolerance of $10^{-12}$, so the two solutions agree.
The test compares the $L^2$ error of the complex solution, $\sqrt{e_r^2 + e_i^2}$ for the errors
$e_r$ and $e_i$ of the real and imaginary parts that `ex22p` prints, against the value from `ex22p`.

!listing canary/test/tests/mfem/ex22p/ex22p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMComplexKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMComplexVariable]: https://mooseframework.inl.gov/source/mfem/variables/MFEMComplexVariable.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMComplexKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMComplexKernel.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMMassKernel.html

[MFEMComplexScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMComplexScalarDirichletBC.html

[MFEMComplexL2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEMComplexL2Error.html

[MFEMSuperLU]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMSuperLU.html
