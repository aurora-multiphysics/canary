# MFEM Example 33

A MOOSE reproduction of MFEM's parallel
[Example 33](https://github.com/mfem/mfem/blob/master/examples/ex33p.cpp) (`ex33p`). It solves the
fractional Poisson problem

!equation
(-\nabla^2)^\alpha u = 1

with $\alpha = 1/2$ and $u = 0$ on all boundaries of MFEM's `star.mesh`. The fractional power is
defined through the eigenfunctions $\phi_k$ of the Laplacian with these boundary conditions: if
$-\nabla^2 \phi_k = \lambda_k \phi_k$, then $(-\nabla^2)^\alpha$ scales each $\phi_k$ by
$\lambda_k^\alpha$. Unlike the Laplacian, this operator is nonlocal. `ex33p` avoids computing
eigenfunctions by approximating $\lambda^{-\alpha}$ with a rational function in partial fractions,

!equation
\lambda^{-\alpha} \approx \sum_{i=0}^{11} \frac{c_i}{\lambda + d_i},

which holds for every eigenvalue at once. Applied to the right-hand side, each term is a
reaction-diffusion problem,

!equation
-\nabla^2 x_i + d_i x_i = 1, \qquad x_i = 0 \text{ on the boundary},

and the solution is the weighted sum $u = \sum_i c_i x_i$. The 12 problems are independent, so
`ex33p` solves them one after another. In the input, each problem is its own weak form, solved by
its own problem composer, and an auxiliary kernel forms the weighted sum. Each step of `ex33p` maps
onto one object in the input:

| `ex33p.cpp` | MOOSE input |
| :- | :- |
| Three `UniformRefinement()` calls | [MFEMFileMesh] with `serial_refine = 3` |
| `ComputePartialFractionApproximation()`, giving $c_i$ and $d_i$ | Values written at the top of the input |
| `DiffusionIntegrator` and `MassIntegrator` with coefficient $d_i$ | [MFEMDiffusionKernel] and [MFEMMassKernel] |
| `DomainLFIntegrator` | [MFEMDomainLFKernel] |
| One bilinear form and solve for each term | One [MFEMWeakForm] and one [MFEMWeakFormProblemComposer] for each term |
| `CGSolver` preconditioned by `HypreBoomerAMG`, for each term | [MFEMCGSolver] preconditioned by [MFEMHypreBoomerAMG], shared |
| `u += c_i x_i` | [MFEMSumAux] with `scale_factors = c_i` |

The coefficients $c_i$ and $d_i$ are the 12 that `ex33p` uses when MFEM is built without LAPACK, as
it is in MOOSE. With LAPACK, `ex33p` computes them when it runs, with the AAA algorithm. The single
solver is set up again for each term's matrix before that term is solved, as `ex33p` sets up a new
one. The test compares the $L^2$ norm of $u$, from [MFEML2Error] with `function = 0`, against the
value computed by `ex33p`.

!listing canary/test/tests/mfem/ex33p/ex33p.i

## Problem Composer Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /ProblemComposers/MFEMWeakFormProblemComposer

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMMassKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMMassKernel.html

[MFEMDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDomainLFKernel.html

[MFEMWeakForm]: https://mooseframework.inl.gov/source/mfem/weakforms/MFEMWeakForm.html

[MFEMWeakFormProblemComposer]: https://mooseframework.inl.gov/source/mfem/problem_composers/MFEMWeakFormProblemComposer.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEMSumAux]: https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMSumAux.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
