# MFEM Example 39

A MOOSE reproduction of MFEM's parallel
[Example 39](https://github.com/mfem/mfem/blob/master/examples/ex39p.cpp) (`ex39p`). It solves the
Poisson problem

!equation
-\nabla \cdot (k \nabla u) = f

on a compass rose, with $u = 0$ on the whole boundary. The mesh, `compass.msh`, is a Gmsh file whose
regions have names: the eight boundary segments are named after compass headings (`NNE`, `ENE`, and
so on), and each of the four points of the rose is split into an `Even` and an `Odd` half, with the
remaining elements named `Base`. The diffusion coefficient is $10^{-6}$ everywhere, plus 1 on `Base`
and 2 on the four `Even` halves. The source is $f = 1$ on the `Even` halves and zero elsewhere.

MFEM reads the names into attribute sets, and MOOSE looks up a name given as a `block` or `boundary`
in these sets, so regions can be selected by name in the input. `ex39p` also creates new sets in
code, `Rose Even` for the four `Even` halves and `Boundary` for the eight boundary segments. MOOSE
cannot define new sets, so the input lists the sets these are made of. It also cannot take a name
that contains a space, such as `N Even`, because it splits lists on spaces, so the `Even` halves are
given by their attribute numbers.

Each step of `ex39p` maps onto objects in the input:

| `ex39p.cpp` | MOOSE input |
| :- | :- |
| Four serial and two parallel `UniformRefinement()` calls | [MFEMFileMesh], `serial_refine = 4`, `parallel_refine = 2` |
| `H1_FECollection` of order 1 | [MFEMScalarFESpace] with `fec_type = H1` and `fec_order = FIRST` |
| `DiffusionIntegrator` everywhere, on `Base` and on `Rose Even` | Three [MFEMDiffusionKernel], the last two with `block` |
| `DomainLFIntegrator` on `Rose Even` | [MFEMDomainLFKernel] with `block` |
| Essential boundary on the set `Boundary` | [MFEMScalarDirichletBC] on the eight boundary names |
| `CGSolver` preconditioned by `HypreBoomerAMG` | [MFEMCGSolver] preconditioned by [MFEMHypreBoomerAMG] |

The diffusion integrators are added in the same order as in `ex39p`, so the assembled matrix and the
CG history are identical. The test compares the $L^2$ norm of the solution, from [MFEML2Error] with
`function = 0`, against the value computed by `ex39p`.

!listing canary/test/tests/mfem/ex39p/ex39p.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/MFEMDiffusionKernel

[MFEMFileMesh]: https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html

[MFEMScalarFESpace]: https://mooseframework.inl.gov/source/mfem/fespaces/MFEMScalarFESpace.html

[MFEMDiffusionKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html

[MFEMDomainLFKernel]: https://mooseframework.inl.gov/source/mfem/kernels/MFEMDomainLFKernel.html

[MFEMScalarDirichletBC]: https://mooseframework.inl.gov/source/mfem/bcs/MFEMScalarDirichletBC.html

[MFEMCGSolver]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMCGSolver.html

[MFEMHypreBoomerAMG]: https://mooseframework.inl.gov/source/mfem/solvers/MFEMHypreBoomerAMG.html

[MFEML2Error]: https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html
