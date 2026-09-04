# Simple Diffusion

A steady-state solve of the Laplace equation on a unit square, driven from left to right by
Dirichlet boundary conditions. It uses only framework objects:
[Diffusion](https://mooseframework.inl.gov/source/kernels/Diffusion.html) on a first-order Lagrange
variable, [DirichletBC](https://mooseframework.inl.gov/source/bcs/DirichletBC.html) on the `left`
and `right` boundaries, and a
[Steady](https://mooseframework.inl.gov/source/executioners/Steady.html) executioner.

!listing test/tests/kernels/simple_diffusion/simple_diffusion.i

## Kernel Parameters id=parameters

The table below is generated from the `moose_test` executable in `MOOSE_DIR`, so it always matches
the MOOSE revision Canary is tested against.

!syntax parameters /Kernels/Diffusion
