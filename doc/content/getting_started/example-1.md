# Example 1: Simple Diffusion

To start let us run a basic steady-state diffusion example.
We will aim to solve the problem detailed [here](https://mooseframework.inl.gov/syntax/MFEM/SteadyStateDiffusion.html).

## Writing an input file

MOOSE-MFEM input files use some extra blocks
which the standard MOOSE diffusion example doesn't.
Here we will talk through all the differences between
the standard
[MOOSE diffusion example](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/examples/ex01_inputfile.html)
and our
[MFEM diffusion example](https://mooseframework.inl.gov/syntax/MFEM/SteadyStateDiffusion.html).

### [Mesh](https://mooseframework.inl.gov/syntax/Mesh/)

Firstly, for any FEM simulation,
we need to mesh over our domain.
This can be done either with an externally generated mesh or,
for simple meshes,
we can use the built in
[`MFEMGeneratedMeshGenerator`](https://mooseframework.inl.gov/source/mfem/meshgenerators/MFEMGeneratedMeshGenerator.html)
object.

!listing test/tests/mfem/kernels/diffusion.i block=Mesh

Here our example opts to us a pre-generated mesh
which it reads in using [`MFEMFileMesh`](https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html)

### [Problem](https://mooseframework.inl.gov/syntax/Problem/)

This block is missing from the MOOSE example
as it defaults to the standard
[`FEProblem`](https://mooseframework.inl.gov/source/problems/FEProblem.html).
As we wish to use the MFEM FE library instead of MOOSE's default libMesh
we chose to set our problem type to
[`MFEMProblem`](https://mooseframework.inl.gov/source/mfem/problem/MFEMProblem.html).

!listing test/tests/mfem/kernels/diffusion.i block=Problem

### [FESpaces](https://mooseframework.inl.gov/moose/syntax/FESpaces/index.html), [Variables](https://mooseframework.inl.gov/moose/syntax/Variables/index.html), and [AuxVariables](https://mooseframework.inl.gov/syntax/AuxVariables/)

Here we arrive at our first MFEM exclusive block.
This is used to create MFEM finite element spaces
in which our problem's variables live.
In standard MOOSE the finite-element family and order are defined within the
[`Variables`](https://mooseframework.inl.gov/syntax/Variables/)
block.
Due to the way MFEM's solvers work these must be defined separately.

There are four key space types to know:
`H1` for continuous finite elements,
`L2` for discontinuous finite elements,
`RT` for Raviart-Thomas H(div) conforming finite elements,
and `ND` for Nedelec H(curl) conforming finite elements.

The choice of finite element collection (FEC) type
determines how neighboring element's polynomial representations
must agree at their boundaries.

- `H1` field values must agree
- `L2` nothing has to agree
- `ND` the tangential vector component agrees
- `RT` the normal vector component agrees

Our example deals with a continuous
[variable](https://mooseframework.inl.gov/syntax/Variables/)
`concentration`
and a vector-field
[auxiliary variable](https://mooseframework.inl.gov/syntax/AuxVariables/)
`concentration_gradient`.
Thus we create two first order FE spaces to solve these variables in.
First order
[shape functions](https://mooseframework.inl.gov/source/problems/FEProblemBase.html#shape_functions)
should be sufficient for our basic diffusion example.

!listing test/tests/mfem/kernels/diffusion.i block=FESpaces Variables AuxVariables

### [Kernels](https://mooseframework.inl.gov/syntax/Kernels/index.html)

In scientific computing kernel is an overloaded term with many definitions.
MOOSE uses a kernel to mean a piece of physics.
They represent terms in the weak form
of our partial differential equations.
These can be selected from a large selection of pre-made kernels,
or with some C++ knowledge [created custom](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/examples/ex02_kernel.html).

Here we are solving a diffusion problem
so we'll call our physics kernel `diff`.
There is an easy choice for the type as we can utilize the prexisting kernel type
[`MFEMDiffusionKernel`](https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html).

!listing test/tests/mfem/kernels/diffusion.i start=[Kernels] end=[Solvers]

### [AuxKernels](https://mooseframework.inl.gov/syntax/AuxKernels/index.html)

Auxiliary kernels are used to compute value that can be defined explicitly
with a known function.
This allows them to be calculated separately from the main solving step.
Here we specify that we want to calculate the gradient of `concentration`
using the
[`MFEMGradAux`](https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMGradAux.html)
kernel.
We further specify to execute this calculation as a post-processing step by
setting `execute_on` to `TIMESTEP_END`.
As we are solving a steady state problem
this is done after the main calculation completes.

!listing test/tests/mfem/kernels/diffusion.i block=AuxKernels

### [BCs](https://mooseframework.inl.gov/syntax/BCs/index.html)

To keep our example simple we will set `u` to 1 at the bottom boundary
and 0 at the top boundary of of our mesh.
We can do so in our input file by simply specifying top and bottom
boundaries in the `BCs` section of the MOOSE input file.

The type to select for a constant value boundary condition is
[`MFEMScalarDirichletBC`](https://mooseframework.inl.gov/source/bcs/DirichletBC.html).
We must also specify the variable the boundary will apply to
along with its value.
You might wonder what happens to the other boundaries
which are unspecified here.
These are set to natural (zero flux)
Neumann boundary conditions by default.
To see a full list of supported BCs see [this page](https://mooseframework.inl.gov/syntax/BCs/index.html).

!listing test/tests/mfem/kernels/diffusion.i block=BCs

### [Solvers](https://mooseframework.inl.gov/syntax/Solvers/index.html)

This is another MFEM exclusive block where
we will select which algorithm
should solve our kernel's equations.
In other words, the [`Kernels`](https://mooseframework.inl.gov/syntax/Kernels/index.html)
block defines the mathematical operators,
and the [`Solvers`](https://mooseframework.inl.gov/syntax/Solvers/index.html)
block defines how those matrices are inverted.

Our primary solver is selected under the `main` option.
Here we decided upon
[`MFEMHypreGMRES`](https://mooseframework.inl.gov/moose/source/mfem/solvers/MFEMHypreGMRES.html)
which applies a Generalized Minimal Residual Method
([GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method))
to approximate the solution.
The understanding the specific solvers selected is not important to our example.

There are other definitions within this block,
the first of which is
A preconditioner is used to lower the required iterations to solve the problem.

This applies a the
[`MFEMHypreBoomerAMG`](https://mooseframework.inl.gov/moose/source/mfem/solvers/MFEMHypreBoomerAMG.html)
preconditioner which applies the
[Hypre BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html)
algebraic multi-grid method.

!listing test/tests/mfem/kernels/diffusion.i block=Solvers

### [Executioner](https://mooseframework.inl.gov/syntax/Executioner/index.html)

Simply put the executioner block defines how the code is executed.
Here we can control the solving method,
approach to time stepping,
and or control [PETSc](https://petsc.org/release/documentation/) behavior.
PETSc is used as the underlying solver in MOOSE
and its options can be specified here.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Executioner

The two main types of solver for most simulations are
[`NEWTON`](https://mooseframework.inl.gov/modules/phase_field/Solving.html)
or
[`PJFNK`](https://mooseframework.inl.gov/modules/phase_field/Solving.html)
(Preconditioned Jacobian Free Network Newton Krylov).

### Outputs

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Outputs
