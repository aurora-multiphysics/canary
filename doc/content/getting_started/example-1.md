# Example 1: Simple Diffusion

The aim of this example is to demystify MOOSE-MFEM scripts and explain the main
differences between MOOSE with MFEM and standard MOOSE. Some basic familiarity
with MOOSE input files with be assumed. Specifically, we will talk through the
differences between the standard static
[MOOSE diffusion example](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/examples/ex01_inputfile.html)
and the
[MOOSE-MFEM diffusion example](https://mooseframework.inl.gov/syntax/MFEM/SteadyStateDiffusion.html)
input files. MFEM solver based input files use some extra blocks and kernels
which may look unfamiliar compared with a standard input file.

## Input File

### [Mesh](https://mooseframework.inl.gov/syntax/Mesh/)

Firstly, for any FEM simulation, we need to mesh over our domain. This can be
done either with an auto-generated mesh or an externally generated mesh file.
Here our example opts to us a pre-generated mesh which is loaded using
`MFEMFileMesh`](https://mooseframework.inl.gov/source/mfem/mesh/MFEMFileMesh.html).

!listing test/tests/mfem/kernels/diffusion.i block=Mesh

### [Problem](https://mooseframework.inl.gov/syntax/Problem/)

The `Problem` block is missing from the MOOSE example as without a definition
it defaults to the standard
[`FEProblem`](https://mooseframework.inl.gov/source/problems/FEProblem.html).
As we wish to use the MFEM finite element library instead of MOOSE's default
libMesh we set our problem type to
[`MFEMProblem`](https://mooseframework.inl.gov/source/mfem/problem/MFEMProblem.html).

!listing test/tests/mfem/kernels/diffusion.i block=Problem

### [FESpaces](https://mooseframework.inl.gov/moose/syntax/FESpaces/index.html), [Variables](https://mooseframework.inl.gov/moose/syntax/Variables/index.html), and [AuxVariables](https://mooseframework.inl.gov/syntax/AuxVariables/)

Here we arrive at our first MFEM exclusive block.
The [`FESpaces`](https://mooseframework.inl.gov/moose/syntax/FESpaces/index.html)
block is used to specify the MFEM finite element spaces in which our problem's
variables are defined. For standard input file these are defined within the
[`Variables`](https://mooseframework.inl.gov/syntax/Variables/)
block, but due to the way MFEM's solvers work they must be defined separately.

There are four key finite element collection (FEC) types available to use:
`H1` for continuous finite elements,
`L2` for discontinuous finite elements,
`RT` for Raviart-Thomas H(div) conforming finite elements,
and `ND` for Nedelec H(curl) conforming finite elements.

The choice of FEC type determines how neighboring element's polynomial
representations must agree at their boundaries.

- For `H1` field values must agree.
- For `L2` nothing has to agree.
- For `ND` the tangential vector component agrees.
- For `RT` the normal vector component agrees.

Our example deals with a continuous variable `concentration` and a vector-field
[auxiliary variable](https://mooseframework.inl.gov/syntax/AuxVariables/)
`concentration_gradient` (discussed in more detail later).
Thus we create two first order finite element spaces to house these variables.
First order
[shape functions](https://mooseframework.inl.gov/source/problems/FEProblemBase.html#shape_functions)
will be sufficient for our basic diffusion example.

!listing test/tests/mfem/kernels/diffusion.i block=FESpaces Variables AuxVariables

### [Kernels](https://mooseframework.inl.gov/syntax/Kernels/index.html)

In scientific computing kernel is an overloaded term with many definitions.
MOOSE uses kernel to mean a piece of physics. MOOSE kernels represent terms in
the weak form of our partial differential equations. These can be selected from
a large number of pre-made kernels, or with some C++ knowledge
[custom kernels](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/examples/ex02_kernel.html)
can be created.

Here we are solving a diffusion problem so we will call our physics kernel
`diff`. For this there is an easy choice for the type as we can utilize the
pre-existing
[`MFEMDiffusionKernel`](https://mooseframework.inl.gov/source/mfem/kernels/MFEMDiffusionKernel.html).

!listing test/tests/mfem/kernels/diffusion.i start=[Kernels] end=[Solvers]

### [AuxKernels](https://mooseframework.inl.gov/syntax/AuxKernels/index.html)

Auxiliary kernels are used to compute value that can be defined explicitly with
a known function. This allows them to be calculated separately from the main
finite element solving step. Here we calculate the gradient of `concentration`
using the
[`MFEMGradAux`](https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMGradAux.html)
kernel. We further specify that this calculation should be executed as a
post-processing step by setting `execute_on` to `TIMESTEP_END`.

!listing test/tests/mfem/kernels/diffusion.i block=AuxKernels

### [BCs](https://mooseframework.inl.gov/syntax/BCs/index.html)

For boundary conditions to keep our example simple we will set our variable `u`
to 1 at the bottom and 0 at the top of our mesh. We can do so in our input file
by simply specifying top and bottom boundaries in the `BCs` section of the
MOOSE input file. The type to select for a constant value boundary condition
used with an MFEM solver is
[`MFEMScalarDirichletBC`](https://mooseframework.inl.gov/source/bcs/DirichletBC.html).

You might wonder what happens to the other boundaries which are unspecified
here. These are set to natural (zero flux) Neumann boundary conditions by
default. To see a full list of supported BCs head to
[this page](https://mooseframework.inl.gov/syntax/BCs/index.html).

!listing test/tests/mfem/kernels/diffusion.i block=BCs

### [Solvers](https://mooseframework.inl.gov/syntax/Solvers/index.html)

[`Solvers`](https://mooseframework.inl.gov/syntax/Solvers/index.html)
is another MFEM exclusive block where we select which algorithm should solve
our kernel's equations. In other words, the
[`Kernels`](https://mooseframework.inl.gov/syntax/Kernels/index.html)
block defines the mathematical operators, and the
[`Solvers`](https://mooseframework.inl.gov/syntax/Solvers/index.html)
block defines how those matrices are utilized.

The primary solver is specified under `main`. Here we decided to use
[`MFEMHypreGMRES`](https://mooseframework.inl.gov/moose/source/mfem/solvers/MFEMHypreGMRES.html)
which applies a generalized minimal residual method
([GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method))
to approximate the solution. Understanding the specific choice of solver is not
important to our example. Also defined within `main` is the linear convergence
tolerance `l_tol`, the max linear iterations `l_max_its`, and a
pre-conditioning kernel.

Pre-conditioners are used to lower the iterations required to converge the
solution of a given problem. Their aim is not to calculate the solution itself
but to transform the matrices so they are easier to solve. The pre-conditioning
solver options defined for us here are, `boomeramg` and `jacobi`.
Here we select the
[`MFEMHypreBoomerAMG`](https://mooseframework.inl.gov/moose/source/mfem/solvers/MFEMHypreBoomerAMG.html)
pre-conditioner which applies the
[Hypre BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html)
algebraic multi-grid method. The second jacobi smoother is defined here as
inactive and only included as an alternative pre-conditioner for users to
experiment with.

!listing test/tests/mfem/kernels/diffusion.i block=Solvers

### [Executioner](https://mooseframework.inl.gov/syntax/Executioner/index.html)

The [`Executioner`](https://mooseframework.inl.gov/syntax/Executioner/index.html)
block controls the approach to time stepping and the device architecture to run
the simulation on. In normal MOOSE input files this would also be where you
specify your solver, and/or your PETSc options absent here. When using MFEM
these have already been defined above in the
[`Solvers`](https://mooseframework.inl.gov/syntax/Solvers/index.html) block.
Instead we simply specify the static exectioner
[`MFEMSteady`](https://mooseframework.inl.gov/source/mfem/executioners/MFEMSteady.html).

!listing test/tests/mfem/kernels/diffusion.i block=Executioner

### [Postprocessors](https://mooseframework.inl.gov/syntax/Postprocessors/index.html)

Under `Postprocessors` we specify calculations to be run after the solve. The
aim in this example is to calculate the $L^2$ norm of the computed
concentration which is done here with
[`MFEML2Error`](https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html).
Setting `function = 0` tells
[`MFEML2Error`](https://mooseframework.inl.gov/source/mfem/postprocessors/MFEML2Error.html)
to compute the error against a constant zero function. For a full list of
post-processors head to
[this page](https://mooseframework.inl.gov/syntax/Postprocessors/index.html).

!listing test/tests/mfem/kernels/diffusion.i block=Postprocessors

### [VectorPostprocessors](https://mooseframework.inl.gov/syntax/VectorPostprocessors/index.html)

We can also specify vector post-processors that act on arrays of data. Here we
use [`MFEMVariableLineValueSampler`](https://mooseframework.inl.gov/docs/PRs/33902/site/source/mfem/vectorpostprocessors/MFEMVariableLineValueSampler.html)
to return samples of concentration within a region of interest with their
coordinates attached. This is a useful feature with the full list of available
vector post-processors [here](https://mooseframework.inl.gov/syntax/Postprocessors/index.html).

!listing test/tests/mfem/kernels/diffusion.i block=VectorPostprocessors

### Outputs

There are several MFEM specific outputs defined in this block with the `CSV`
type selected as the active output type. Users are free to experiment with the
other types. Notably
[`ParaviewDataCollection`](https://docs.mfem.org/4.8/classmfem_1_1ParaViewDataCollection.html)
provides output suitable for viewing via ParaView (seen below).

!listing test/tests/mfem/kernels/diffusion.i block=Outputs

## Running the Problem

This example can be run by using the following commands from your MOOSE
repository root:

```bash
# Go to test directory, simulation data will be output here
cd test
./moose_test-opt -i tests/mfem/kernels/diffusion.i
```

By default this will generate a csv results file. If `active` in the `Output`
block is changed to `ParaViewDataCollection` the output can be viewed from
Paraview.

!style halign=center
!media media/example-1-paraview.png style=width:70% alt=The concentration solution to diffusion.i displayed in paraview. The mug has a smooth gradient from the values of 1 at the top to 0 at the bottom.
