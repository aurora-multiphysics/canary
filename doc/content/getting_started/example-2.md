# Example 2: Transient Heat Transfer

Now lets look at time stepping solvers with the MOOSE-MFEM
[heat transfer problem](https://mooseframework.inl.gov/syntax/MFEM/HeatTransfer.html),
also solved using the mug mesh. Here we will focus on the differences moving
from a static to a time stepping solver with MFEM options. This problem remains
similar to the previous diffusion example but introduces an optional convective
boundary at the top side of the mug mesh modeling heat loss to a surrounding
medium.

## Input File

### Time Derivative Kernel

An immediate difference to the static diffusion input file is that within the
`Kernel` block we have defined an extra kernel. The diffusion kernel seen in
example 1 remains but there is an extra temperature time derivative kernel
`dT_dt` selected as
[`MFEMTimeDerivativeMassKernel`](https://mooseframework.inl.gov/source/mfem/kernels/MFEMTimeDerivativeMassKernel.html).
This adds the weak form of the equation element $\rho c_p\dot{u}$ for our
solver to work with, where the coefficient $\rho c_p$ is set to 1 by default.

!listing test/tests/mfem/kernels/heattransfer.i start=[Kernels] end=[BCs]

### Convective Boundary Condition

As discussed we have defined an extra boundary condition to represent heat
exchange between solid surface and a surrounding fluid. This is governed by the
following equation

!equation id=convection equation
q_\text{n} = h(T - T_\infty)

where $q_\text{n}$ is the outward heat flux from the surface, $h$ is the
convective heat transfer coefficient, and $T_\infty$ is the surrounding
temperature. Users have the option to switch the active top boundary to a fixed
dirichlet boundary by changing the BCs listed as active to
`active = 'bottom top_dirichlet'`.

Unlike Example 1, this input identifies boundaries by their numeric IDs rather
than their names. For the mug mesh used in both examples, `boundary = '1'`
selects the bottom surface and `boundary = '2'` selects the top, equivalent to
`boundary = 'bottom'` and `boundary = 'top'`, respectively. These names and IDs
are defined by the mesh, so the correspondence may differ for other meshes.

!listing test/tests/mfem/kernels/heattransfer.i block=BCs

### Time Stepping

The transience in our equation requires the executioner to solve across time.
This is trivially done by selecting
[`MFEMTransient`](https://mooseframework.inl.gov/source/mfem/executioners/MFEMTransient.html)
as the type in the
[`Executioner`](https://mooseframework.inl.gov/syntax/Executioner/index.html)
block and setting a time step `dt`, a `start_time`, and an `end_time`.
This applies the default implicit Euler scheme to our problem.

!listing test/tests/mfem/kernels/heattransfer.i block=Executioner

### Assembly Level

In the
[`Executioner`](https://mooseframework.inl.gov/syntax/Executioner/index.html)
block we expose the `assembly_level` setting which defaults to `legacy`. This
tells MFEM how to assemble the
[FEM operator matrices](https://mfem.org/howto/assembly_levels/).
This can be helpful for advanced users and allows for optmisation decisions
balancing memory usage vs compute along side archetectural decisions (e.g. if
using GPU acceleration is is advisable to select a level other than `legacy`).
For now we will not utilize this setting but it is useful to have an awareness
off.

!listing test/tests/mfem/kernels/heattransfer.i block=Executioner start=assembly_level end=dt

### Inactive Auxiliary Variable

Also present in this input file is the inactive auxiliary variable aimed to
calculate the running average temperature. This is calculated using the
[`MFEMScalarTimeAverageAux` ](https://mooseframework.inl.gov/source/mfem/auxkernels/MFEMScalarTimeAverageAux.html)
kernel but will not run unless the user removes the inactive label from both
the `AuxVariables` and `AuxKernels` block.

!listing test/tests/mfem/kernels/heattransfer.i block=AuxVariables AuxKernels

## Running the Problem

This example can be run by using the following commands from your MOOSE
repository root:

```bash
# Go to test directory
cd test
./moose_test-opt -i tests/mfem/kernels/heattransfer.i
# output expected in:
# <moose-root>/test/OutputData/HeatTransfer
```

The current configuration produces the following CSV line samples:

```
heattransfer_line_sample_0001.csv
heattransfer_line_sample_0002.csv
heattransfer_line_sample_0003.csv
```

!content pagination previous=getting_started/example-1.md next=getting_started/example-3.md
