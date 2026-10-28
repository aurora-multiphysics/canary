# Example 2: Transient Heat Transfer

Now lets look at time stepping solvers with a
[heat transfer problem](https://mooseframework.inl.gov/syntax/MFEM/HeatTransfer.html),
also solved on the mug mesh.
Here we will focus on the differences moving from
a static to a time stepping solver.
This problem introduces an optional convective boundary at the top of the mug
to simulate

## Input File

### Time Derivative Kernel

An immediate difference is that we have defined two kernels
within the `Kernel` block.
One remains the diffusion kernel seen in example 1
but there is an extra temperature time derivative kernel `dT_dt` selected as
[MFEMTimeDerivativeMassKernel](https://mooseframework.inl.gov/source/mfem/kernels/MFEMTimeDerivativeMassKernel.html).
This adds the weak form of $\rho c_p\dot{u}$ for our solver to work with,
where the coefficient $\rho c_p$ is set to 1 by default.

!listing test/tests/mfem/kernels/heattransfer.i start=[Kernels] end=[BCs]

### Convective Boundary Condition

For this problem we have defined an extra boundary condition
to represent heat exchange between solid surface and a surrounding fluid.
This is governed by the following equation

!equation id=convection equation
q_\text{n} = h(T - T_\infty)

Where $q_\text{n}$ is the outward heat flux from the surface,
$h$ is the convective heat transfer coefficient,
and $T_\infty$ is the surrounding temperature.

!listing test/tests/mfem/kernels/heattransfer.i block=BCs

### Assembly Level

The transience in our equation requires the executioner to solve across time.
This is trivially done by selecting the
[`MFEMTransient`](https://mooseframework.inl.gov/source/mfem/executioners/MFEMTransient.html)
executioner and setting a time step, start time, and end time.

Here we expose the `assembly_level` setting which defaults to `legacy`.
This tells MFEM how to assemble the
[FEM operator matrices](https://mfem.org/howto/assembly_levels/).
This can be helpful for advanced users
and allows for optmisation decisions balancing memory usage vs compute
along side archetectural decisions
(e.g. use of gpus vs cpus as they benifit from different matrix assemblies).
For now we will not utilize this setting
but it is useful to have an awareness off.

!listing test/tests/mfem/kernels/heattransfer.i block=Executioner

### Inactive Auxiliary Variable

Also present in this example is the inactive auxiliary variable
