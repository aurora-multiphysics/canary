# Example 3: Linear Elastic Problem

Here we will talk through the MFEM
[Linear Elasticity Problem](https://mooseframework.inl.gov/syntax/MFEM/LinearElasticity.html)
which is based on [MFEM example 2](https://mfem.org/examples/).
This is to highlight benefits of using different solvers and explain some more
input file syntax. This problem also introduces how to define a mesh over
multiple materials with different parameters.

## Input File

### Mesh and Boundaries

The problem setup is a bar composed of of two imaginary materials Rigidium and
Bendium with self-explanatory material properties. The bar will be fixed at one
end face and a uniform force applied across the opposite end face.

To define this problem in a moose input file we start by familiarly loading a
mesh file and applying relevant boundary conditions. This comes with boundary
labels defined in the `.mesh` file, on which we will define the following
boundary conditions:

| Boundary label | Location              | Boundary condition in this example              |
| -------------- | --------------------- | ----------------------------------------------- |
| `1`            | End face at $x = 0$   | Fixed displacement: $u_x = u_y = u_z = 0$       |
| `2`            | End face at $x = 8$   | Uniform load in $-z$ direction: $(0, 0, -0.01)$ |
| `3`            | Remaining outer faces | Traction-free (the natural boundary condition)  |

Some settings under the `Mesh` block are so far unseen in these examples.

- `uniform_refine` defines how many refining passes are made to subdivide each
  element into a smaller element. This allows for the user to select a finer
  mesh without regenerating the mesh file.
- `displacement` defines which variable describes the mesh's deformation so
  that the mesh can deform appropriately.

!listing test/tests/mfem/kernels/linearelasticity.i block=Mesh BCs

### Functor Materials

A functor in programming is an object you can call like a function. In this
MOOSE input file they are used to provide values to our solver when given a
location and time state. We use them within the
[`FunctorMaterials`](https://mooseframework.inl.gov/syntax/FunctorMaterials/)
block to describe our two material properties that apply to different sections
of the mesh. These properties are defined via three settings:

- `prop_names` labels the properties being created
- `prop_values` assigns a value to the properties named
- `block` assigns the defined properties to a specific mesh block

!listing test/tests/mfem/kernels/linearelasticity.i block=FunctorMaterials

### Finite Element Space Settings

There are two new settings within the `FESpaces` block that previous examples
have not exposed.

- `range_dim` sets the number of components expected in a vector field. Here
  we have three dimensions to define our vectors in.
- `ordering` controls how the solver arranges the solution vectors. This is
  an advanced option used to specify the
  [major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
  of the vectors. For the most part the default will suffice but check which
  ordering the preconditioning solver expects.

!listing test/tests/mfem/kernels/linearelasticity.i block=FESpaces

### Selecting the Physics Kernel and Solver

Luckily for us there is the existing physics kernel
[`MFEMLinearElasticityKernel`](https://mooseframework.inl.gov/source/mfem/kernels/MFEMLinearElasticityKernel.html)
for solving linear elastic problems with MFEM spaces. This does require the
extra definition of two material properties `lambda` and `mu` which here we set
to the corresponding functor material property names.

!listing test/tests/mfem/kernels/linearelasticity.i block=Kernels

The `Solvers` block here preserves the last example's
[BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html)
preconditioner but swaps the `main` solver for
[`MFEMHyperPCG`](https://mooseframework.inl.gov/source/mfem/solvers/MFEMHyprePCG.html).
While a generalized minimal residual method (GMRES) approach used previously
would still solve this example, a preconditioned conjugate gradient (PCG)
approach is better suited to this problem.

With the positive elastic constants in this example and one end clamped, the
stiffness matrix is symmetric positive definite. Positive definite means every
displacement allowed by the beam stores positive strain energy. For this beam,
that holds because:

- Both materials have elastic constants that give positive energy when strained.
- Clamping one end prevents rigid body translations and rotations, which would
  otherwise allow displacement without storing energy.

PCG solvers exploit this structure and allows for less computational overhead
than the more general GMRES solver.

!listing test/tests/mfem/kernels/linearelasticity.i block=Solvers

### Print Level

Also set within the `Solvers` block is the `print_level` option. This dictates
what information about the simulation is output to the terminal as it runs.
The meaning of each option is as follows:

| print_level | BoomerAMG output                 |
| ----------- | -------------------------------- |
| `0`         | None                             |
| `1`         | Setup information                |
| `2`         | Solve information                |
| `3`         | Both setup and solve information |
