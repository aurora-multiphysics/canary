# Example 1: Simple Diffusion

To start let us run a basic steady-state diffusion example.
We will aim to solve the diffusion equation

\begin{equation}
-\nabla^{2}u=0 \label{diffusion_eq}
\end{equation}

over a 2D square domain, where $u$ is the unknown scalar field.
MOOSE allows us to solve this differential equations
if the problem's behavior is fully captured by a MOOSE input file.

## Writing an input file

These input files are formed of six sections:

- Mesh
- Variables
- Kernels
- Boundary Conditions
- Executioner
- Outputs

### [Mesh](https://mooseframework.inl.gov/syntax/Mesh/)

Firstly, we need to mesh over our domain.
This can be done either with an externally generated mesh or,
for simple meshes,
we can use the built in [`GenerateMesh`](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html) object.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Mesh

Here we specify our mesh's dimension with
[`dim`](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html)
along with the number of elements to generate each direction
[`nx`](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html)
and
[`ny`](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html).
If you clicked through to the
[docs page](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html)
for
[`GenerateMesh`](https://mooseframework.inl.gov/moose/source/mesh/GeneratedMesh.html)
you would notice that we are specifying only a fraction of the parameters.
MOOSE allows users to specify the minimum parameters to describe an object
with the remaining being optional set as sensible defaults.

### [Variables](https://mooseframework.inl.gov/moose/syntax/Variables/index.html)

The variables block of our input file
allows us to specify what the unknowns of our problem.
Often referred to as nonlinear variables in moose documentation,
the variables block allows us to specify the order and the family
of the
[shape functions](https://mooseframework.inl.gov/source/problems/FEProblemBase.html#shape_functions)
used to solve for the variable.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Variables

Here we have kept it simple with a first order or linear shape function,
and specified that they should be from the lagrangian family.
For a list of available shape function orders and families
see the
[variable system](https://mooseframework.inl.gov/moose/syntax/Variables/index.html)
docs.

### [Kernels](https://mooseframework.inl.gov/syntax/Kernels/index.html)

In computing kernel is an overloaded term with many definitions.
MOOSE uses a kernel to mean a piece of physics.
They represent terms in the weak form
of our partial differential equations.
These can be selected from a large selection of premade kernels,
or with some C++ knowledge [created custom](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/examples/ex02_kernel.html).

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Kernels

Here we are solving a diffusion problem
so we'll call our physics kernel `diff`.
There is an easy choice for the type
as we can utilize the kernel type
[`Diffusion`](https://mooseframework.inl.gov/source/kernels/Diffusion.html)
provided by MOOSE.
We then specify the variable we want this kernel to act on as `u`.

### [Boundary Conditions (BCs)](https://mooseframework.inl.gov/syntax/BCs/index.html)

To keep our example simple we will set `u` to 0 at the left boundary and 1 at the right boundary of of our 2D square mesh.
We can do so in our input file by simply specifying left and right
boundaries in the `BCs` section of the MOOSE input file.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=BCs

The type to select for a constant value boundary condition is
[`DirichletBC`](https://mooseframework.inl.gov/source/bcs/DirichletBC.html).
We must also specify the variable the boundary will apply to
along with its value.
You might wonder what happens to the top and bottom boundaries
which are unspecified here.
These are set to natural (zero flux)
Neumann boundary conditions by default.
To see a full list of supported BCs see [this page](https://mooseframework.inl.gov/syntax/BCs/index.html).

### [Executioner](https://mooseframework.inl.gov/syntax/Executioner/index.html)

Simply put the executioner block defines how the code is executed.
Here we can control the solving method,
approach to time stepping,
and or control [PETSc](https://petsc.org/release/documentation/) behavior.
PETSc is used as the underlying solver in MOOSE
and its options can be specified here.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Executioner

### Outputs

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Outputs
