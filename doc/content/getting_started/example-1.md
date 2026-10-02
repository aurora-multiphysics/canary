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

### Mesh

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

### Variables

The variables block of our input file
demonstrates this use of sensible defaults
only specifying the variable symbol `u`.
MAKE THIS BETTER THIS IS NOT EXPLICIT ENOUGH.

!listing canary/test/tests/kernels/simple_diffusion/simple_diffusion.i block=Variables

### Kernels

### Boundary Conditions (BCs)

### Executioner
