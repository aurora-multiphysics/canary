# MFEM Example 27 (examples/ex27p.cpp)
#
# Solves the Laplace problem -div(grad(u)) = 0 on a rectangle with two circular
# holes, with a different boundary condition on each part of the boundary:
#
#   attribute 1 (lower edge):  Neumann,   n.grad(u) = 1
#   attribute 2 (upper edge):  Robin,     n.grad(u) + u = 1
#   attribute 3 (left hole):   Dirichlet, u = 0
#   attribute 4 (right hole):  natural,   n.grad(u) = 0
#
# The narrow ends of the rectangle are identified with each other, so the
# solution is periodic in x. ex27p builds this mesh in code, refines it twice
# and bends the hole boundaries onto circles; periodic-rectangle-two-holes.mesh
# is that mesh, written by ex27p itself. As in ex27p, it is refined once more in
# parallel.
# The system is solved with first order Lagrange elements and hypre's PCG
# preconditioned by BoomerAMG, and the L2 norm of the solution is compared
# against the value computed by ex27p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/periodic-rectangle-two-holes.mesh
  parallel_refine = 1
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [neumann]
    type = MFEMBoundaryIntegratedBC
    variable = u
    boundary = 1
    coefficient = 1
  []
  [robin]
    type = MFEMConvectiveHeatFluxBC
    variable = u
    boundary = 2
    heat_transfer_coefficient = 1
    T_infinity = 1
  []
  [dirichlet]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 3
    coefficient = 0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 1
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 200
    print_level = 2
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [solution_l2_norm]
    type = MFEML2Error
    variable = u
    function = 0
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex27p
  []
[]
