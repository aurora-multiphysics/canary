# MFEM Example 1 (examples/ex1p.cpp)
#
# Solves the Poisson problem -\nabla^2 u = 1 on the star mesh with zero
# Dirichlet conditions on all external boundaries. As in ex1p, the mesh is
# refined serially as many times as keeps it at or under 10,000 elements (four
# times for the 20 elements of the star mesh), then twice more in parallel. The
# solution is found with first order H1 elements using PCG with a BoomerAMG
# preconditioner, matching the settings of ex1p. The L2 norm of the solution
# is compared against the value computed by ex1p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/star.mesh
  serial_refine = 4
  parallel_refine = 2
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
  [dirichlet]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = 0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
  []
  [source]
    type = MFEMDomainLFKernel
    variable = u
    coefficient = 1
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 1
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 2000
    print_level = 1
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
    file_base = OutputData/ex1p
  []
[]
