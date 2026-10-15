# MFEM Example 14 (examples/ex14p.cpp)
#
# Solves the Poisson problem -div(grad(u)) = 1 on the star mesh with u = 0 on all
# boundaries, using a discontinuous Galerkin discretization. As in ex14p, the
# mesh is refined serially as many times as keeps it at or under 10,000 elements
# (four times for the 20 elements of the star mesh), then twice more in
# parallel. The solution is found with first order discontinuous elements, using
# the symmetric interior penalty method with the penalty parameter of ex14p. The
# Dirichlet condition is imposed weakly through boundary face terms. The system
# is solved with CG preconditioned by BoomerAMG, and the L2 norm of the solution
# is compared against the value computed by ex14p.

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
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = FIRST
    basis = GaussLegendre
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[BCs]
  [dirichlet]
    type = MFEMDGDiffusionBC
    variable = u
    boundary_coefficient = 0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
  []
  [dg_diffusion]
    type = MFEMDGDiffusionKernel
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
    l_tol = 0
    l_max_its = 1
    print_level = 1
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 500
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
    file_base = OutputData/ex14p
  []
[]
