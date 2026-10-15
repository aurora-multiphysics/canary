# MFEM Example 4 (examples/ex4p.cpp)
#
# Solves the H(div) diffusion problem -grad(div F) + F = f on the star mesh,
# with the right hand side f chosen so that the exact solution is
# F = (cos(kappa x) sin(kappa y), cos(kappa y) sin(kappa x)) with kappa = pi.
# The normal component of F is set to that of the exact solution on all
# boundaries. As in ex4p, the mesh is refined serially as many times as keeps it
# at or under 1,000 elements (twice for the 20 elements of the star mesh), then
# twice more in parallel. The solution is found with lowest order
# Raviart-Thomas elements using PCG with an AMS preconditioner, matching the
# settings of ex4p in two dimensions. The L2 error of the solution is compared
# against the value computed by ex4p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/star.mesh
  serial_refine = 2
  parallel_refine = 2
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[Variables]
  [F]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[Functions]
  [F_exact]
    type = ParsedVectorFunction
    expression_x = 'cos(kappa * x) * sin(kappa * y)'
    expression_y = 'cos(kappa * y) * sin(kappa * x)'
    symbol_names = kappa
    symbol_values = 3.141592653589793
  []
  [f_exact]
    type = ParsedVectorFunction
    expression_x = '(1 + 2 * kappa * kappa) * cos(kappa * x) * sin(kappa * y)'
    expression_y = '(1 + 2 * kappa * kappa) * cos(kappa * y) * sin(kappa * x)'
    symbol_names = kappa
    symbol_values = 3.141592653589793
  []
[]

[BCs]
  [normal_F]
    type = MFEMVectorNormalDirichletBC
    variable = F
    vector_coefficient = F_exact
  []
[]

[Kernels]
  [divdiv]
    type = MFEMDivDivKernel
    variable = F
  []
  [mass]
    type = MFEMVectorFEMassKernel
    variable = F
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = F
    vector_coefficient = f_exact
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HDivFESpace
    print_level = 1
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = ams
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
  [l2_error]
    type = MFEMVectorL2Error
    variable = F
    function = F_exact
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex4p
  []
[]
