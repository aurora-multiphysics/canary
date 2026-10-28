# MFEM Example 3 (examples/ex3p.cpp)
#
# Solves the definite Maxwell problem curl curl E + E = f on the beam mesh, with
# the right hand side f chosen so that the exact solution is
# E = (sin(kappa y), sin(kappa z), sin(kappa x)) with kappa = pi. The tangential
# component of E is set to that of the exact solution on all boundaries. As in
# ex3p, the mesh is refined serially as many times as keeps it at or under 1,000
# elements (once for the 48 elements of the beam mesh), then twice more in
# parallel. The solution is found with first order Nedelec elements using
# hypre's PCG with an AMS preconditioner, matching the settings of ex3p. The L2
# error of the solution is compared against the value computed by ex3p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/beam-tet.mesh
  serial_refine = 1
  parallel_refine = 2
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
[]

[Variables]
  [E]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[Functions]
  [E_exact]
    type = ParsedVectorFunction
    expression_x = 'sin(kappa * y)'
    expression_y = 'sin(kappa * z)'
    expression_z = 'sin(kappa * x)'
    symbol_names = kappa
    symbol_values = 3.141592653589793
  []
  [f_exact]
    type = ParsedVectorFunction
    expression_x = '(1 + kappa * kappa) * sin(kappa * y)'
    expression_y = '(1 + kappa * kappa) * sin(kappa * z)'
    expression_z = '(1 + kappa * kappa) * sin(kappa * x)'
    symbol_names = kappa
    symbol_values = 3.141592653589793
  []
[]

[BCs]
  [tangential_E]
    type = MFEMVectorTangentialDirichletBC
    variable = E
    vector_coefficient = E_exact
  []
[]

[Kernels]
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = E
  []
  [mass]
    type = MFEMVectorFEMassKernel
    variable = E
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = E
    vector_coefficient = f_exact
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    print_level = 1
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = ams
    # Start from the boundary values, as mfem::HyprePCG does by default in ex3p
    use_initial_guess = true
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 500
    print_level = 2
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [l2_error]
    type = MFEMVectorL2Error
    variable = E
    function = E_exact
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex3p
  []
[]
