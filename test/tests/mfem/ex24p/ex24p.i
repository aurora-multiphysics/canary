# MFEM Example 24 (examples/ex24p.cpp)
#
# Projects the gradient of p = sin(x) sin(y) sin(z), interpolated in a first
# order H1 space, onto first order Nedelec elements on the beam-hex mesh: E in
# H(curl) solves (E, v) = (grad p, v) for all v in H(curl). As in ex24p
# (problem 0), the mesh is refined serially as many times as keeps it at or under
# 1,000 elements (twice for the 8 elements of beam-hex), then once more in
# parallel. The gradient is also computed with the discrete gradient operator,
# and the exact gradient is projected onto the Nedelec space. The L2 errors of
# all three against the exact gradient are compared against the values computed
# by ex24p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/beam-hex.mesh
  serial_refine = 2
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

[AuxVariables]
  [p]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [E_interpolant]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [E_projection]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[Functions]
  [p_exact]
    type = ParsedFunction
    expression = 'sin(x) * sin(y) * sin(z)'
  []
  [grad_p_exact]
    type = ParsedVectorFunction
    expression_x = 'cos(x) * sin(y) * sin(z)'
    expression_y = 'sin(x) * cos(y) * sin(z)'
    expression_z = 'sin(x) * sin(y) * cos(z)'
  []
[]

[ICs]
  [p_ic]
    type = MFEMScalarIC
    variable = p
    coefficient = p_exact
  []
  [E_projection_ic]
    type = MFEMVectorIC
    variable = E_projection
    vector_coefficient = grad_p_exact
  []
[]

[Kernels]
  [mass]
    type = MFEMVectorFEMassKernel
    variable = E
  []
  [gradient]
    type = MFEMMixedVectorGradientKernel
    variable = E
    trial_variable = p
    coefficient = -1
  []
[]

[AuxKernels]
  [discrete_gradient]
    type = MFEMGradAux
    variable = E_interpolant
    source = p
  []
[]

[Solvers]
  [jacobi]
    type = MFEMOperatorJacobiSmoother
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = jacobi
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 1000
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [solution_l2_error]
    type = MFEMVectorL2Error
    variable = E
    function = grad_p_exact
  []
  [interpolant_l2_error]
    type = MFEMVectorL2Error
    variable = E_interpolant
    function = grad_p_exact
  []
  [projection_l2_error]
    type = MFEMVectorL2Error
    variable = E_projection
    function = grad_p_exact
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex24p
  []
[]
