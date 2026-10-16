# MFEM Example 35 (examples/ex35p.cpp)
#
# Solves the damped harmonic oscillator -div(grad(u)) - omega^2 u + i omega sigma u = 0
# for a complex scalar field u on the fichera mesh, with omega = 2 pi and
# sigma = 2. The oscillation is driven through the boundary condition: on the
# port, boundaries 7, 8, 11 and 12, u is the second eigenmode of the Laplacian
# on the port, computed by the sub-application port_mode.i; on the rest of the
# boundary u = 0. As in ex35p, the mesh is refined once in serial and once in
# parallel and the problem is discretized with first order Lagrange elements.
# ex35p solves the complex system with FGMRES to a relative tolerance of 1e-6,
# preconditioned block by block with BoomerAMG on a separate real operator,
# which MOOSE cannot form; with -slu it solves directly with SuperLU, as here.
# The L2 norm of the complex solution is compared against the value computed by
# ex35p -slu.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/fichera-mixed.mesh
  serial_refine = 1
  parallel_refine = 1
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
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
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # A complex problem cannot hold real variables, so the real port mode is received as the real
  # part of a complex field
  [port_bc]
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [port]
    type = MFEMComplexScalarDirichletBC
    variable = u
    coefficient_real = port_bc_real
    coefficient_imag = 0
  []
[]

[Kernels]
  # -omega^2 = -(2 pi)^2 and omega sigma = 4 pi, written out to full precision
  [stiffness]
    type = MFEMComplexKernel
    variable = u
    [RealComponent]
      type = MFEMDiffusionKernel
      coefficient = 1
    []
    [ImagComponent]
      type = MFEMDiffusionKernel
      coefficient = 0
    []
  []
  [mass_and_loss]
    type = MFEMComplexKernel
    variable = u
    [RealComponent]
      type = MFEMMassKernel
      coefficient = -39.47841760435743
    []
    [ImagComponent]
      type = MFEMMassKernel
      coefficient = 12.566370614359172
    []
  []
[]

[Solvers]
  [direct]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [port]
    type = FullSolveMultiApp
    input_files = port_mode.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_port]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = port_bc
    variables = port_bc
    from_multi_app = port
  []
[]

[Postprocessors]
  [solution_l2_norm]
    type = MFEMComplexL2Error
    variable = u
    function_real = 0
    function_imag = 0
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex35p
  []
[]
