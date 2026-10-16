# MFEM Example 34 (examples/ex34p.cpp): current density
#
# Sub-application of ex34p.i. Receives the electric potential phi on the
# conductor from electric_potential.i, and solves for the current density
# J = -grad(phi) in H(div) on the conductor, with J.n = 0 on its walls: the
# symmetry planes (boundaries 9 to 16) and its interface with the rest of the
# domain (boundary 25, created with the submesh). The mass system is solved with
# CG preconditioned by Jacobi; ex34p uses hypre's diagonal scaling, which MOOSE
# does not have. J is then transferred to an H(div) field on the full mesh,
# zero outside the conductor.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/fichera-mixed-conductor.mesh
  serial_refine = 1
  parallel_refine = 1
[]

[Problem]
  type = MFEMProblem
[]

[SubMeshes]
  [conductor]
    type = MFEMDomainSubMesh
    block = 2
  []
[]

[FESpaces]
  [ConductorH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = conductor
  []
  [ConductorHDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
    submesh = conductor
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[Variables]
  [conductor_J]
    type = MFEMVariable
    fespace = ConductorHDivFESpace
  []
[]

[AuxVariables]
  [phi]
    type = MFEMVariable
    fespace = ConductorH1FESpace
  []
  [J]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[BCs]
  [zero_normal_current]
    type = MFEMVectorNormalDirichletBC
    variable = conductor_J
    boundary = '9 10 11 12 13 14 15 16 25'
    vector_coefficient = '0 0 0'
  []
[]

[Kernels]
  # (J, v) + (grad(phi), v) = 0, with the known phi moved to the right-hand side
  [mass]
    type = MFEMVectorFEMassKernel
    variable = conductor_J
  []
  [gradient]
    type = MFEMMixedVectorGradientKernel
    variable = conductor_J
    trial_variable = phi
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
    l_max_its = 2000
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [potential]
    type = FullSolveMultiApp
    input_files = electric_potential.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_potential]
    type = MultiAppMFEMCopyTransfer
    source_variables = phi
    variables = phi
    from_multi_app = potential
  []
  [to_full_mesh]
    type = MFEMSubMeshTransfer
    from_variable = conductor_J
    to_variable = J
    execute_on = TIMESTEP_END
  []
[]
