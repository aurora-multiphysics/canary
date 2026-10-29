# MFEM Example 34 (examples/ex34p.cpp): electric potential
#
# Sub-application of current_density.i. Solves -div(grad(phi)) = 0 on the
# conductor, a submesh of the fichera mesh, with phi = 0 on boundary 2 and
# phi = 1 on boundary 23, using CG preconditioned by BoomerAMG as in ex34p.

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
[]

[Variables]
  [phi]
    type = MFEMVariable
    fespace = ConductorH1FESpace
  []
[]

[BCs]
  [phi_0]
    type = MFEMScalarDirichletBC
    variable = phi
    boundary = 2
    coefficient = 0
  []
  [phi_1]
    type = MFEMScalarDirichletBC
    variable = phi
    boundary = 23
    coefficient = 1
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = phi
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = -1
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = boomeramg
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
