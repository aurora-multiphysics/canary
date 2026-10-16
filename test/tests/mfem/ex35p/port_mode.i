# MFEM Example 35 (examples/ex35p.cpp): port boundary condition
#
# Sub-application of ex35p.i. Computes the five lowest eigenmodes of the
# Laplacian, -div(grad(phi)) = lambda phi with phi = 0 on its edges, on the port:
# a submesh made of boundaries 7, 8, 11 and 12 of the fichera mesh. As in ex35p,
# the eigenmodes are found with LOBPCG preconditioned by BoomerAMG, with random
# seed 75. The second mode (index 1) is transferred to a field on the full mesh,
# zero away from the port, for use as the boundary condition in ex35p.i.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/fichera-mixed.mesh
  serial_refine = 1
  parallel_refine = 1
[]

[Problem]
  type = MFEMEigenproblem
  num_modes = 5
[]

[SubMeshes]
  [port]
    type = MFEMBoundarySubMesh
    boundary = '7 8 11 12'
  []
[]

[FESpaces]
  [PortH1FESpace]
    type = MFEMGenericFESpace
    fec_name = H1_2D_P1
    submesh = port
  []
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [port_mode]
    type = MFEMVariable
    fespace = PortH1FESpace
  []
[]

[AuxVariables]
  [port_bc]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [port_edges]
    type = MFEMScalarDirichletBC
    variable = port_mode
    coefficient = 0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = port_mode
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    l_tol = 0
    l_max_its = 1
    print_level = 0
  []
  [lobpcg]
    type = MFEMHypreLOBPCG
    preconditioner = boomeramg
    random_seed = 75
    l_tol = 1e-8
    l_max_its = 200
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Transfers]
  [to_full_mesh]
    type = MFEMSubMeshTransfer
    from_variable = port_mode_1
    to_variable = port_bc
    execute_on = TIMESTEP_END
  []
[]
