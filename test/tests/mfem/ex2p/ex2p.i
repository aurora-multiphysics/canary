# MFEM Example 2 (examples/ex2p.cpp)
#
# Solves a linear elasticity problem, -div(sigma(u)) = 0, for a cantilever beam
# made of two materials. The beam is fixed (u = 0) on boundary 1 and pulled down
# by a constant traction on boundary 2. As in ex2p, the mesh is refined serially
# as many times as keeps it at or under 1,000 elements (twice for the 16
# elements of the beam mesh), then once more in parallel. The solution is found
# with first order vector H1 elements using hypre's PCG with a BoomerAMG
# preconditioner set up for elasticity, matching the settings of ex2p run with
# -elast. The L2 norm of the displacement is compared against the value
# computed by ex2p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/beam-tri.mesh
  serial_refine = 2
  parallel_refine = 1
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
    ordering = VDIM
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [fixed]
    type = MFEMVectorDirichletBC
    variable = displacement
    boundary = 1
    vector_coefficient = '0 0'
  []
  [pull_down]
    type = MFEMVectorBoundaryIntegratedBC
    variable = displacement
    boundary = 2
    vector_coefficient = '0 -0.01'
  []
[]

[FunctorMaterials]
  [material_1]
    type = MFEMGenericFunctorMaterial
    prop_names = 'lambda mu'
    prop_values = '50 50'
    block = 1
  []
  [material_2]
    type = MFEMGenericFunctorMaterial
    prop_names = 'lambda mu'
    prop_values = '1 1'
    block = 2
  []
[]

[Kernels]
  [elasticity]
    type = MFEMLinearElasticityKernel
    variable = displacement
    lambda = lambda
    mu = mu
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    fespace = H1FESpace
    l_tol = 0
    l_max_its = 1
    print_level = 1
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-8
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
  [solution_l2_norm]
    type = MFEMVectorL2Error
    variable = displacement
    function = '0 0'
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex2p
  []
[]
