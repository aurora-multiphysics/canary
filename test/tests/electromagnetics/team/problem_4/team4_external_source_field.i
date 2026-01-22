# Solve for the magnetic field around a closed conductor subject to
# global current constraint.

normal_b_boundaries = '1 2 5 6' 
vacuum_permeability = '${fparse (4*pi*1e-7)}'

[Problem]
  type = MFEMProblem
[]

[Mesh]
    type = MFEMMesh
    file = ./team4_symmetrized.e
[]

[Functions]
  # Here, externally applied B field = grad (magnetic_potential)
  [boundary_magnetic_potential]
    type = ParsedFunction
    expression = z
  []
[]

[FunctorMaterials]
  [Domain]
    type = MFEMGenericFunctorMaterial
    prop_names = permeability
    prop_values = ${vacuum_permeability}
  []
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
  [magnetic_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [external_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  [update_external_h_field]
    type = MFEMGradAux
    variable = external_h_field
    source = magnetic_potential
    scale_factor = '${fparse (0.1/vacuum_permeability)}' # -B0 * mu
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  # Set zero of magnetic potential on symmetry plane
  [MagneticInsulatingBoundaries]
    type = MFEMScalarDirichletBC
    variable = magnetic_potential
    boundary = ${normal_b_boundaries}
    coefficient = boundary_magnetic_potential
  []
[]

[Kernels]
  [mu.gradphi,gradphi]
    type = MFEMDiffusionKernel
    variable = magnetic_potential
    coefficient = permeability
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = boomeramg
  l_tol = 1e-8
  l_max_its = 100
[]

[Executioner]
  type = MFEMSteady
[]

[Postprocessors]
  [MagneticEnergy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = ${fparse 0.5*vacuum_permeability}
    dual_variable = external_h_field
    primal_variable = external_h_field
    execution_order_group = 1
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/HPhiMagnetostaticClosedCoilCSV
  []
  [VacuumParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/HPhiMagnetostaticClosedCoil
    vtk_format = ASCII
  []
[]
