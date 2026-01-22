# Solve TEAM Problem 4 for the induced current in the FELIX brick due to a decaying background field

conductor_domains = 'Brick'
vacuum_domain = 'Vacuum'
vacuum_resistivity = 1.0 # S/m
vacuum_permeability = '${fparse (4*pi*1e-7)}' # T m/A
conductor_conductivity = 2.538e7 # S/m 
conductor_resistivity = '${fparse 1.0/(conductor_conductivity)}'

tangential_induced_h_boundaries = '1 2 5 6' # H x n = 0, no field from current in system, so H_ind = 0

[Problem]
  type = MFEMProblem
[]

[Mesh]
  type = MFEMMesh
    file = ./team4_symmetrized.e
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[Variables]
  [induced_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [external_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [j_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[AuxKernels]
  [update_j_field]
    type = MFEMCurlAux
    variable = j_field
    source = induced_h_field
    scale_factor = 1.0
    execute_on = TIMESTEP_END
    execution_order_group = 2
  []
[]

[Functions]
  [zero_vector]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = '0'
  []
  [dbext_dt_magnitude]
    type = ParsedFunction
    expression = (-mu0/tau)*exp(-t/tau)
    symbol_names = 'mu0 tau'
    symbol_values = '${vacuum_permeability} 0.0119'
  [] 
[]

[FunctorMaterials]
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity permeability'
    prop_values = '${conductor_resistivity} ${vacuum_permeability}'
    block = ${conductor_domains}
  []
  [Vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity permeability'
    prop_values = '${vacuum_resistivity} ${vacuum_permeability}'
    block = ${vacuum_domain}
  []
[]

[BCs]
  [tangential_induced_h_bc]
    type = MFEMVectorTangentialDirichletBC
    variable = induced_h_field
    vector_coefficient = zero_vector
    boundary = ${tangential_induced_h_boundaries}
  []
[]

[Kernels]
  [mu.dhind_dt,hind'] 
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = induced_h_field
    coefficient = permeability
  []
  [mu.dhext_dt,hind'] 
    type = MFEMMixedVectorMassKernel
    variable = induced_h_field
    trial_variable = external_h_field
    coefficient = dbext_dt_magnitude
  []
  [rho.curlhind,curlh']
    type = MFEMCurlCurlKernel
    variable = induced_h_field
    coefficient = resistivity
    # block = ${conductor_domains}
  []
[]

[Preconditioner]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
[]

[Solver]
  type = MFEMHyprePCG
  preconditioner = ams
  l_tol = 1e-16
  l_abs_tol = 1e-16
  l_max_its = 500
[]

[Executioner]
  type = MFEMTransient
  dt = 0.005
  start_time = 0.0
  end_time = 0.02
[]

[MultiApps]
  [source_divergence_cleaning]
    type = FullSolveMultiApp
    input_files = team4_external_source_field.i
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Transfers]
  [from_external_field]
    type = MultiAppMFEMCopyTransfer
    source_variable = external_h_field
    variable = external_h_field
    from_multi_app = source_divergence_cleaning
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/TEAM4CSV
  []
  [VacuumParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TEAM4
  []
[]
