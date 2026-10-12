# Definite Maxwell problem solved with Nedelec elements of the first kind
# based on MFEM Example 3.

omega=${fparse 2.0*3.14159265358979323846*50.0}  # angular frequency 2*PI
sigma_vac = 0
sigma_coil = 0
sigma_target = 0.3278e8 # Siemens per meter (S/m) of the conductivity plate
nu=795774.715 #  (meters/Henry) = 1/magentic permiablity of free space
epsilon= 8.8541878176e-12 #Farads/m of free space

[Mesh]
  type = MFEMFileMesh
  file = ./team_coil_two_vols_plate_exterior_tet_m_offset_box.msh
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[SubMeshes]
  [coil_complement]
    type = MFEMDomainSubMesh
    block = 'Free_space Plate'
  []
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
  [SubmeshHCurlFESpace]
      type = MFEMVectorFESpace
      fec_type = ND
      fec_order = FIRST
      submesh = coil_complement
  []
  [SubmeshHDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
    submesh = coil_complement
  []
[]

[Variables]
  [a_field]
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [b_field]
    type = MFEMComplexVariable
    fespace = HDivFESpace
  []
  [source_j_field]
    type = MFEMComplexVariable
    fespace = HDivFESpace
  []
[]

[AuxKernels]
  [curl]
    type = MFEMComplexCurlAux
    variable = b_field
    source = a_field
    execute_on = TIMESTEP_END
  []
[]

[BCs]
    # A = iE/w on coil surface
    # [coil_surface_a_field]
    #     type = MFEMComplexVectorTangentialDirichletBC
    #     variable = a_field
    #     vector_coefficient_real = 0.0
    #     vector_coefficient_imag = 0.0
    #     boundary = 7
    # []
    [tangential_a_bdr]
        type = MFEMComplexVectorTangentialDirichletBC
        variable = a_field
        boundary = 'Boundary' #free space
    []
[]

[Constraints]
  [tree_cotree_gauge]
    type = MFEMComplexTreeCotreeGaugeEssentialConstraint
    variable = a_field
    # Gauge only the non-conducting region: block 3 is the vacuum surrounding the
    # TorusCore/TorusSheath conductors, where the sigma*dA/dt term already fixes
    # the gauge. Conductor edges seed the spanning forest but are not gauged.
    block = 'Free_space First_half Second_half'
    # Boundaries where a tangential Dirichlet condition is applied to a_field, so
    # the interior gauge is seeded consistently with that boundary condition.
    boundary = 'Boundary'
  []
[]

[Functions]
    # (i * \omega * \sigma - \omega^2 * \epsilon0)* A represented as (massCoef + i*loss_coef)*A
    # where massCoef = -omega^2 * epsilon0, lossCoef = \omega * sigma
    [mass_coef]
        type = ParsedFunction
        expression = -${epsilon}*${omega}^2
    []
    [loss_coef_vac]
        type = ParsedFunction
        expression = ${omega}*${sigma_vac}
    []
    [loss_coef_coil]
        type = ParsedFunction
        expression = ${omega}*${sigma_coil}
    []
    [loss_coef_target]
        type = ParsedFunction
        expression = ${omega}*${sigma_target}
    []
[]

[FunctorMaterials]
    #expose \sigma, nu, mass/loss for j*\omega*\sigma
    [vacuum]
        type = MFEMGenericFunctorMaterial
        prop_names = 'massCoef lossCoef sigma nu'
        prop_values = 'mass_coef loss_coef_vac ${sigma_vac} ${nu}'
        block = 'Free_space'
    []
    [coil]
        type = MFEMGenericFunctorMaterial
        prop_names = 'massCoef lossCoef sigma nu'
        prop_values = 'mass_coef loss_coef_coil ${sigma_coil} ${nu}'
        block = 'First_half Second_half'
    []
    [target]
        type = MFEMGenericFunctorMaterial
        prop_names = 'massCoef lossCoef sigma nu '
        prop_values = 'mass_coef loss_coef_target ${sigma_target} ${nu} '
        block = 'Plate'
    []
[]

[Kernels]
  [curlcurl]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMCurlCurlKernel
      coefficient = ${nu}
    []
  []
  [conductive_mass_complex_plate_freespace]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMVectorFEMassKernel
      coefficient = massCoef # = - (omega**2)*epsilon
    []
  []
  [conductive_mass_complex_plate]
    type = MFEMComplexKernel
    variable = a_field
    block = 'Plate'
    [ImagComponent]
      type = MFEMVectorFEMassKernel
      coefficient = loss_coef_target
    []
  []
  [source]
    type = MFEMComplexKernel
    variable = a_field
    block = 'First_half Second_half'
    [RealComponent]
      type = MFEMVectorFEDomainLFKernel
      vector_coefficient = source_j_field_real # = J
    []
  []
[]

[Solvers]
  [superlu]
    type = MFEMMUMPS
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[VectorPostprocessors]
   [line_sample]
     type = MFEMComplexVariableLineValueSampler
     variable = 'b_field'
     start_point = '0.03 0.0 0.00685'
     end_point = '0.03 0.11 0.00685'
     num_points = 20
   []
 []

[MultiApps]
  [subapp]
    type = FullSolveMultiApp
    input_files = closed_coil_submeshlaplacesolve.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_sub_source_j_field]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = source_j_field
    variables = source_j_field
    from_multi_app = subapp
  []
[]

[Outputs]
 [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/TEAM3CSV
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TEAM3
  []
[]
