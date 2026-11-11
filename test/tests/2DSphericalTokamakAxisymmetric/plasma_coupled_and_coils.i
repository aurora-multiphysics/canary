[Mesh]
  type = MFEMFileMesh
  file = corrected_pf_mesh_quad4_2d_bc12.e
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

  [RTFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []

  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
    basis = GaussLegendre
  []
[]

[Variables]
  [Atheta]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # PF current density for visualization
  [Jtheta_pf]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # Spatial plasma current density for visualization
  [Jtheta_plasma]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # PF source mu0 * R * Jtheta
  [source_pf]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # Plasma source mu0 * R * Jtheta
  [source_plasma]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # Magnetic field
  [B]
    type = MFEMVariable
    fespace = RTFESpace
  []
[]

[Functions]
  [zero]
    type = MFEMParsedFunction
    expression = 0
  []

  # ======================================================
  # Axisymmetric coordinate system
  #
  # x = R
  # y = Z
  # ======================================================

  [cylindrical]
    type = MFEMCoordinateTransformations
    coord_type = RZ
    inv_r_eps = 1e-6
  []


  # ======================================================
  # PF COILS
  #
  # Blocks 1-10
  # ======================================================

  [coil_dx]
    type = MFEMPiecewiseConstantFromCSV
    csv_file = pf_currents_pf_only.csv
    header = false
    value_column_index = 2
    attributes = '1 2 3 4 5 6 7 8 9 10'

    # Safe outside PF attributes
    default_value = 1.0
  []

  [coil_dz]
    type = MFEMPiecewiseConstantFromCSV
    csv_file = pf_currents_pf_only.csv
    header = false
    value_column_index = 3
    attributes = '1 2 3 4 5 6 7 8 9 10'

    # Safe outside PF attributes
    default_value = 1.0
  []

  [coil_current]
    type = MFEMPiecewiseConstantFromCSV
    csv_file = pf_currents_pf_only.csv
    header = false
    value_column_index = 4
    attributes = '1 2 3 4 5 6 7 8 9 10'

    default_value = 0.0
  []

  # Bluemira dx and dz are half-dimensions:
  #
  # Area = (2 dx)(2 dz) = 4 dx dz
  #
  # Jtheta = I / Area
  #
  [Jtheta_pf_fun]
    type = MFEMParsedFunction
    expression = 'I/(4*dx*dz)'
    symbol_names = 'I dx dz'
    symbol_values = 'coil_current coil_dx coil_dz'
  []

  # SI source:
  #
  # mu0 * R * Jtheta
  #
  [source_pf_fun]
    type = MFEMParsedFunction
    expression = 'mu0*r*J'
    symbol_names = 'mu0 r J'
    symbol_values =
      '1.25663706212e-6
       cylindrical_r
       Jtheta_pf_fun'
  []


  # ======================================================
  # PLASMA
  #
  # Block 11
  #
  # Spatial Jtheta(R,Z) from Bluemira
  # ======================================================

  [Jtheta_plasma_fun]
    type = PiecewiseBilinear
    data_file = plasma_jtheta_bilinear.csv

    # MFEM mesh coordinates:
    #
    # x = R
    # y = Z
    #
    xaxis = 0
    yaxis = 1
  []

  # SI plasma source:
  #
  # mu0 * R * Jtheta(R,Z)
  #
  [source_plasma_fun]
    type = MFEMParsedFunction
    expression = 'mu0*r*J'
    symbol_names = 'mu0 r J'
    symbol_values =
      '1.25663706212e-6
       cylindrical_r
       Jtheta_plasma_fun'
  []
[]

[FunctorMaterials]
  # ======================================================
  # Axisymmetric A_theta operator
  #
  # -div(R grad Atheta) + Atheta/R
  #
  # RHS = mu0 R Jtheta
  # ======================================================

  [cyl_coeffs]
    type = MFEMGenericFunctorMaterial
    prop_names = 'diffCoef massCoef'
    prop_values = 'cylindrical_r cylindrical_inv_r'
  []
[]

[Kernels]
  # ======================================================
  # A_theta operator
  # ======================================================

  [diffusion]
    type = MFEMDiffusionKernel
    variable = Atheta
    coefficient = diffCoef
  []

  [mass]
    type = MFEMMassKernel
    variable = Atheta
    coefficient = massCoef
  []


  # ======================================================
  # PF source
  #
  # ONLY blocks 1-10
  # ======================================================

  [source_coils]
    type = MFEMDomainLFKernel
    variable = Atheta
    coefficient = source_pf_fun
    block = '1 2 3 4 5 6 7 8 9 10'
  []


  # ======================================================
  # Plasma source
  #
  # ONLY block 11
  # ======================================================

  [source_plasma_current]
    type = MFEMDomainLFKernel
    variable = Atheta
    coefficient = source_plasma_fun
    block = '11'
  []
[]

[AuxKernels]
  # ======================================================
  # PF Jtheta diagnostic
  #
  # No block parameter here:
  # MFEMScalarProjectionAux does not support block.
  #
  # Outside PF blocks:
  # I = 0, dx = 1, dz = 1
  # => Jtheta = 0
  # ======================================================

  [project_pf_Jtheta]
    type = MFEMScalarProjectionAux
    variable = Jtheta_pf
    coefficient = Jtheta_pf_fun
    execute_on = INITIAL
  []


  # ======================================================
  # Plasma Jtheta diagnostic
  #
  # PiecewiseBilinear is projected over the MFEM mesh.
  # The actual plasma source kernel is still restricted
  # to block 11.
  # ======================================================

  [project_plasma_Jtheta]
    type = MFEMScalarProjectionAux
    variable = Jtheta_plasma
    coefficient = Jtheta_plasma_fun
    execute_on = INITIAL
  []


  # ======================================================
  # PF source diagnostic
  # ======================================================

  [project_pf_source]
    type = MFEMScalarProjectionAux
    variable = source_pf
    coefficient = source_pf_fun
    execute_on = INITIAL
  []


  # ======================================================
  # Plasma source diagnostic
  # ======================================================

  [project_plasma_source]
    type = MFEMScalarProjectionAux
    variable = source_plasma
    coefficient = source_plasma_fun
    execute_on = INITIAL
  []


  # ======================================================
  # Magnetic field
  #
  # B_r = -dA_theta/dZ
  #
  # B_z = dA_theta/dR + A_theta/R
  # ======================================================

  [B_from_Atheta]
    type = MFEMMeridionalCurlAux
    variable = B
    source = Atheta
  []
[]

[BCs]
  # Temporary test using boundary attribute 12.
  #
  # Once the dedicated far-field boundary is created,
  # replace 12 with that boundary attribute.
  [essential]
    type = MFEMScalarDirichletBC
    variable = Atheta
    boundary = 12
    coefficient = zero
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []

  [PCG]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-8
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection

    file_base = OutputData/AxisymmetricTokamakPlasmaCoupledAndCoils

    scalar_coefficients =
      'cylindrical_r
       cylindrical_inv_r
       cylindrical_p
       cylindrical_z
       coil_current
       coil_dx
       coil_dz'

    vtk_format = ASCII
  []
[]