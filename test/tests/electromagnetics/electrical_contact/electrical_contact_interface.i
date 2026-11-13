# Electrical contact with an interface kernel
# ===========================================
#
# The problem of electrical_contact.i, a copper rod and an aluminium rod touching over
# a filmed spot on their axis, with the film imposed as a jump in potential across a
# surface of no thickness instead of meshed as a layer. Build the mesh with
#
#   electrical_contact_mesh.py --interface --msh electrical_contact_interface.msh \
#                              --exodus electrical_contact_interface.e
#
# Formulation
# -----------
# The current through the film is J.n = C [phi], with [phi] = phi_copper - phi_aluminium
# the jump in potential across it. The weak form, weighted by 2 pi r for the axisymmetric
# problem as in electrical_contact.i, is
#
#   (2 pi r sigma grad phi, grad phi')_rods + (2 pi r C [phi], [phi'])_spot = 0.
#
# The potential is represented by two variables, one for each rod, both defined over
# the whole mesh. Each is governed by a diffusion kernel restricted to its own rod, and
# is fixed at zero in the other rod everywhere but on the face it shares with its own:
# the constraints below. On the spot, the interface kernel couples them through the
# film. On the gap outside it nothing couples them, so no current crosses it.

# Geometry. These must match the mesh: see electrical_contact_mesh.py.
rod_radius = 5e-3 # m, b
spot_radius = 1e-3 # m, a
rod_half_length = 20e-3 # m, from the contact plane to each electrode

# Materials, at room temperature.
copper_conductivity = 5.8e7 # S/m
aluminium_conductivity = 3.77e7 # S/m
contact_conductance = 1e11 # S/m^2, C, of the film over the spot

applied_voltage = 1e-3 # V

# The resistance of the rods with a uniform current along them.
bulk_resistance = '${fparse (1 / copper_conductivity + 1 / aluminium_conductivity)
                            * rod_half_length / (pi * rod_radius^2)}'

[Mesh]
  type = MFEMFileMesh
  file = ./electrical_contact_interface.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  [copper_weighted_conductivity]
    type = ParsedFunction
    expression = '2 * pi * x * sigma'
    symbol_names = 'sigma'
    symbol_values = '${copper_conductivity}'
  []
  [aluminium_weighted_conductivity]
    type = ParsedFunction
    expression = '2 * pi * x * sigma'
    symbol_names = 'sigma'
    symbol_values = '${aluminium_conductivity}'
  []
  [weighted_contact_conductance]
    type = ParsedFunction
    expression = '2 * pi * x * c'
    symbol_names = 'c'
    symbol_values = '${contact_conductance}'
  []
  # The power dissipated in the film per unit area of the half-plane, 2 pi r C [phi]^2.
  [film_power_density]
    type = MFEMParsedFunction
    expression = 'k * (u - v)^2'
    symbol_names = 'k u v'
    symbol_values = 'weighted_contact_conductance potential_copper potential_aluminium'
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = SECOND
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = SECOND
  []
[]

[Variables]
  [potential_copper]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [potential_aluminium]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [potential_copper_gradient]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [potential_aluminium_gradient]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[Kernels]
  [copper_current_conservation]
    type = MFEMDiffusionKernel
    variable = potential_copper
    coefficient = copper_weighted_conductivity
    block = 'Copper'
  []
  [aluminium_current_conservation]
    type = MFEMDiffusionKernel
    variable = potential_aluminium
    coefficient = aluminium_weighted_conductivity
    block = 'Aluminium'
  []
[]

[InterfaceKernels]
  [film]
    type = MFEMJumpMassInterfaceKernel
    variable = potential_copper
    neighbor_var = potential_aluminium
    boundary = 'spot'
    coefficient = weighted_contact_conductance
  []
[]

[Constraints]
  # Each variable is zero in the other rod, except on the contact plane, which it
  # shares with its own.
  [copper_potential_in_aluminium]
    type = MFEMScalarEssentialConstraint
    variable = potential_copper
    block = 'Aluminium'
    constrain_subdomain_boundary = false
    constrain_exterior_boundary = true
  []
  [aluminium_potential_in_copper]
    type = MFEMScalarEssentialConstraint
    variable = potential_aluminium
    block = 'Copper'
    constrain_subdomain_boundary = false
    constrain_exterior_boundary = true
  []
[]

[BCs]
  [copper_electrode]
    type = MFEMScalarDirichletBC
    variable = potential_copper
    boundary = 'copper_electrode'
    coefficient = ${applied_voltage}
  []
  [aluminium_electrode]
    type = MFEMScalarDirichletBC
    variable = potential_aluminium
    boundary = 'aluminium_electrode'
    coefficient = 0
  []
[]

[AuxKernels]
  [update_copper_gradient]
    type = MFEMGradAux
    variable = potential_copper_gradient
    source = potential_copper
    execute_on = TIMESTEP_END
  []
  [update_aluminium_gradient]
    type = MFEMGradAux
    variable = potential_aluminium_gradient
    source = potential_aluminium
    execute_on = TIMESTEP_END
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 1000
    print_level = 0
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Postprocessors]
  # The power dissipated in each rod and in the film, W.
  [CopperPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = potential_copper_gradient
    dual_variable = potential_copper_gradient
    coefficient = copper_weighted_conductivity
    block = 'Copper'
  []
  [AluminiumPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = potential_aluminium_gradient
    dual_variable = potential_aluminium_gradient
    coefficient = aluminium_weighted_conductivity
    block = 'Aluminium'
  []
  [FilmPower]
    type = MFEMBoundaryIntegralPostprocessor
    coefficient = film_power_density
    boundary = 'spot'
    quadrature_order = 6 # [phi]^2 is of degree 4, and r of degree 1
  []
  [Power]
    type = ParsedPostprocessor
    expression = 'p_cu + p_al + p_film'
    pp_names = 'CopperPower AluminiumPower FilmPower'
    pp_symbols = 'p_cu p_al p_film'
  []
  [Resistance]
    type = ParsedPostprocessor
    expression = 'v^2 / p'
    pp_names = 'Power'
    pp_symbols = 'p'
    constant_names = 'v'
    constant_expressions = '${applied_voltage}'
  []
  [ContactResistance]
    type = ParsedPostprocessor
    expression = 'r - r_bulk'
    pp_names = 'Resistance'
    pp_symbols = 'r'
    constant_names = 'r_bulk'
    constant_expressions = '${bulk_resistance}'
  []
  [FilmPowerFraction]
    type = ParsedPostprocessor
    expression = 'p_film / p'
    pp_names = 'FilmPower Power'
    pp_symbols = 'p_film p'
  []
  [CopperElectrodeCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = potential_copper_gradient
    coefficient = copper_weighted_conductivity
    boundary = 'copper_electrode'
  []
  [AluminiumElectrodeCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = potential_aluminium_gradient
    coefficient = aluminium_weighted_conductivity
    boundary = 'aluminium_electrode'
  []
[]

[VectorPostprocessors]
  # Both potentials across the contact plane; the current density through the spot is
  # C times their difference there.
  [contact_plane_potential]
    type = MFEMVariableLineValueSampler
    variables = 'potential_copper potential_aluminium'
    start_point = '1e-9 0 0'
    end_point = '${fparse spot_radius - 1e-9} 0 0'
    num_points = 201
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/ElectricalContactInterfaceCSV
    execute_on = FINAL
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ElectricalContactInterface
  []
[]
