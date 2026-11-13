# Electrical contact: current through a filmed contact spot between two rods
# ==========================================================================
#
# A copper rod and an aluminium rod, 10 mm across, are pressed together end to end.
# Real surfaces are rough, and the rods touch only where the asperities on their
# faces meet: here, over a single circular spot 2 mm across on their axis, the
# "a-spot" of Holm's theory of electric contacts. The spot is covered by a thin
# film of oxide or tarnish, and elsewhere the faces are separated by a gap that
# carries no current. A potential difference between electrodes 20 mm either side
# of the contact drives a current through it. The example and its results are
# documented in doc/content/problems/electrical_contact.md. Build the mesh first with
#
#   electrical_contact_mesh.py
#
# The contact adds two resistances in series to those of the rods. The current
# has to crowd through the spot, which takes a constriction resistance; Holm's
# value for an equipotential spot of radius a between half-spaces is
# (rho_1 + rho_2) / (4 a), reduced here by the finite radius of the rods. And it
# has to cross the film, whose resistance is 1 / (C pi a^2) for a contact
# conductance C, in S/m^2, if the current crosses it uniformly. Holm's formula is
# their sum, but the two interact: the film evens out the current over the spot,
# which raises the constriction resistance above that of an equipotential spot.
# The resistance is compared in the documentation with a semi-analytic solution,
# contact_reference.py, which also shows that Holm's sum is a lower bound.
#
# Formulation
# -----------
# Electric potential phi in H1, with J = -sigma grad phi and div J = 0. The
# problem is axisymmetric and meshed in the (r, z) half-plane, with r along x and z
# along y. MFEM has no axisymmetric coordinate system of its own, so the weak form
# is written out in full,
#
#   (2 pi r sigma grad phi, grad phi') = 0,
#
# with the 2 pi r in the conductivity. The weight makes every integral over the
# half-plane an integral over the volume: the power below is that dissipated in
# the rods, not per radian of them.
#
# The film has no thickness to speak of, and the potential jumps across it by
# J.n / C. MFEM-MOOSE has no interface kernel to impose the jump, so the film is
# meshed as a layer of elements of thickness t, a micron, over the spot, with
# conductivity sigma_f = C t. Its resistance through the thickness, t / sigma_f per
# unit area, is then 1 / C whatever t, while the current it can carry along its
# length, sigma_f t = C t^2, vanishes as t does. The gap outside the spot, of the
# same thickness, is not meshed at all.
#
# Boundary conditions
# -------------------
# The electrodes are held at the applied potential, on the copper, and at zero, on
# the aluminium. They are far enough from the contact, four rod radii, that the
# current through them is uniform to a part in a million. Everywhere else, the
# axis, the sides of the rods and their faces either side of the gap, carries no
# current: the natural boundary condition.

# Geometry. These must match the mesh: see electrical_contact_mesh.py.
rod_radius = 5e-3 # m, b
spot_radius = 1e-3 # m, a
rod_half_length = 20e-3 # m, from the contact plane to each electrode
film_thickness = 1e-6 # m, t

# Materials, at room temperature.
copper_conductivity = 5.8e7 # S/m, annealed copper, 100% IACS
aluminium_conductivity = 3.77e7 # S/m, pure aluminium
contact_conductance = 1e11 # S/m^2, C, of the film over the spot
film_conductivity = '${fparse contact_conductance * film_thickness}' # S/m

applied_voltage = 1e-3 # V

# The resistance the rods would have with a uniform current along them, from the
# electrodes to the faces of the film, which the contact adds to.
bulk_resistance = '${fparse (1 / copper_conductivity + 1 / aluminium_conductivity)
                            * (rod_half_length - film_thickness / 2) / (pi * rod_radius^2)}'

[Mesh]
  type = MFEMFileMesh
  file = ./electrical_contact.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  # The conductivity of each block, weighted by 2 pi r for the axisymmetric weak form.
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
  [film_weighted_conductivity]
    type = ParsedFunction
    expression = '2 * pi * x * sigma'
    symbol_names = 'sigma'
    symbol_values = '${film_conductivity}'
  []
  # |J| = sigma |grad phi|.
  [current_density_magnitude]
    type = MFEMParsedFunction
    expression = 'sigma * g'
    symbol_names = 'sigma g'
    symbol_values = 'conductivity potential_grad_mag'
  []
[]

[FunctorMaterials]
  [copper]
    type = MFEMGenericFunctorMaterial
    prop_names = 'conductivity weighted_conductivity'
    prop_values = '${copper_conductivity} copper_weighted_conductivity'
    block = 'Copper'
  []
  [aluminium]
    type = MFEMGenericFunctorMaterial
    prop_names = 'conductivity weighted_conductivity'
    prop_values = '${aluminium_conductivity} aluminium_weighted_conductivity'
    block = 'Aluminium'
  []
  [film]
    type = MFEMGenericFunctorMaterial
    prop_names = 'conductivity weighted_conductivity'
    prop_values = '${film_conductivity} film_weighted_conductivity'
    block = 'Film'
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = SECOND
  []
  # The gradient of the H1 space lies in the H(curl) space of the same order, so
  # grad phi is held exactly, and so is the power computed from it.
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = SECOND
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = FIRST
  []
[]

[Variables]
  [potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [potential_gradient]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [current_density]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[Kernels]
  [current_conservation]
    type = MFEMDiffusionKernel
    variable = potential
    coefficient = weighted_conductivity
  []
[]

[BCs]
  [copper_electrode]
    type = MFEMScalarDirichletBC
    variable = potential
    boundary = 'copper_electrode'
    coefficient = ${applied_voltage}
  []
  [aluminium_electrode]
    type = MFEMScalarDirichletBC
    variable = potential
    boundary = 'aluminium_electrode'
    coefficient = 0
  []
[]

[AuxKernels]
  [update_potential_gradient]
    type = MFEMGradAux
    variable = potential_gradient
    source = potential
    execute_on = TIMESTEP_END
  []
  [update_current_density]
    type = MFEMScalarProjectionAux
    variable = current_density
    coefficient = current_density_magnitude
    execute_on = TIMESTEP_END
  []
[]

[Solvers]
  # The system is symmetric positive definite. AMG preconditioned CG takes 23
  # iterations with the default film, and 150 for a clean contact, when the layer is
  # far more conductive than the rods. The tolerance is tight because the resistance
  # is wanted to a part in a million.
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
  # The power dissipated in the rods and the film, W.
  [Power]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = potential_gradient
    dual_variable = potential_gradient
    coefficient = weighted_conductivity
  []
  [FilmPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = potential_gradient
    dual_variable = potential_gradient
    coefficient = weighted_conductivity
    block = 'Film'
  []
  # The resistance between the electrodes, from the power. Of all the ways the
  # discrete potential could carry the current, the solution dissipates least, but
  # it cannot carry it in quite the way the true one would; so it dissipates more,
  # and the resistance converges from below as the mesh is refined.
  [Resistance]
    type = ParsedPostprocessor
    expression = 'v^2 / p'
    pp_names = 'Power'
    pp_symbols = 'p'
    constant_names = 'v'
    constant_expressions = '${applied_voltage}'
  []
  # What the contact adds to the resistance of the rods: the constriction and the
  # film together, the quantity Holm's formula is for.
  [ContactResistance]
    type = ParsedPostprocessor
    expression = 'r - r_bulk'
    pp_names = 'Resistance'
    pp_symbols = 'r'
    constant_names = 'r_bulk'
    constant_expressions = '${bulk_resistance}'
  []
  # The share of the power dissipated in the film, which would be its resistance
  # 1 / (C pi a^2) over the whole resistance, if the current crossed it uniformly.
  [FilmPowerFraction]
    type = ParsedPostprocessor
    expression = 'p_film / p'
    pp_names = 'FilmPower Power'
    pp_symbols = 'p_film p'
  []
  # The current through each electrode, A, from the flux of the potential gradient
  # across it. Both should be V / Resistance: they are a check on the solve, but
  # less accurate than the power, as they depend on the gradient only at the
  # boundary.
  [CopperElectrodeCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = potential_gradient
    coefficient = copper_weighted_conductivity
    boundary = 'copper_electrode'
  []
  [AluminiumElectrodeCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = potential_gradient
    coefficient = aluminium_weighted_conductivity
    boundary = 'aluminium_electrode'
  []
[]

[VectorPostprocessors]
  # |J| through the spot, across the middle of the film. The ends are drawn in by
  # a nanometre, off the axis and the edge of the spot.
  [spot_current_density]
    type = MFEMVariableLineValueSampler
    variables = 'current_density'
    start_point = '1e-9 0 0'
    end_point = '${fparse spot_radius - 1e-9} 0 0'
    num_points = 201
  []
  # The potential along the axis, through both rods and the film.
  [axial_potential]
    type = MFEMVariableLineValueSampler
    variables = 'potential'
    start_point = '0 ${fparse -rod_half_length} 0'
    end_point = '0 ${rod_half_length} 0'
    num_points = 401
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/ElectricalContactCSV
    execute_on = FINAL
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ElectricalContact
  []
[]
