# TEAM Problem 4 (FELIX brick): eddy currents in a decaying background field
# ==========================================================================
#
# Run this input; it pulls in team4_external_source_field.i as a sub-app. See
# README.md in this directory for the benchmark description and reference data.
#
# A hollow aluminium brick sits in a uniform axial field that is switched off with
# an exponential decay,
#
#   B_ext(t) = B0 * exp(-t/tau) * z_hat,
#
# and the decay drives eddy currents that circulate around the hole. Only one octant
# is modelled, using the symmetry planes at x = 0, y = 0 and z = 0.
#
# Formulation
# -----------
# The total field is split as H = H_ext + H_ind. H_ext comes from the sub-app and is
# a discrete gradient, so curl(H_ext) = 0 and the whole current density is
#
#   J = curl(H_ind).
#
# Substituting the split into Faraday's law with E = rho * J gives the weak form
# solved here for H_ind in H(curl):
#
#   (mu dH_ind/dt, H') + (rho curl(H_ind), curl(H')) = -(dB_ext/dt, H').
#
# Boundary conditions
# -------------------
# Under reflection in z the induced field has H_x and H_y odd, so both vanish on the
# z = 0 plane; those are exactly the tangential components there, so the homogeneous
# tangential condition is the correct symmetry condition rather than an approximation.
# On the far-field faces it holds because the induced field has decayed away.
#
# The x = 0 and y = 0 planes are deliberately left free. Their tangential components
# include H_z, which is not zero, so constraining them would be wrong. The natural
# condition of this weak form, n x E = 0, is the physically correct one: it stops
# current crossing the symmetry plane, which is what J_x = 0 on x = 0 and J_y = 0 on
# y = 0 already require.

conductor_domains = 'Brick'
vacuum_domain = 'Vacuum'

vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m
conductor_conductivity = 2.538e7 # S/m, aluminium alloy
conductor_resistivity = '${fparse 1.0 / conductor_conductivity}' # ohm m

# The curl-curl operator is singular where the resistivity vanishes, so free space is
# given a large but finite resistivity instead of an infinite one. At eight orders of
# magnitude above the brick it leaves a negligible current outside the conductor
# while keeping the system solvable.
vacuum_resistivity = 1.0 # ohm m

# B0 enters through H_ext, which the sub-app builds; only tau is needed here.
decay_time = 0.0119 # s, tau

tangential_induced_h_boundaries = 'z_low z_high x_high y_high'

[Mesh]
  type = MFEMFileMesh
  file = ./team4_symmetrized.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  # Scales H_ext, which has magnitude B0/mu0, to give dB_ext/dt = -(B0/tau) exp(-t/tau).
  [db_ext_dt_over_h_ext]
    type = ParsedFunction
    expression = '-(mu0 / tau) * exp(-t / tau)'
    symbol_names = 'mu0 tau'
    symbol_values = '${vacuum_permeability} ${decay_time}'
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
  # Received from the sub-app; the finite element space must match the one it is
  # defined on there.
  [external_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  # J = curl(H_ind), which lives naturally in H(div).
  [j_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[Kernels]
  [dh_ind_dt]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = induced_h_field
    coefficient = permeability
  []
  [curl_curl_h_ind]
    type = MFEMCurlCurlKernel
    variable = induced_h_field
    coefficient = resistivity
  []
  [db_ext_dt]
    type = MFEMMixedVectorMassKernel
    variable = induced_h_field
    trial_variable = external_h_field
    coefficient = db_ext_dt_over_h_ext
  []
[]

[BCs]
  [tangential_induced_h_bc]
    type = MFEMVectorTangentialDirichletBC
    variable = induced_h_field
    vector_coefficient = '0 0 0'
    boundary = ${tangential_induced_h_boundaries}
  []
[]

[AuxKernels]
  [update_j_field]
    type = MFEMCurlAux
    variable = j_field
    source = induced_h_field
  []
[]

[MultiApps]
  # H_ext has no time dependence of its own, so it is solved for once and reused.
  [external_source_field]
    type = FullSolveMultiApp
    input_files = team4_external_source_field.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_external_field]
    type = MultiAppMFEMCopyTransfer
    source_variables = external_h_field
    variables = external_h_field
    from_multi_app = external_source_field
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-12
    l_max_its = 500
  []
[]

[Executioner]
  type = MFEMTransient
  dt = 0.001
  start_time = 0.0
  end_time = 0.02
[]

[Postprocessors]
  # Ohmic power dissipated in the brick, P = integral of rho * |J|^2.
  [OhmicHeating]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = resistivity
    primal_variable = j_field
    dual_variable = j_field
    block = ${conductor_domains}
  []
  # Current crossing the two limbs of the brick wall that the modelled octant cuts
  # through, on the y = 0 and x = 0 symmetry planes respectively.
  [InnerLimbCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = j_field
    boundary = semi_xsection_in
  []
  [OuterLimbCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    variable = j_field
    boundary = semi_xsection_out
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/TEAM4CSV
  []
  # Each snapshot of the field is around 18 MB, so only every fifth step is written.
  # The postprocessor histories above keep the full time resolution.
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TEAM4
    time_step_interval = 5
  []
[]
