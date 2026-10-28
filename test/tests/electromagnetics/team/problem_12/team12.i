# TEAM Problem 12 (FELIX cantilevered beam): eddy currents coupled to deflection
# ===============================================================================
#
# Run this input; it pulls in team12_beam.i as a sub-app. The benchmark and the
# comparison with experiment are documented in doc/content/problems/team_problem_12.md.
# Build the mesh first with
#
#   moose_test-opt -i team12_mesh.i --mesh-only team12.e
#
# An aluminium cantilever lies along x, with its 3.175 mm thickness along y. It sits in a
# steady field B_x along its length and a transient field B_y across its thickness,
#
#   B_ext(t) = (B_x, B0 * exp(-t/tau), 0).
#
# The decay of B_y induces eddy currents circulating in the plane of the beam. Their
# component J_z across the width crosses B_x, giving a force density J_z * B_x along y
# that bends the beam. As the beam moves through B_x it induces a motional EMF, v x B,
# which drives currents opposing the motion, so the deflection and the current are
# coupled: the larger B_x, the more strongly the motion is damped.
#
# This app solves the eddy current problem on the whole domain; team12_beam.i solves for
# the beam's elastic motion. Every time step this app hands the sub-app the Lorentz force
# density and receives the beam velocity back, which it uses for the motional EMF of the
# next step.
#
# Formulation
# -----------
# The field is split as H = H_ext + H_ind, with H_ext = B_ext / mu0 uniform, so the whole
# current density is J = curl(H_ind). In the beam Ohm's law for a moving conductor reads
# J = sigma * (E + v x B), so E = rho * J - v x B, and Faraday's law gives the weak form
# solved here for H_ind in H(curl),
#
#   (mu dH_ind/dt, H') + (rho curl(H_ind), curl(H')) = -(dB_ext/dt, H') + (v x B, curl(H')).
#
# The motional EMF v x B is discontinuous across the surface of the beam, so it is applied
# in this weak form, through (v x B, curl(H')), rather than projected onto H(curl).
#
# The self-field of the eddy currents is a few mT at most against B_x >= 0.2 T, so both
# the Lorentz force and the motional EMF are evaluated with the applied field alone.
#
# Boundary conditions
# -------------------
# The induced field has decayed away at the outer boundary, so its tangential component
# is set to zero there.

beam_conductivity = 25.32e6 # S/m, aluminium, Table 5 of the LSTC test case document
vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m

# The curl-curl operator is singular where the resistivity vanishes, so free space is
# given a large but finite resistivity instead of an infinite one. At more than seven
# orders of magnitude above the beam it leaves a negligible current outside it.
vacuum_resistivity = 1.0 # ohm m

steady_field = 0.5 # T, B_x; the benchmark takes 0.2, 0.5, 0.7 and 0.9 T
transient_field = 0.055 # T, B0 in B_y = B0 * exp(-t/tau)
decay_time = 0.0066 # s, tau

[Mesh]
  type = MFEMFileMesh
  file = ./team12.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  [applied_b_field]
    type = ParsedVectorFunction
    expression_x = 'bx'
    expression_y = 'b0 * exp(-t / tau)'
    expression_z = '0'
    symbol_names = 'bx b0 tau'
    symbol_values = '${steady_field} ${transient_field} ${decay_time}'
  []
  # 1/dx over the two element layers either side of mid-length, on the z > 0 half of the
  # beam. Integrating J_x against it gives the current through half of the cross-section,
  # averaged over the 19.4 mm slab; with div(J) = 0 that is exactly the mean of the
  # section current there, which is flat about mid-length.
  [current_window]
    type = ParsedFunction
    expression = 'if(x > x0 - left & x < x0 + right & z > 0, 1 / (left + right), 0)'
    symbol_names = 'x0 left right'
    symbol_values = '0.2435 ${fparse 0.1675 / 17} ${fparse 0.0765 / 8}'
  []
  [unit_x]
    type = ParsedVectorFunction
    expression_x = '1'
    expression_y = '0'
    expression_z = '0'
  []
  # -dB_ext/dt; B_x is constant, so only the decaying component contributes.
  [applied_field_source]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '(b0 / tau) * exp(-t / tau)'
    expression_z = '0'
    symbol_names = 'b0 tau'
    symbol_values = '${transient_field} ${decay_time}'
  []
[]

[FunctorMaterials]
  # conductor_indicator restricts the Lorentz force and the motional EMF to the beam.
  [Beam]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity permeability conductor_indicator'
    prop_values = '${fparse 1.0 / beam_conductivity} ${vacuum_permeability} 1.0'
    block = Beam
  []
  [Vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity permeability conductor_indicator'
    prop_values = '${vacuum_resistivity} ${vacuum_permeability} 0.0'
    block = Vacuum
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
  # The spaces of the fields exchanged with team12_beam.i, which must match the ones
  # they are defined on there.
  [VelocityFESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = SECOND
    range_dim = 3
  []
  [ForceFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = FIRST
    range_dim = 3
  []
[]

[Variables]
  [induced_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  # J = curl(H_ind), which lives naturally in H(div).
  [j_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  # The unit vector along the beam, which RT0 represents exactly.
  [x_direction]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  # Received from team12_beam.i.
  [beam_velocity]
    type = MFEMVariable
    fespace = VelocityFESpace
  []
  [motional_emf]
    type = MFEMVariable
    fespace = ForceFESpace
  []
  # Sent to team12_beam.i.
  [lorentz_force_density]
    type = MFEMVariable
    fespace = ForceFESpace
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
    type = MFEMVectorFEDomainLFKernel
    variable = induced_h_field
    vector_coefficient = applied_field_source
  []
  [motional_emf]
    type = MFEMVectorFEDomainLFCurlKernel
    variable = induced_h_field
    vector_coefficient = motional_emf
  []
[]

[BCs]
  [tangential_induced_h_bc]
    type = MFEMVectorTangentialDirichletBC
    variable = induced_h_field
    vector_coefficient = '0 0 0'
    boundary = 'left right bottom top back front'
  []
[]

[AuxKernels]
  [x_direction]
    type = MFEMVectorProjectionAux
    variable = x_direction
    vector_coefficient = unit_x
    execute_on = INITIAL
  []
  # v x B from the velocity at the end of the previous step, ahead of this step's solve.
  [update_motional_emf]
    type = MFEMCrossProductAux
    variable = motional_emf
    first_source_vec = beam_velocity
    second_source_vec = applied_b_field
    coefficient = conductor_indicator
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [update_j_field]
    type = MFEMCurlAux
    variable = j_field
    source = induced_h_field
  []
  [update_lorentz_force_density]
    type = MFEMCrossProductAux
    variable = lorentz_force_density
    first_source_vec = induced_h_field_curl
    second_source_vec = applied_b_field
    coefficient = conductor_indicator
  []
[]

[MultiApps]
  [beam]
    type = TransientMultiApp
    input_files = team12_beam.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [to_beam]
    type = MultiAppMFEMCopyTransfer
    source_variables = lorentz_force_density
    variables = lorentz_force_density
    to_multi_app = beam
  []
  [from_beam]
    type = MultiAppMFEMCopyTransfer
    source_variables = beam_velocity
    variables = beam_velocity
    from_multi_app = beam
  []
  # A single sub-app, so the reduction over sub-apps just passes its value through.
  [tip_deflection]
    type = MultiAppPostprocessorTransfer
    from_postprocessor = TipDeflection
    to_postprocessor = TipDeflection
    from_multi_app = beam
    reduction_type = average
  []
  [station_deflection]
    type = MultiAppPostprocessorTransfer
    from_postprocessor = StationDeflection
    to_postprocessor = StationDeflection
    from_multi_app = beam
    reduction_type = average
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
  # Implicit Euler is the only scheme available, so the step has to be small against
  # both tau and the 65 ms period of the beam's first bending mode; see the time step
  # study in the documentation.
  dt = 1e-4
  start_time = 0.0
  end_time = 0.15
[]

[Postprocessors]
  # Current crossing half of the beam's cross-section at mid-length: the current
  # circulating around the beam.
  [CirculatingCurrent]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = j_field
    dual_variable = x_direction
    coefficient = current_window
    block = Beam
  []
  # Displacements received from team12_beam.i, at the free end and 167 mm from it.
  [TipDeflection]
    type = Receiver
  []
  [StationDeflection]
    type = Receiver
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/TEAM12CSV
  []
[]
