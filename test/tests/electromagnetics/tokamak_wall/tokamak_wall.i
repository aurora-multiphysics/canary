# Tokamak wall induction: eddy currents driven by a changing plasma current
# ========================================================================
#
# Run this input; it pulls in plasma_current_source.i as a sub-app. The example,
# its verification and its results are documented in
# doc/content/problems/tokamak_wall_induction.md.
#
# A crude model of a tokamak: a torus of plasma current, inside a toroidal
# stainless steel wall, in free space. Whenever the plasma current changes, the
# change in its magnetic flux induces currents in the wall, which oppose the change
# and then decay resistively. The plasma current density J(x, y, z, t) is any vector
# function of position and time, set in plasma_current_density.i; by default the
# plasma current is ramped up, held, and then lost in a fast current quench, as in
# a disruption.
#
# Formulation
# -----------
# Magnetic vector potential, B = curl(A), in H(curl), in the temporal gauge, so
# that the electric field is E = -dA/dt. Ampere's law with the eddy current
# J_eddy = sigma E and the applied plasma current J gives the weak form
#
#   (sigma dA/dt, A') + (nu curl A, curl A') = (J, A'),
#
# with nu = 1/mu0 everywhere. No time-harmonic or other form is assumed for J; the
# equation is integrated in time with implicit Euler.
#
# The conductivity is that of the wall in the wall. Everywhere else it would be
# zero, which leaves every gradient supported outside the wall in the null space of
# the operator; a small conductivity there regularises it instead. It is far too
# small to carry a current that matters (see nonconductor_conductivity below), and
# because J is made discretely divergence free by the sub-app it does not drive any
# of those gradients.
#
# J is supplied by plasma_current_source.i as the projection of the user's current
# density onto H(curl), applied over the source blocks, together with a correction
# potential, applied everywhere, that removes the current leaking across the faceted
# edge of the source region; see there for why. The sub-app is stepped with this
# app, so that J is evaluated at the time each step is taken to.
#
# The far-field surface carries n x A = 0, which is B.n = 0: no flux leaves the
# air cylinder. It is far enough away that this barely affects the field near the
# wall; see the documentation for the size of the effect.

# 316L stainless steel at room temperature: resistivity 74 micro-ohm cm.
wall_conductivity = '${fparse 1.0 / 7.4e-7}' # S/m

vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m
vacuum_reluctivity = '${fparse 1.0 / vacuum_permeability}' # m/H

# The regularising conductivity outside the wall, eight orders of magnitude below
# the wall's. Its effect on the reported quantities is linear in it: raising it to
# 1 S/m changes the wall current by 5e-5 of its peak, and lowering it to 1e-4 S/m by
# 5e-7, with the same number of conjugate gradient iterations either way.
nonconductor_conductivity = 1e-2 # S/m

# Where the plasma current density is taken from. Widen this to let J reach into
# the vacuum between the plasma and the wall, for a plasma that moves.
source_blocks = 'plasma'

# Order of the H(curl) space the vector potential lives in, and of the H(div) and
# L2 spaces the derived fields are written to. FIRST/CONSTANT and SECOND/FIRST are
# the two pairings that make sense.
fe_order = FIRST
flux_order = CONSTANT

dt = 2.5e-4 # s
end_time = 0.060 # s

[Mesh]
  type = MFEMFileMesh
  file = tokamak_wall.e
[]

[Problem]
  type = MFEMProblem
[]

[FunctorMaterials]
  [Wall]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${wall_conductivity}
    block = wall
  []
  [NonConductors]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${nonconductor_conductivity}
    block = 'plasma vacuum air'
  []
[]

[Functions]
  # With phi_hat = (-y, x, 0) / R the toroidal unit vector, these are -phi_hat, for
  # the toroidal current density sigma E . phi_hat = sigma dA/dt . (-phi_hat), and
  # phi_hat / (2 pi R). Integrating a current density against the second over a
  # volume gives the toroidal current through a poloidal cross-section of it,
  # averaged around the torus: for J = J_phi(R, z) phi_hat, the integral of
  # J_phi / (2 pi R) over the volume is the integral of J_phi over the
  # cross-section. z phi_hat / (2 pi R) likewise gives the vertical moment of that
  # current, the integral of J_phi z over the cross-section, which is zero for
  # currents that are symmetric about the midplane. The floor on R only guards the
  # z axis, where none of them is used.
  [reverse_toroidal_direction]
    type = ParsedVectorFunction
    expression_x = 'y / sqrt(max(x^2 + y^2, 1e-30))'
    expression_y = '-x / sqrt(max(x^2 + y^2, 1e-30))'
    expression_z = '0'
  []
  [toroidal_current_weight]
    type = ParsedVectorFunction
    expression_x = '-y / (2 * pi * max(x^2 + y^2, 1e-30))'
    expression_y = 'x / (2 * pi * max(x^2 + y^2, 1e-30))'
    expression_z = '0'
  []
  [toroidal_current_moment_weight]
    type = ParsedVectorFunction
    expression_x = '-y * z / (2 * pi * max(x^2 + y^2, 1e-30))'
    expression_y = 'x * z / (2 * pi * max(x^2 + y^2, 1e-30))'
    expression_z = '0'
  []
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = ${flux_order}
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = ${flux_order}
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  # The two pieces of the plasma current, from the sub-app. The first is the user's
  # current density, applied over the source blocks; the second is the correction
  # that makes their sum divergence free, and acts over the whole mesh.
  [applied_current_density]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [source_correction_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []

  # E = -dA/dt, B = curl A, and the toroidal component of the eddy current density,
  # J_phi = sigma E . phi_hat, which is the whole of it while everything is
  # axisymmetric.
  [e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [b_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [toroidal_current_density]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # phi_hat / (2 pi R) and z phi_hat / (2 pi R), for the toroidal currents reported
  # below.
  [toroidal_weight]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [toroidal_moment_weight]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[Kernels]
  [eddy_current]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = a_field
    coefficient = conductivity
  []
  [curl_curl]
    type = MFEMCurlCurlKernel
    variable = a_field
    coefficient = ${vacuum_reluctivity}
  []
  # -(J_applied - grad(phi), A'): kernels are all on the left hand side.
  [plasma_current]
    type = MFEMMixedVectorMassKernel
    variable = a_field
    trial_variable = applied_current_density
    coefficient = -1.0
    block = ${source_blocks}
  []
  [plasma_current_correction]
    type = MFEMMixedVectorGradientKernel
    variable = a_field
    trial_variable = source_correction_potential
    coefficient = 1.0
  []
[]

[BCs]
  [tangential_a_far_field]
    type = MFEMVectorTangentialDirichletBC
    variable = a_field
    vector_coefficient = '0 0 0'
    boundary = far_field
  []
[]

[AuxKernels]
  [project_toroidal_weight]
    type = MFEMVectorProjectionAux
    variable = toroidal_weight
    vector_coefficient = toroidal_current_weight
    execute_on = INITIAL
  []
  [project_toroidal_moment_weight]
    type = MFEMVectorProjectionAux
    variable = toroidal_moment_weight
    vector_coefficient = toroidal_current_moment_weight
    execute_on = INITIAL
  []
  [update_e_field]
    type = MFEMSumAux
    variable = e_field
    source_variables = da_field_dt
    scale_factors = -1.0
  []
  [update_b_field]
    type = MFEMCurlAux
    variable = b_field
    source = a_field
  []
  # Built from dA/dt rather than e_field, since MFEM auxkernels are not ordered by
  # the coefficients they read, and this could otherwise run before e_field is
  # updated.
  [update_toroidal_current_density]
    type = MFEMInnerProductAux
    variable = toroidal_current_density
    first_source_vec = da_field_dt
    second_source_vec = reverse_toroidal_direction
    coefficient = conductivity
  []
[]

[MultiApps]
  # Stepped with this app, so that J is projected and corrected at the time each
  # step is taken to.
  [plasma_current_source]
    type = TransientMultiApp
    input_files = plasma_current_source.i
    cli_args = "fe_order=${fe_order};source_blocks='${source_blocks}'"
    execute_on = TIMESTEP_BEGIN
  []
[]

[Transfers]
  [from_plasma_current_source]
    type = MultiAppMFEMCopyTransfer
    source_variables = 'applied_current_density source_correction_potential'
    variables = 'applied_current_density source_correction_potential'
    from_multi_app = plasma_current_source
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    # Not singular: the conductivity is nowhere zero.
    singular = false
    print_level = 0
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-10
    l_max_its = 500
    print_level = 0
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  dt = ${dt}
  start_time = 0.0
  end_time = ${end_time}
[]

[Postprocessors]
  # The toroidal plasma current the solve received, from the applied current density.
  # The correction carries no net toroidal current, being a gradient.
  [plasma_current]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = applied_current_density
    dual_variable = toroidal_weight
    block = ${source_blocks}
  []
  # The toroidal current induced in the wall, sigma E integrated over a poloidal
  # cross-section of it and averaged around the torus.
  [wall_current]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = e_field
    dual_variable = toroidal_weight
    coefficient = conductivity
    block = wall
  []
  # Its vertical moment, the integral of J_phi z over the cross-section, in A m. A
  # plasma that stays on the midplane induces currents symmetric about it, for which
  # this vanishes; one that moves vertically induces currents that oppose the motion,
  # running one way above the midplane and the other way below.
  [wall_current_moment]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = e_field
    dual_variable = toroidal_moment_weight
    coefficient = conductivity
    block = wall
  []
  # Ohmic power dissipated in the wall, the integral of sigma |E|^2.
  [wall_power]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = e_field
    dual_variable = e_field
    coefficient = conductivity
    block = wall
  []

  # B_z at two magnetic probes: on the axis at the centre of the machine, and on the
  # outboard midplane between the plasma and the wall. The outboard probe is turned
  # pi/64 off the y = 0 plane, which puts it inside a toroidal layer of elements
  # rather than on a face between two, where the tangential components of the lowest
  # order B jump, for meshes of 32, 64 or 96 layers.
  [bz_centre]
    type = VectorPostprocessorComponent
    vectorpostprocessor = probes
    vector_name = b_field_2
    index = 0
  []
  [bz_outboard]
    type = VectorPostprocessorComponent
    vectorpostprocessor = probes
    vector_name = b_field_2
    index = 1
  []
[]

[VectorPostprocessors]
  [probes]
    type = MFEMVariablePointValueSampler
    variable = b_field
    points = '0 0 0
              ${fparse 1.35 * cos(pi / 64)} ${fparse 1.35 * sin(pi / 64)} 0'
    outputs = none
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/TokamakWall
  []
  # Every 2 ms is enough to follow the wall current in ParaView; the postprocessor
  # histories above keep the full time resolution.
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TokamakWall
    vtk_format = BINARY
    time_step_interval = 8
  []
[]
