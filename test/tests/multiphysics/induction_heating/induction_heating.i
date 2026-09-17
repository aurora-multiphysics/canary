# Induction heating of a water-cooled billet: the parent app
# ==========================================================
#
# A MultiApp problem in which the three physics are each solved by the code best
# suited to them, and exchanged between apps every timestep:
#
#   induction_heating.i   this file. Heat conduction in the billet, on a libMesh
#                         mesh with the heat transfer module.
#   induction_em.i        time-harmonic eddy currents in the billet and the air
#                         around it, solved with MFEM in H(curl).
#   induction_coolant.i   the coolant annulus, solved with the thermal
#                         hydraulics module as a 1-D single phase flow channel.
#
# See doc/content/problems/induction_heating.md for the physics, the geometry and
# the results. Run mesh.i first if the two meshes are not present.
#
#             +-----------------------------------------+
#             |          induction_heating.i            |
#             |      billet temperature, libMesh        |
#             +-----------------------------------------+
#                 |   ^                     |   ^
#      temperature|   |joule_heating T_wall |   |T_fluid, htc
#                 v   |                     v   |
#      +----------------------+   +--------------------------+
#      |    induction_em.i    |   |    induction_coolant.i   |
#      |    eddy currents,    |   |    coolant channel,      |
#      |         MFEM         |   |   thermal hydraulics     |
#      +----------------------+   +--------------------------+
#
# The loop is explicit: within a timestep each app sees the other apps' state
# from the end of the previous step. Both couplings are stiff enough in only one
# direction for that to be stable at the timestep used here, but a problem run
# closer to thermal runaway would want Picard iteration (fixed_point_max_its on
# the executioner) instead.

billet_half_length = 0.1
channel_n_elems = 10

# Stainless steel.
thermal_conductivity = 16 # W/m/K
specific_heat = 500 # J/kg/K
density = 7900 # kg/m^3

initial_temperature = 300 # K

[Mesh]
  file = induction_billet.e
[]

[Variables]
  [temperature]
    family = LAGRANGE
    order = FIRST
    initial_condition = ${initial_temperature}
  []
[]

[AuxVariables]
  # Cycle-averaged Ohmic dissipation, from the electromagnetic app. It is
  # elemental because MFEM computes it in a constant L2 space, which is the only
  # discontinuous space the MFEM to libMesh transfer accepts.
  [joule_heating]
    family = MONOMIAL
    order = CONSTANT
  []
  # Coolant state, from the thermal hydraulics app via [CoupledHeatTransfers].
  [T_fluid]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${initial_temperature}
  []
  [htc]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0
  []
[]

[Kernels]
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = temperature
  []
  [conduction]
    type = HeatConduction
    variable = temperature
  []
  # The eddy current heating enters as a volumetric source. CoupledForce adds
  # -(joule_heating, psi), which is the source on the left hand side.
  [ohmic_heating]
    type = CoupledForce
    variable = temperature
    v = joule_heating
  []
[]

[Materials]
  [billet]
    type = HeatConductionMaterial
    thermal_conductivity = ${thermal_conductivity}
    specific_heat = ${specific_heat}
  []
  [billet_density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = ${density}
  []
[]

# The curved side of the billet is cooled by the channel; its two flat ends are
# left adiabatic, which is the natural condition.
#
# This block does all of the solid side work for the coupling: it adds a
# CoupledConvectiveHeatFluxBC on 'billet_side' driven by T_fluid and htc, a
# NearestPointLayeredSideAverage that reduces the surface temperature to one
# value per channel element, and the three transfers that move T_wall out and
# T_fluid and htc back. The axis description must match the FlowChannel1Phase in
# induction_coolant.i.
[CoupledHeatTransfers]
  [coolant]
    boundary = billet_side
    T = temperature
    T_wall = T_wall
    T_fluid = T_fluid
    htc = htc
    multi_app = coolant
    T_fluid_user_objects = T_fluid_uo
    htc_user_objects = Hw_uo

    position = '0 0 -${billet_half_length}'
    orientation = '0 0 1'
    length = ${fparse 2 * billet_half_length}
    n_elems = ${channel_n_elems}
  []
[]

[MultiApps]
  [electromagnetics]
    type = TransientMultiApp
    input_files = induction_em.i
    execute_on = TIMESTEP_END
    # induction_em.i names its own output files so that it can be run standalone
    # and gold checked on its own; redirect them here so that running it as a
    # sub-app does not overwrite the results of that run.
    cli_args = 'Outputs/Values/file_base=OutputData/InductionEM_sub;'
               'Outputs/ParaViewDataCollection/file_base=OutputData/InductionEM_sub'
  []
  [coolant]
    type = TransientMultiApp
    input_files = induction_coolant.i
    execute_on = TIMESTEP_END
    # The flow channel needs a much smaller timestep than the thermal diffusion
    # in the solid does, so it takes several of its own per parent step.
    sub_cycling = true
  []
[]

[Transfers]
  # Temperature out to the electromagnetic app, which uses it to set the
  # billet conductivity, and the resulting dissipation back.
  [temperature_to_em]
    type = MultiApplibMeshToMFEMShapeEvaluationTransfer
    to_multi_app = electromagnetics
    source_variables = temperature
    variables = temperature
  []
  [joule_heating_from_em]
    type = MultiAppMFEMTolibMeshShapeEvaluationTransfer
    from_multi_app = electromagnetics
    source_variables = joule_heating
    variables = joule_heating
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  # The billet surface and the coolant settle in about twenty seconds, but the
  # interior is governed by conduction across the radius, r^2 rho c_p / k, which
  # is nearer a hundred. Three hundred seconds reaches a steady state in every
  # quantity. Halving dt changes the final temperatures by 2e-5 relative, so the
  # two second step is well inside the accurate range.
  dt = 2
  end_time = 300

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6

  # Additive Schwarz with an incomplete LU on each block. This is a linear,
  # symmetric, 1305 unknown conduction problem, so the preconditioner barely
  # matters for convergence; what it must not do is reach for a GPU. On a PETSc
  # built --with-cuda, PCCreate_HYPRE initialises a CUDA device, and every MPI
  # rank then takes a context of several hundred megabytes. BoomerAMG would buy
  # nothing on a problem this size and would make the run die with
  # 'cudaErrorMemoryAllocation' above about ten ranks on an 8 GB card.
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm ilu'
[]

[Postprocessors]
  [total_power]
    type = ElementIntegralVariablePostprocessor
    variable = joule_heating
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_temperature]
    type = NodalExtremeValue
    variable = temperature
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [surface_temperature]
    type = SideAverageValue
    variable = temperature
    boundary = billet_side
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [coolant_temperature]
    type = ElementAverageValue
    variable = T_fluid
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [Values]
    type = CSV
    file_base = OutputData/InductionHeating
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Exodus]
    type = Exodus
    file_base = OutputData/InductionHeating
  []
[]
