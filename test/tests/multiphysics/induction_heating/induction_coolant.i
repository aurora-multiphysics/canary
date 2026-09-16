# Induction heating: the coolant channel sub-app (thermal hydraulics)
# ===================================================================
#
# A single phase water annulus running the length of the billet, which carries
# away the power the eddy currents deposit in it. This is a sub-app of
# induction_heating.i; see doc/content/problems/induction_heating.md.
#
# The coupling is the one HeatTransferFromExternalAppTemperature1Phase is built
# for. The parent app owns the solid, and sends this app the wall temperature
# T_wall, layer by layer along the channel. This app returns the fluid
# temperature and the heat transfer coefficient, which the parent turns back into
# a convective flux on the billet surface. All of that is set up on the parent
# side by the [CoupledHeatTransfers] block, so this file only has to declare the
# component and the two user objects the transfers read from.
#
#      m_dot, T_in                                          p_out
#         --->  ==============================================  --->
#               ^^^^^^^^^^^^ q_wall from the billet ^^^^^^^^^^
#
# The channel is an annulus wrapped around the billet:
#
#   flow area           A     = pi (0.025^2 - 0.020^2) = 7.069e-4 m^2
#   hydraulic diameter  D_h   = 4 A / (2 pi (0.025 + 0.020)) = 0.01 m
#   heated perimeter    P_hf  = 2 pi 0.020 = 0.1257 m

channel_length = 0.2
n_elems = 10

flow_area = 7.069e-4 # m^2
hydraulic_diameter = 0.01 # m
heated_perimeter = 0.1257 # m

inlet_mass_flow_rate = 0.03 # kg/s
inlet_temperature = 300 # K
outlet_pressure = 5e5 # Pa

# Held fixed rather than taken from a correlation, so that the example does not
# turn on the choice of closure relation.
heat_transfer_coefficient = 2000 # W/m^2/K
friction_factor = 0.02

[GlobalParams]
  initial_p = ${outlet_pressure}
  initial_vel = 0
  initial_T = ${inlet_temperature}
  closures = simple_closures
[]

[FluidProperties]
  # Stiffened gas equation of state fitted to liquid water. At 5 bar the
  # saturation temperature is about 425 K, comfortably above the temperatures
  # this channel reaches, so the flow stays single phase.
  [water]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    q = -1167e3
    q_prime = 0
    p_inf = 1e9
    cv = 1816
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [channel]
    type = FlowChannel1Phase
    position = '0 0 -${fparse channel_length / 2}'
    orientation = '0 0 1'
    length = ${channel_length}
    n_elems = ${n_elems}

    A = ${flow_area}
    D_h = ${hydraulic_diameter}
    f = ${friction_factor}

    fp = water
  []

  [billet_wall]
    type = HeatTransferFromExternalAppTemperature1Phase
    flow_channel = channel
    Hw = ${heat_transfer_coefficient}
    P_hf = ${heated_perimeter}
    initial_T_wall = ${inlet_temperature}
    var_type = elemental
  []

  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'channel:in'
    m_dot = ${inlet_mass_flow_rate}
    T = ${inlet_temperature}
  []
  [outlet]
    type = Outlet1Phase
    input = 'channel:out'
    p = ${outlet_pressure}
  []
[]

[AuxVariables]
  [Hw]
    family = MONOMIAL
    order = CONSTANT
    block = channel
  []
[]

[AuxKernels]
  [Hw_aux]
    type = ADMaterialRealAux
    variable = Hw
    property = 'Hw'
  []
[]

[UserObjects]
  # The parent app's transfers read the fluid state back through these. The layer
  # count must match the 'n_elems' the parent declares for the channel.
  [T_fluid_uo]
    type = LayeredAverage
    direction = z
    variable = T
    num_layers = ${n_elems}
    block = channel
  []
  [Hw_uo]
    type = LayeredAverage
    direction = z
    variable = Hw
    num_layers = ${n_elems}
    block = channel
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2

  dt = 0.1
  dtmin = 1e-6

  solve_type = NEWTON
  line_search = basic
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-4
  nl_max_its = 20

  l_tol = 1e-3
  l_max_its = 300

  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu'
[]

[Postprocessors]
  [T_fluid_avg]
    type = ElementAverageValue
    variable = T
    block = channel
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_wall_avg]
    type = ElementAverageValue
    variable = T_wall
    block = channel
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [htc_avg]
    type = ElementAverageValue
    variable = Hw
    block = channel
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # The first and last elements of the channel, which bracket the rise the
  # coolant takes across the billet.
  [T_fluid_in]
    type = PointValue
    variable = T
    point = '0 0 -0.09'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_fluid_out]
    type = PointValue
    variable = T
    point = '0 0 0.09'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/InductionCoolant
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
