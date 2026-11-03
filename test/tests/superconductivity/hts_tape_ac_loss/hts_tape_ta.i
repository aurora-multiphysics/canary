# HTS tape AC loss with the T-A formulation
# =========================================
#
# Benchmark #1 of the HTS Modelling Workgroup, solved with the T-A formulation of
# Zhang et al. and Liang et al., which is documented in
# doc/content/problems/hts_tape_ac_loss_ta.md. hts_tape.i solves the same problem on the
# same mesh with the H formulation, and the problem itself is set out in
# doc/content/problems/hts_tape_ac_loss.md. Build the mesh first with
#
#   hts_tape_mesh.py
#
# A single high temperature superconducting tape, 12 mm wide and 1 um thick, carries a
# transport current I(t) = I_m sin(2 pi f t) at 50 Hz in free space, and follows the
# power law E = E_c (|J|/J_c)^n. Only its cross-section is modelled, with the current
# along z and the field in the plane.
#
# Formulation
# -----------
# The current density in the tape is the curl of a current vector potential T, and the
# magnetic flux density everywhere the curl of a magnetic vector potential A z:
#
#   J_z = curl T,   B = curl (A z).
#
# A solves magnetostatics in the whole cross-section with the tape's current as its
# source, and T solves Faraday's law in the tape with E = rho(|J|) J,
#
#   (nu grad A, grad A') - (curl T, A')_tape = 0,
#   (rho(|curl T|) curl T, curl T')_tape + (dA/dt, curl T')_tape = 0,
#
# with nu = 1/mu0. The second is E + dA/dt = -dV/dz tested against curl T': the voltage
# gradient along the tape is uniform across it, and drops out because every T' has no
# tangential component on the tape's boundary. Free space is pure magnetostatics, so it
# needs no resistivity, carries no current, and has nothing to discretise in time.
#
# The tape is thin, so the current density is taken as uniform through its thickness,
# T = T_y(x) y. That is the thin strip approximation of the T-A formulation, made here
# by meshing the tape one lowest order H(curl) element thick and holding the tangential
# T to zero on its faces: T is then left free only on the element edges that cross the
# tape, and each element carries the constant current density (T_y,right - T_y,left)/h.
# hts_tape_mesh.py meshes the tape one element thick by default, and this input relies
# on it: with more, T also has components along the tape on the edges inside it, it is
# no longer of the thin strip form, and the T equation has a null space of gradients
# there that nothing fixes.
#
# Boundary conditions
# -------------------
# Only the quarter of the cross-section with x >= 0 and y >= 0 is meshed. A is even in
# both x and y, so has no normal derivative on either axis, which is the natural
# condition. On the arc of radius R that bounds free space, A = 0: the field there is
# tangential to the arc, as that of a line current I(t) on the axis would be.
#
# The transport current is imposed through T, which is that of a current vector
# potential confined to the tape: the current through the quarter of the tape is the
# circulation of T round its boundary, I(t)/4. The tangential T is zero on the face
# and the end of the tape and on its centre plane y = 0, and is -I(t)/(4b) on its
# centre line x = 0, where b is the half-thickness. The circulation is exact, and
# TapeCurrent reports I(t) to rounding.
#
# The time-dependent value goes on the centre line rather than the end for Newton's
# sake. Each step starts from the previous solution with the new boundary value, so the
# whole of the step's increment in T first appears as current in the single element
# next to that boundary. At the end of the tape, that element is 9.4 um wide and already
# at J_c: the increment adds 5 J_c to it, and Newton then spends twenty iterations
# working its way down the power law, at (1 - 1/n)^n per iteration. In the middle it is
# 187 um wide and well below J_c, and the increment adds a quarter of J_c to it.

critical_current = 300 # A, I_c
current_ratio = 0.5 # I_m / I_c
frequency = 50 # Hz
period = '${fparse 1 / frequency}'

tape_half_width = 6e-3 # m
tape_half_thickness = 0.5e-6 # m, b
tape_quarter_area = '${fparse tape_half_width * tape_half_thickness}' # m^2
critical_current_density = '${fparse critical_current / (4 * tape_quarter_area)}' # A/m^2, J_c
critical_field = 1e-4 # V/m, E_c
n_value = 25

vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m

# Backward Euler, the only scheme MFEMTransient has, is first order in time, so the loss
# converges as 1/steps_per_cycle; see the convergence table in the documentation.
steps_per_cycle = 400
timestep = '${fparse period / steps_per_cycle}'

# The first fortieth of a cycle is taken in steps ramp_refinement times shorter. From
# the virgin state the current grows by half or more in each of the first few steps,
# the front up to which it has penetrated the tape moves further than the elements at
# the edge are wide, and Newton stalls. The loss, which is taken from the second
# quarter of the cycle on, is the same to eleven figures either way.
ramp_refinement = 10
ramp_end = '${fparse period / 40}'

# The loss is taken over the half cycle between the positive and negative peaks of the
# current, by which time the tape is in the cyclic state, and doubled.
loss_start = '${fparse period / 4}'
loss_end = '${fparse 3 * period / 4}'

# Newton's method stops when the residual falls below one, with each equation divided
# through by the tolerance to which it is to be solved. A row of the T equation is the
# difference in E + dA/dt between neighbouring elements of the tape, times the
# thickness b and the step; a row of the A equation is a current imbalance at a node,
# times the step. Unscaled, the A equation's residual is twelve orders of magnitude
# the larger, and is all a relative tolerance sees.
#
# field_tolerance is well clear of the rounding floor of the T equation, which is near
# 1e-5 E_c, and the loss is the same to eleven figures at 1e-2, 1e-3 and 1e-4 E_c.
# current_tolerance is far below the current any step adds, so that no step can pass
# before Newton has iterated. In the shorter steps of the ramp both are effectively
# looser by the factor by which the step is shorter.
field_tolerance = '${fparse 1e-3 * critical_field}' # V/m
current_tolerance = '${fparse 1e-6 * critical_current}' # A
t_equation_scale = '${fparse 1 / (tape_half_thickness * field_tolerance * timestep)}'
a_equation_scale = '${fparse 1 / (current_tolerance * timestep)}'

# T has no meaning outside the tape, but its finite element space spans the whole mesh,
# and its edges in free space need an equation. They are given (T, T') = 0, which with
# the tangential T zero on every boundary free space shares with the tape makes T, and
# the current density, exactly zero there. The edges that carry T in the tape belong to
# no free space element, so this has no effect on the tape. It is scaled as the rows of
# a dimensionless T.
air_t_scale = '${fparse 4 * tape_half_thickness / (critical_current * tape_half_width^2)}' # 1/(A m)

[Mesh]
  type = MFEMFileMesh
  file = ./hts_tape.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  [tape_resistivity]
    type = MFEMParsedFunction
    expression = '(e_c / j_c) * (j / j_c)^(n - 1)'
    symbol_names = 'j e_c j_c n'
    symbol_values = 't_field_curl_mag ${critical_field} ${critical_current_density} ${n_value}'
  []
  # |J| drho/d|J|, for the Jacobian.
  [tape_j_dresistivity_dj]
    type = MFEMParsedFunction
    expression = '(n - 1) * rho'
    symbol_names = 'n rho'
    symbol_values = '${n_value} tape_resistivity'
  []
  # The two, scaled with the rest of the T equation.
  [t_equation_resistivity]
    type = MFEMParsedFunction
    expression = 's * rho'
    symbol_names = 's rho'
    symbol_values = '${t_equation_scale} tape_resistivity'
  []
  [t_equation_j_dresistivity_dj]
    type = MFEMParsedFunction
    expression = 's * j_drho_dj'
    symbol_names = 's j_drho_dj'
    symbol_values = '${t_equation_scale} tape_j_dresistivity_dj'
  []
  # E.J, the power dissipated per unit volume.
  [tape_power_density]
    type = MFEMParsedFunction
    expression = 'rho * j^2'
    symbol_names = 'j rho'
    symbol_values = 't_field_curl_mag tape_resistivity'
  []
  # The tangential T on the centre line of the tape, which carries the transport
  # current. Only its y component is tangential there.
  [centre_line_t]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '-i_m * sin(2 * pi * f * t) / (4 * b)'
    symbol_names = 'i_m f b'
    symbol_values = '${fparse current_ratio * critical_current} ${frequency} ${tape_half_thickness}'
  []
[]

[FESpaces]
  # A is second order. At first order the T equation, which sees A only through its
  # average over each element of the tape, barely sees a current density that
  # alternates from one element to the next, and the current density in the middle of
  # the tape, where the resistivity is negligible and nothing else holds it, comes out
  # ragged by 5% of J_c. At second order it is as smooth as the H formulation's.
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = SECOND
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  # The curl of the lowest order H(curl) space is constant on each element, so these
  # hold J and E.J exactly.
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = H1FESpace
    time_derivative = da_dt
  []
  [t_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  # J_z = curl T, signed.
  [j_field]
    type = MFEMVariable
    fespace = L2FESpace
  []
  [power_density]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[Kernels]
  # The A equation.
  [grad_a]
    type = MFEMDiffusionKernel
    variable = a_field
    coefficient = '${fparse a_equation_scale / vacuum_permeability}'
  []
  [tape_current]
    type = MFEMMixedScalarCurlKernel
    variable = a_field
    trial_variable = t_field
    coefficient = '${fparse -a_equation_scale}'
    block = 'Tape'
  []
  # The T equation.
  [tape_e]
    type = MFEMNLCurlCurlKernel
    variable = t_field
    k_coefficient = t_equation_resistivity
    curlu_dk_dcurlu_coefficient = t_equation_j_dresistivity_dj
    block = 'Tape'
  []
  [tape_da_dt]
    type = MFEMMixedScalarWeakCurlKernel
    variable = t_field
    trial_variable = da_dt
    coefficient = ${t_equation_scale}
    block = 'Tape'
  []
  # T in free space, which is zero.
  [air_t]
    type = MFEMVectorFEMassKernel
    variable = t_field
    coefficient = ${air_t_scale}
    block = 'Air'
  []
[]

[BCs]
  [far_field]
    type = MFEMScalarDirichletBC
    variable = a_field
    coefficient = 0
    boundary = 'far_field'
  []
  # No current crosses the face or the end of the tape, or its centre plane. The part
  # of the x axis in free space is included for T there, which is zero.
  [tape_edges]
    type = MFEMVectorTangentialDirichletBC
    variable = t_field
    vector_coefficient = '0 0'
    boundary = 'x_axis tape_boundary'
  []
  # The single element edge where the tape meets the y axis.
  [transport_current]
    type = MFEMVectorTangentialDirichletBC
    variable = t_field
    vector_coefficient = centre_line_t
    boundary = 'tape_centre_line'
  []
[]

[AuxKernels]
  [update_j_field]
    type = MFEMCurlAux
    variable = j_field
    source = t_field
  []
  # Evaluated everywhere, but only ever averaged over the tape.
  [update_power_density]
    type = MFEMScalarProjectionAux
    variable = power_density
    coefficient = tape_power_density
  []
[]

[Solvers]
  # A direct solver: the system couples a magnetostatic A to a T equation whose
  # resistivity ranges over tens of orders of magnitude, and is small, about 10000
  # unknowns.
  [mumps]
    type = MFEMMUMPS
    print_level = 0
  []
  # See field_tolerance and current_tolerance above for the residual. At n = 25 Newton
  # takes one to eight iterations a step, and the backtracking line search is used in a
  # handful of steps, where the full update overshoots up the power law. At n = 101 a few
  # steps run out of iterations; see the documentation. MFEMTransient carries on
  # regardless, so set print_level = 1 to check a new case.
  [newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 50
    rel_tol = 0
    abs_tol = 1
    print_level = 0
    line_search = backtracking
  []
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = ${loss_end}
  [TimeStepper]
    type = FunctionDT
    time_t = '0 ${ramp_end}'
    time_dt = '${fparse timestep / ramp_refinement} ${timestep}'
    interpolate = false
  []
[]

[Postprocessors]
  # The transport current through the whole tape, from the current density in the
  # quarter that is modelled. It is I(t) to rounding: see the boundary conditions.
  [TapeCurrent]
    type = ParsedPostprocessor
    expression = '4 * area * j_avg'
    pp_names = 'TapeCurrentDensity'
    pp_symbols = 'j_avg'
    constant_names = 'area'
    constant_expressions = '${tape_quarter_area}'
  []
  [TapeCurrentDensity]
    type = MFEMElementAverageValue
    variable = j_field
    block = 'Tape'
    outputs = none
  []
  # Power dissipated per unit length of the whole tape, W/m.
  [Power]
    type = ParsedPostprocessor
    expression = '4 * area * p_avg'
    pp_names = 'TapePowerDensity'
    pp_symbols = 'p_avg'
    constant_names = 'area'
    constant_expressions = '${tape_quarter_area}'
  []
  [TapePowerDensity]
    type = MFEMElementAverageValue
    variable = power_density
    block = 'Tape'
    outputs = none
  []
  # The loss per cycle, J/m, accumulated over the half cycle between the peaks. The
  # integral is the sum of the power at the end of each step, which is what backward
  # Euler dissipates over that step.
  [PowerInLossWindow]
    type = ParsedPostprocessor
    expression = 'if(t > t0 + dt / 2 & t < t1 + dt / 2, p, 0)'
    pp_names = 'Power'
    pp_symbols = 'p'
    constant_names = 't0 t1 dt'
    constant_expressions = '${loss_start} ${loss_end} ${timestep}'
    use_t = true
    outputs = none
  []
  [HalfCycleLoss]
    type = TimeIntegratedPostprocessor
    value = PowerInLossWindow
    time_integration_scheme = implicit-euler
    outputs = none
  []
  [LossPerCycle]
    type = ParsedPostprocessor
    expression = '2 * q'
    pp_names = 'HalfCycleLoss'
    pp_symbols = 'q'
  []
[]

[VectorPostprocessors]
  # J across the half-width of the tape, halfway through its half-thickness, which the
  # benchmark asks for at the peak of the current. The ends are drawn in by 10 nm, off
  # the axis and the edge of the tape, where J is discontinuous.
  [current_density_profile]
    type = MFEMVariableLineValueSampler
    variable = j_field
    start_point = '1e-8 ${fparse tape_half_thickness / 2} 0'
    end_point = '${fparse tape_half_width - 1e-8} ${fparse tape_half_thickness / 2} 0'
    num_points = 241
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/HTSTapeTACSV
    execute_vector_postprocessors_on = NONE
  []
  [PeakProfile]
    type = CSV
    file_base = OutputData/HTSTapeTAPeakProfile
    execute_postprocessors_on = NONE
    sync_times = ${loss_start}
    sync_only = true
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/HTSTapeTA
    time_step_interval = 10
  []
[]
