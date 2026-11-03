# HTS tape AC loss: a superconducting tape carrying an alternating transport current
# ==================================================================================
#
# Benchmark #1 of the HTS Modelling Workgroup. The benchmark and the results are
# documented in doc/content/problems/hts_tape_ac_loss.md. hts_tape_ta.i solves the same
# problem on the same mesh with the T-A formulation. Build the mesh first with
#
#   hts_tape_mesh.py
#
# A single high temperature superconducting tape, 12 mm wide and 1 um thick, carries a
# transport current I(t) = I_m sin(2 pi f t) at 50 Hz in free space. The tape is long
# enough that only its cross-section has to be modelled, with the current along z and
# the field in the plane. Current first flows in at the edges of the tape at the
# critical current density, and only fills the tape as I_m approaches the critical
# current I_c; the field that has penetrated the edges reverses every half cycle and
# dissipates energy as it does so. The benchmark asks for that energy, the AC loss per
# cycle, for I_m/I_c from 0.1 to 0.9, and for the current density across the tape.
#
# The superconductor follows the power law E = E_c (|J|/J_c)^n, the usual description of
# its steep but finite rise in resistivity above J_c. It is compared in the documentation
# with the Norris formula, the exact loss of a thin strip in the critical state, which the
# power law approaches as n grows.
#
# Formulation
# -----------
# Magnetic field H in H(curl), with J = curl H the current density along z. Faraday's law
# with E = rho(|J|) J gives the weak form
#
#   (mu0 dH/dt, H') + (rho(|curl H|) curl H, curl H') = 0,
#
# with the power law resistivity rho = (E_c/J_c) (|J|/J_c)^(n-1) in the tape, solved by
# Newton's method, and a resistivity in free space large enough that it carries no
# current to speak of (see air_resistivity_gradient below). There is no null space to
# remove: the mass term from the time derivative makes the operator definite everywhere.
#
# Boundary conditions
# -------------------
# Only the quarter of the cross-section with x >= 0 and y >= 0 is meshed. The current is
# even in both x and y, so H_x is odd in y and H_y odd in x: the tangential field vanishes
# on both axes, and is set to zero there.
#
# The transport current is imposed through the field on the arc of radius R that bounds
# free space, which is that of a line current I(t) on the axis of the tape:
#
#   H = I(t) / (2 pi r^2) (-y, x).
#
# Its circulation round the boundary of the quarter is I(t)/4, so the tangential
# condition fixes the current through the quarter, whatever the field looks like
# inside. Only the shape of the field on the arc is approximate, by the difference
# between a line current and a strip 12 mm wide seen from R = 60 mm. The circulation
# is not quite exact once discretised: the lowest order H(curl) space takes the
# tangential field at the midpoint of each boundary edge, which overestimates the
# integral along an edge subtending an angle d by a fraction d^2/12. The arc is meshed
# finely enough to put that at 2e-5, and TapeCurrent reports what actually flows.

critical_current = 300 # A, I_c
current_ratio = 0.5 # I_m / I_c
frequency = 50 # Hz
period = '${fparse 1 / frequency}'

tape_half_width = 6e-3 # m
tape_half_thickness = 0.5e-6 # m
tape_quarter_area = '${fparse tape_half_width * tape_half_thickness}' # m^2
critical_current_density = '${fparse critical_current / (4 * tape_quarter_area)}' # A/m^2, J_c
critical_field = 1e-4 # V/m, E_c
n_value = 25

vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m

# Free space has no resistivity of its own, but the curl-curl operator needs one there.
# It has to be large enough that free space carries no current that matters, and no
# larger: the curl-curl term over an element of size h is dt rho / (mu0 h^2) times its
# mass term, and next to the tip of the tape h is the 0.5 um half-thickness. A uniform
# 1 ohm m makes that ratio 1e16 there, the tape's own equations drown in the rounding
# error of free space's, and Newton stalls with the tape unconverged: the loss then
# moves by a percent with the tolerance or the mesh. So the resistivity grows as the
# square of the distance d from the tape, rho = k d^2, as the elements do, and is
# capped far away. The current free space carries goes as the integral of 1/rho, which
# stays small: see TapeCurrent, and the documentation.
air_resistivity_gradient = 1e4 # ohm/m, k
air_resistivity_max = 1.0 # ohm m
# Right against the tape, where d vanishes. Still five orders of magnitude above the
# tape at J_c.
air_resistivity_min = 1e-9 # ohm m

# Backward Euler, the only scheme MFEMTransient has, is first order in time, so the loss
# converges as 1/steps_per_cycle; see the convergence table in the documentation.
steps_per_cycle = 400
timestep = '${fparse period / steps_per_cycle}'

# The loss is taken over the half cycle between the positive and negative peaks of the
# current. The tape has started from the virgin state, but by the first peak it is in
# the cyclic state: the field everywhere in it is what it would be at that point of any
# later cycle. The loss in a half cycle is half the loss per cycle.
loss_start = '${fparse period / 4}'
loss_end = '${fparse 3 * period / 4}'

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
    symbol_values = 'h_field_curl_mag ${critical_field} ${critical_current_density} ${n_value}'
  []
  # |J| drho/d|J|, for the Jacobian.
  [tape_j_dresistivity_dj]
    type = MFEMParsedFunction
    expression = '(n - 1) * rho'
    symbol_names = 'n rho'
    symbol_values = '${n_value} tape_resistivity'
  []
  # E.J, the power dissipated per unit volume.
  [tape_power_density]
    type = MFEMParsedFunction
    expression = 'e_c * j_c * (j / j_c)^(n + 1)'
    symbol_names = 'j e_c j_c n'
    symbol_values = 'h_field_curl_mag ${critical_field} ${critical_current_density} ${n_value}'
  []
  [air_resistivity]
    type = ParsedFunction
    expression = 'min(rho_max, max(rho_min, k * (max(x - a, 0)^2 + max(y - b, 0)^2)))'
    symbol_names = 'k rho_min rho_max a b'
    symbol_values = '${air_resistivity_gradient} ${air_resistivity_min} ${air_resistivity_max}
                     ${tape_half_width} ${tape_half_thickness}'
  []
  # The field of the transport current on the arc bounding free space. It is written
  # with the distance from the axis rather than the radius of the arc: the boundary
  # is the polygon of the arc's chords, and the tangential field is sampled at their
  # midpoints, inside the arc.
  [far_field_h]
    type = ParsedVectorFunction
    expression_x = '-i_m * sin(2 * pi * f * t) / (2 * pi * (x^2 + y^2)) * y'
    expression_y = 'i_m * sin(2 * pi * f * t) / (2 * pi * (x^2 + y^2)) * x'
    symbol_names = 'i_m f'
    symbol_values = '${fparse current_ratio * critical_current} ${frequency}'
  []
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  # The curl of the lowest order H(curl) space is constant on each element, in two
  # dimensions as well as three, so these hold J and E.J exactly.
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
[]

[Variables]
  [h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  # J_z = curl H, signed.
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
  [dh_dt]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = h_field
    coefficient = ${vacuum_permeability}
  []
  [curl_e_tape]
    type = MFEMNLCurlCurlKernel
    variable = h_field
    k_coefficient = tape_resistivity
    curlu_dk_dcurlu_coefficient = tape_j_dresistivity_dj
    block = 'Tape'
  []
  [curl_e_air]
    type = MFEMCurlCurlKernel
    variable = h_field
    coefficient = air_resistivity
    block = 'Air'
  []
[]

[BCs]
  [symmetry]
    type = MFEMVectorTangentialDirichletBC
    variable = h_field
    vector_coefficient = '0 0'
    boundary = 'x_axis y_axis'
  []
  [transport_current]
    type = MFEMVectorTangentialDirichletBC
    variable = h_field
    vector_coefficient = far_field_h
    boundary = 'far_field'
  []
[]

[AuxKernels]
  [update_j_field]
    type = MFEMCurlAux
    variable = j_field
    source = h_field
  []
  # Evaluated everywhere, but only ever averaged over the tape.
  [update_power_density]
    type = MFEMScalarProjectionAux
    variable = power_density
    coefficient = tape_power_density
  []
[]

[Solvers]
  # A direct solver. The linear systems are small, about 4400 unknowns, and the tape,
  # which is all mass term at the scale of its elements, sits beside free space, which is
  # mostly curl-curl. AMS preconditioned PCG does converge on them, and gives the same
  # loss, but takes around 300 iterations a solve and fourteen times as long. With a
  # uniform free space resistivity it does not converge at all.
  [mumps]
    type = MFEMMUMPS
    print_level = 0
  []
  # The residual is dominated by free space, which is linear, so it falls by eight orders
  # of magnitude well before the tape has converged: at rel_tol = 1e-8 the loss is still
  # 0.2% out. From 1e-10 it no longer moves. Here Newton takes two or three iterations a
  # step and never needs to backtrack, but at n = 101, or with fewer steps per cycle, the
  # full step can overshoot so far up the power law that the residual overflows and the
  # run fails; backtracking brings it back. print_level = 1 lists the iterations.
  [newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 50
    rel_tol = 1e-10
    print_level = 0
    line_search = backtracking
  []
[]

[Executioner]
  type = MFEMTransient
  dt = ${timestep}
  start_time = 0.0
  end_time = ${loss_end}
[]

[Postprocessors]
  # The transport current through the whole tape, from the current density in the
  # quarter that is modelled. It should follow I(t): the remainder is carried by free
  # space.
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
    file_base = OutputData/HTSTapeCSV
    execute_vector_postprocessors_on = NONE
  []
  [PeakProfile]
    type = CSV
    file_base = OutputData/HTSTapePeakProfile
    execute_postprocessors_on = NONE
    sync_times = ${loss_start}
    sync_only = true
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/HTSTape
    time_step_interval = 10
  []
[]
