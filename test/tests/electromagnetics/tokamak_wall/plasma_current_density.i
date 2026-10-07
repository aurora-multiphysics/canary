# Tokamak wall induction: the plasma current density
# ==================================================
#
# This is the file to edit to change the source. It defines
# plasma_current_density, the current density J(x, y, z, t) in A/m^2, as any
# vector function of position and time. plasma_current_source.i includes it and
# evaluates J at every time step, at the time the step is being taken to; nothing
# is assumed about its time dependence, and it need not be separable, harmonic or
# axisymmetric. Only its values in the source blocks are used (source_blocks in
# tokamak_wall.i, the plasma by default), so it does not have to vanish outside
# them.
#
# J should be divergence free, as any real current density is: a current that
# piles up charge cannot be represented by the magnetoquasistatic model. The
# sub-app only removes the small discrete leak through the faceted edge of the
# source region, not a divergence J has in its own right; tokamak_wall.i reports
# the toroidal current the solve actually received as plasma_current, which shows
# when the two differ.
#
# The default is a crude model of a plasma discharge ending in a disruption: the
# plasma current I_p(t) is ramped up, held on a flat top, and then lost in a fast
# current quench. It flows toroidally, peaked on the magnetic axis with a
# parabolic profile over a circular cross-section of minor radius a,
#
#   J = (2 I_p(t) / (pi a^2)) (1 - rho^2 / a^2) phi_hat    for rho < a, 0 beyond,
#
# with phi_hat = (-y, x, 0) / R and rho the distance from the magnetic axis,
# R = R0, z = z_axis(t). This is divergence free, and it vanishes at the plasma
# edge, so it stays continuous wherever the edge falls in the mesh.
#
# By default the plasma stays put, z_axis = 0, and J is just I_p(t) times a fixed
# profile. Setting vertical_displacement makes it a vertical displacement event:
# the plasma drifts upwards by that distance, at a steady rate, between
# displacement_start and displacement_end, and is then quenched where it ends up.
# J is then not separable in space and time. A displaced plasma leaves the plasma
# block, so source_blocks in tokamak_wall.i has to take in the vacuum too:
#
#   tokamak_wall.i source_blocks='plasma vacuum' plasma_current_source:vertical_displacement=0.08
#
# I_p(t) is a piecewise linear function, so it can be replaced by a measured
# waveform, for example with PiecewiseLinear reading a CSV file of time and
# current.

major_radius = 1.0 # m, R0; must match the mesh
plasma_radius = 0.3 # m, a; must match the mesh

peak_plasma_current = 1e6 # A
ramp_end = 0.010 # s, end of the linear ramp up from zero
flat_top_end = 0.030 # s, start of the current quench
quench_end = 0.032 # s, end of the linear current quench, after which I_p = 0

vertical_displacement = 0.0 # m
displacement_start = 0.020 # s
displacement_end = 0.030 # s

[Functions]
  [plasma_current]
    type = PiecewiseLinear
    x = '0 ${ramp_end} ${flat_top_end} ${quench_end}'
    y = '0 ${peak_plasma_current} ${peak_plasma_current} 0'
  []
  [magnetic_axis_height]
    type = PiecewiseLinear
    x = '0 ${displacement_start} ${displacement_end}'
    y = '0 0 ${vertical_displacement}'
  []
  # J_phi / I_p, in 1/m^2: the parabolic profile, normalised to unit total current.
  [current_profile]
    type = ParsedFunction
    expression = '2 / (pi * a^2) * max(1 - ((sqrt(x^2 + y^2) - R0)^2 + (z - za)^2) / a^2, 0)'
    symbol_names = 'R0 a za'
    symbol_values = '${major_radius} ${plasma_radius} magnetic_axis_height'
  []
  # The floor on R only guards the z axis, which is outside the plasma and where
  # the numerators vanish anyway.
  [plasma_current_density]
    type = ParsedVectorFunction
    expression_x = '-y / sqrt(max(x^2 + y^2, 1e-30)) * Ip * j'
    expression_y = 'x / sqrt(max(x^2 + y^2, 1e-30)) * Ip * j'
    expression_z = '0'
    symbol_names = 'Ip j'
    symbol_values = 'plasma_current current_profile'
  []
[]
