# TEAM Problem 7: asymmetrical conductor with a hole
# ==================================================
#
# A racetrack coil carrying a sinusoidal current sits 30 mm above a thick
# aluminium plate with an off-centre square hole, and induces eddy currents in it.
# The flux density above the plate and the eddy current density on its surfaces
# are compared with the measurements of the benchmark at 50 Hz and 200 Hz. The
# benchmark and the results are documented in doc/content/problems/team_problem_7.md.
#
# The mesh is built by team7_mesh.py and is not edited by hand.
#
# Formulation
# -----------
# Time-harmonic magnetic vector potential, with fields varying as exp(i omega t)
# and B = curl(A), in H(curl). In the conductor the electric scalar potential is
# absorbed into A (the A* or modified vector potential formulation), so that
# E = -i omega A there and the eddy current density is J = -i omega sigma A:
#
#   (nu curl A, curl A') + (i omega sigma A, A')_plate = (J_s, A')
#
# nu is the reluctivity of free space everywhere: nothing in the model is
# magnetic. The coil is a stranded winding, so it carries the source current J_s
# and no eddy current of its own; like the air, it has no conductivity.
#
# The hole needs no special treatment. It is simply free space, and the plate
# around it is multiply connected, which the A* formulation handles as it stands.
#
# The null space outside the plate
# --------------------------------
# Inside the plate the conductivity term removes the gradients from the null
# space of the curl-curl operator, but in the air and the coil nothing does, and
# the system is singular. A tree-cotree gauge restricted to those blocks fixes the
# edge degrees of freedom on a spanning tree of them to zero, which removes the
# gradients there without touching the plate. The free space region is multiply
# connected, since it threads the hole, and the cotree takes care of that too.
#
# A gauged system only reproduces the physical solution if its right hand side
# is orthogonal to the gradients it removed, that is, if J_s is discretely
# divergence free. The benchmark's winding current is specified analytically, and
# projecting it onto the mesh is not enough: the meshed coil is faceted where its
# corners are curved, and a little current leaks out across the facets. So J_s is
# built in two parts, as for TEAM Problem 13,
#
#   J_s = J_applied restricted to the coil  -  grad(phi),
#
# where phi, the correction_potential below, solves
#
#   (grad(phi), grad(psi)) = (J_applied, grad(psi))_coil    for all psi,
#
# making (J_s, grad(psi)) = 0 for every psi: the smallest correction that makes
# the source discretely divergence free. Here phi is solved for alongside A, as a
# second complex variable, rather than in a sub-app, since a complex problem
# cannot couple to real variables. The phi equation does not depend on A, so the
# coupled system is block triangular and the direct solver handles it directly.
#
# Boundary conditions
# -------------------
# The outer faces of the free space box are where the benchmark puts them (Fig. 2
# of the problem statement), and carry n x A = 0, which is B.n = 0, as it asks.
# phi is zero there too, which leaves exactly the gradients of potentials that
# vanish on that boundary in the null space the gauge removes.

frequency = 50 # Hz; the benchmark is also run at 200 Hz
angular_frequency = '${fparse 2 * pi * frequency}' # rad/s

plate_conductivity = 3.526e7 # S/m
vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m
vacuum_reluctivity = '${fparse 1.0 / vacuum_permeability}' # m/H

ampere_turns = 2742 # AT
coil_area = '${fparse 0.025 * 0.100}' # m2, the winding's 25 mm x 100 mm section
current_density = '${fparse ampere_turns / coil_area}' # A/m2

# The coil's centre in plan, and the half-width of the square its corner arcs are
# centred on: the outlines are that square offset by 25 mm and 50 mm.
coil_centre_x = 0.194 # m
coil_centre_y = 0.100 # m
coil_core = 0.05 # m

conductor_block = 'Plate'
non_conducting_blocks = 'Air Coil'

[Mesh]
  type = MFEMFileMesh
  file = ./team7.e
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[Variables]
  [a_field]
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [correction_potential]
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # The specified winding current, projected onto H(curl). It is real: the phase of
  # the coil current is the reference, and the benchmark's omega t = 0 is the
  # instant it peaks.
  [applied_current]
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [b_field]
    type = MFEMComplexVariable
    fespace = HDivFESpace
  []
  # The eddy current density -i omega sigma A. It lives in the same space as A and
  # is only meaningful in the plate, which is the only place it is sampled.
  [j_field]
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
[]

[Functions]
  # The specified winding current. d is the distance from the square core the
  # corner arcs are centred on, so grad(d) points outwards across the winding and
  # has unit length, and z_hat x grad(d) runs anticlockwise around the coil seen
  # from above, at uniform magnitude. That is the sense that gives a positive B_z
  # inside the coil, as the measurements have it. The floor on the denominator
  # only guards the core, where the numerators vanish too; outside the coil the
  # expression is not used.
  [applied_current_density]
    type = ParsedVectorFunction
    expression_x = '-J0 * (max(y-cy-a,0)+min(y-cy+a,0)) / max(sqrt((max(x-cx-a,0)+min(x-cx+a,0))^2 + (max(y-cy-a,0)+min(y-cy+a,0))^2), 1e-30)'
    expression_y = 'J0 * (max(x-cx-a,0)+min(x-cx+a,0)) / max(sqrt((max(x-cx-a,0)+min(x-cx+a,0))^2 + (max(y-cy-a,0)+min(y-cy+a,0))^2), 1e-30)'
    expression_z = '0'
    symbol_names = 'J0 cx cy a'
    symbol_values = '${current_density} ${coil_centre_x} ${coil_centre_y} ${coil_core}'
  []
[]

[ICs]
  [applied_current]
    type = MFEMComplexVectorIC
    variable = applied_current
    vector_coefficient_real = applied_current_density
    vector_coefficient_imag = '0 0 0'
  []
[]

[FunctorMaterials]
  [FreeSpace]
    type = MFEMGenericFunctorMaterial
    prop_names = reluctivity
    prop_values = ${vacuum_reluctivity}
  []
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = omega_sigma
    prop_values = '${fparse angular_frequency * plate_conductivity}'
    block = ${conductor_block}
  []
[]

[Kernels]
  [curl_curl]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMCurlCurlKernel
      coefficient = reluctivity
    []
  []
  [eddy_current]
    type = MFEMComplexKernel
    variable = a_field
    block = ${conductor_block}
    [ImagComponent]
      type = MFEMVectorFEMassKernel
      coefficient = omega_sigma
    []
  []
  # The source. Coupled variables that are not solved for are moved to the right
  # hand side, which is what the signs here account for: the two terms together
  # contribute (J_applied - grad(phi), A').
  [coil_current]
    type = MFEMMixedSesquilinearFormKernel
    variable = a_field
    trial_variable = applied_current
    block = 'Coil'
    [RealComponent]
      type = MFEMMixedVectorMassKernel
      coefficient = -1.0
    []
  []
  [coil_current_correction]
    type = MFEMMixedSesquilinearFormKernel
    variable = a_field
    trial_variable = correction_potential
    [RealComponent]
      type = MFEMMixedVectorGradientKernel
      coefficient = 1.0
    []
  []
  # (grad(phi), grad(psi)) - (J_applied, grad(psi))_coil = 0. The weak divergence
  # term is integrated over the coil alone, since that is where the leak comes
  # from; the projected current is not zero in the air elements next to it.
  [correction_diffusion]
    type = MFEMComplexKernel
    variable = correction_potential
    [RealComponent]
      type = MFEMDiffusionKernel
    []
  []
  [applied_current_divergence]
    type = MFEMMixedSesquilinearFormKernel
    variable = correction_potential
    trial_variable = applied_current
    block = 'Coil'
    [RealComponent]
      type = MFEMMixedVectorWeakDivergenceKernel
    []
  []
[]

[BCs]
  [tangential_a_far_field]
    type = MFEMComplexVectorTangentialDirichletBC
    variable = a_field
    boundary = 'far_field'
  []
  [correction_potential_far_field]
    type = MFEMComplexScalarDirichletBC
    variable = correction_potential
    boundary = 'far_field'
  []
[]

[Constraints]
  [tree_cotree_gauge]
    type = MFEMComplexTreeCotreeGaugeEssentialConstraint
    variable = a_field
    block = ${non_conducting_blocks}
    boundary = 'far_field'
  []
[]

[AuxKernels]
  [update_b_field]
    type = MFEMComplexCurlAux
    variable = b_field
    source = a_field
  []
  [update_j_field]
    type = MFEMComplexSumAux
    variable = j_field
    source_variables = 'a_field'
    scale_factors_real = '0.0'
    scale_factors_imag = '${fparse -angular_frequency * plate_conductivity}'
  []
[]

[Solvers]
  [mumps]
    type = MFEMMUMPS
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
  assembly_level = legacy
[]

[VectorPostprocessors]
  # Table 2 of the problem statement: B_z along the lines A1-B1 (y = 72 mm) and
  # A2-B2 (y = 144 mm), both 15 mm above the plate at z = 34 mm, at x = 0 to 288 mm
  # in steps of 18 mm. The z component of b_field is the third column.
  [bz_a1_b1]
    type = MFEMComplexVariablePointValueSampler
    variable = b_field
    points = '0.000 0.072 0.034  0.018 0.072 0.034  0.036 0.072 0.034  0.054 0.072 0.034
              0.072 0.072 0.034  0.090 0.072 0.034  0.108 0.072 0.034  0.126 0.072 0.034
              0.144 0.072 0.034  0.162 0.072 0.034  0.180 0.072 0.034  0.198 0.072 0.034
              0.216 0.072 0.034  0.234 0.072 0.034  0.252 0.072 0.034  0.270 0.072 0.034
              0.288 0.072 0.034'
  []
  [bz_a2_b2]
    type = MFEMComplexVariablePointValueSampler
    variable = b_field
    points = '0.000 0.144 0.034  0.018 0.144 0.034  0.036 0.144 0.034  0.054 0.144 0.034
              0.072 0.144 0.034  0.090 0.144 0.034  0.108 0.144 0.034  0.126 0.144 0.034
              0.144 0.144 0.034  0.162 0.144 0.034  0.180 0.144 0.034  0.198 0.144 0.034
              0.216 0.144 0.034  0.234 0.144 0.034  0.252 0.144 0.034  0.270 0.144 0.034
              0.288 0.144 0.034'
  []
  # Table 3: J_y along y = 72 mm on the top (A3-B3, z = 19 mm) and bottom (A4-B4,
  # z = 0) faces of the plate, at x = 0, 18 and 126 to 288 mm; the line crosses the
  # hole in between. The y component of j_field is the second column. Every one of
  # these points is on the surface of the plate, and three of them on an edge of
  # it or of the hole, where the eddy current density is discontinuous, so each is
  # moved a micrometre into the plate to make sure it is sampled from the plate's
  # side.
  [jy_a3_b3]
    type = MFEMComplexVariablePointValueSampler
    variable = j_field
    points = '0.000001 0.072 0.018999  0.017999 0.072 0.018999  0.126001 0.072 0.018999
              0.144    0.072 0.018999  0.162    0.072 0.018999  0.180    0.072 0.018999
              0.198    0.072 0.018999  0.216    0.072 0.018999  0.234    0.072 0.018999
              0.252    0.072 0.018999  0.270    0.072 0.018999  0.288    0.072 0.018999'
  []
  [jy_a4_b4]
    type = MFEMComplexVariablePointValueSampler
    variable = j_field
    points = '0.000001 0.072 0.000001  0.017999 0.072 0.000001  0.126001 0.072 0.000001
              0.144    0.072 0.000001  0.162    0.072 0.000001  0.180    0.072 0.000001
              0.198    0.072 0.000001  0.216    0.072 0.000001  0.234    0.072 0.000001
              0.252    0.072 0.000001  0.270    0.072 0.000001  0.288    0.072 0.000001'
  []
[]

[Postprocessors]
  # Time-averaged Joule loss in the plate, (1/2) integral of |J|^2 / sigma, which is
  # (omega^2 sigma / 2) integral of |A|^2. Not asked for by the benchmark, but a
  # single number that summarises the whole eddy current distribution.
  [JouleLoss]
    type = MFEMComplexVectorPeriodAveragedPostprocessor
    primal_variable = j_field
    dual_variable = j_field
    coefficient = '${fparse 1.0 / plate_conductivity}'
    block = ${conductor_block}
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/TEAM7CSV
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TEAM7
  []
[]
