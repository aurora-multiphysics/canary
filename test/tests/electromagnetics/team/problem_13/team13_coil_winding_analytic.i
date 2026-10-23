# TEAM Problem 13: the winding current, written down
# ==================================================
#
# Sub-app of team13_coil_source.i, which launches it when coil_winding = analytic,
# the default; team13.i is the input to run. team13_coil_winding_laplace.i is the
# alternative, and either hands back the same thing: applied_current, the winding
# current before it is made divergence free.
#
# The benchmark's coil is wound on a former whose corner arcs are concentric, so
# every turn follows a curve at constant distance d from the square core the arcs
# are centred on. The winding direction is z_hat x grad(d), exactly, and the
# current is that times the ampere-turn density J0. It is divergence free before it
# is discretised; only the faceting of the meshed coil boundary spoils that, which
# is what team13_coil_source.i corrects for.
#
# Nothing is solved for here. The projection is onto the whole mesh: the expression
# is continuous across the coil boundary, so the air elements set the coil's
# boundary edges to the same values the coil elements do, and the current is only
# ever integrated over the coil.

ampere_turns = 1000 # AT; overridden by the main app
fe_order = FIRST # likewise; it has to match the main app's H(curl) space

coil_core = 0.05 # m, half-width of the square the corner arcs are centred on
coil_area = '${fparse 0.025 * 0.100}' # m2
current_density = '${fparse ampere_turns / coil_area}' # A/m2

[Mesh]
  type = MFEMFileMesh
  file = ./team13.e
[]

[Problem]
  type = MFEMProblem
  solve = false
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
[]

[Functions]
  # grad(d) points radially outwards from the core and has unit length, so
  # z_hat x grad(d) runs around the coil at uniform magnitude. The floor on the
  # denominator only guards the core, where the numerators vanish too and the
  # current is zero; outside the coil the expression is not used.
  [applied_current_density]
    type = ParsedVectorFunction
    expression_x = '-J0 * (max(y-a,0)+min(y+a,0)) / max(sqrt((max(x-a,0)+min(x+a,0))^2 + (max(y-a,0)+min(y+a,0))^2), 1e-30)'
    expression_y = 'J0 * (max(x-a,0)+min(x+a,0)) / max(sqrt((max(x-a,0)+min(x+a,0))^2 + (max(y-a,0)+min(y+a,0))^2), 1e-30)'
    expression_z = '0'
    symbol_names = 'J0 a'
    symbol_values = '${current_density} ${coil_core}'
  []
[]

[AuxVariables]
  [applied_current]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  [project_applied_current]
    type = MFEMVectorProjectionAux
    variable = applied_current
    vector_coefficient = applied_current_density
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = MFEMSteady
[]

# No Outputs block: the applied current is transferred to the parent app.
