# TEAM Problem 13: the winding current, found from the mesh
# =========================================================
#
# Sub-app of team13_coil_source.i, which launches it when coil_winding = laplace;
# team13.i is the input to run. It is the alternative to
# team13_coil_winding_analytic.i, and hands back the same thing: applied_current, the
# winding current before it is made divergence free.
#
# The benchmark's coil is a winding carrying a uniform ampere-turn density J0 around
# a closed loop. This app finds the direction of that winding from the mesh rather
# than having it written down, which is what a coil with no closed-form path needs:
# the potential V that satisfies Laplace's equation in
# the coil and jumps by a fixed amount across a cut through it drives a current
# around the loop, so grad(V) runs along the winding. Only its direction is kept,
#
#   J_applied = J0 grad(V) / |grad(V)|,
#
# because the conduction current itself would crowd to the inside of the bends and
# stop representing a winding.
#
# For this coil it is an approximation. In the straight sections and in the corner
# arcs it is the direction of the turns, but near the ends of the arcs, where the
# uniform field of a straight section meets the 1/r field of an arc, it picks up a
# component across them, and once normalised it is no longer divergence free before
# it is discretised. The reported flux densities move by at most 0.5%.
#
# The jump is imposed as for any closed conductor (see MOOSE's
# test/tests/mfem/submeshes/cut_closed_coil.i). A potential V_ext that is 1 on the cut
# and 0 everywhere else is kept only on a transition region, the elements touching
# the cut from one side, and V = V_ext + V_ind with V_ind solving
#
#   (grad(V_ind), grad(v))_coil = -(grad(V_ext), grad(v))_transition    for all v.
#
# The whole calculation is done on a submesh holding only the coil, and the result
# transferred to the full mesh afterwards. Projecting onto the full mesh directly
# would let the air elements set the coil's boundary edges as well, and they would
# set them to zero.
#
# Which side of the cut the transition region falls on is up to the orientation the
# mesh gives the cut, and that sets the sign of grad(V). The sign is fixed here
# instead: the winding runs anticlockwise seen from +z, as in the problem statement,
# so that the reported flux densities come out positive.

ampere_turns = 1000 # AT; overridden by the main app
fe_order = FIRST # likewise; it has to match the main app's H(curl) space

coil_area = '${fparse 0.025 * 0.100}' # m2
current_density = '${fparse ampere_turns / coil_area}' # A/m2

[Mesh]
  type = MFEMFileMesh
  file = ./team13.e
[]

[Problem]
  type = MFEMProblem
[]

[SubMeshes]
  [cut]
    type = MFEMCutTransitionSubMesh
    cut_boundary = 'coil_cut'
    block = 'Coil'
    transition_subdomain = 'coil_transition'
    transition_subdomain_boundary = 'coil_transition_boundary'
    closed_subdomain = 'coil_closed'
  []
  [coil]
    type = MFEMDomainSubMesh
    block = 'coil_closed'
  []
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
  [CoilH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
    submesh = coil
  []
  [CoilHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
    submesh = coil
  []
  [TransitionH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
    submesh = cut
  []
  [TransitionHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
    submesh = cut
  []
[]

[Functions]
  # Anticlockwise about z. Only its sign against grad(V) is used, and that is the
  # same everywhere in the coil, which surrounds the z axis.
  [anticlockwise]
    type = ParsedVectorFunction
    expression_x = '-y'
    expression_y = 'x'
    expression_z = '0'
  []
  # The factor turning -grad(V) into the applied current. The floor on |grad(V)|
  # only guards against division by zero; the field does not vanish in the coil.
  [winding_scale]
    type = MFEMParsedFunction
    expression = 'if(orientation > 0, ${current_density}, -${current_density}) / max(E, 1e-30)'
    symbol_names = 'orientation E'
    symbol_values = 'winding_orientation coil_e_field_mag'
  []
[]

[Variables]
  [coil_induced_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
[]

[AuxVariables]
  [coil_external_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
  [transition_external_potential]
    type = MFEMVariable
    fespace = TransitionH1FESpace
  []
  [transition_external_e_field]
    type = MFEMVariable
    fespace = TransitionHCurlFESpace
  []
  [coil_external_e_field]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
  [coil_induced_e_field]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
  # -grad(V), the field that drives the current around the loop.
  [coil_e_field]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
  [coil_anticlockwise]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
  [coil_applied_current]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
  # The result, on the full mesh, for the parent app.
  [applied_current]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[ICs]
  [coil_external_potential_ic]
    type = MFEMScalarBoundaryIC
    variable = coil_external_potential
    boundary = 'coil_cut'
    coefficient = 1.0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = coil_induced_potential
  []
  [cut_source]
    type = MFEMMixedGradGradKernel
    variable = coil_induced_potential
    trial_variable = coil_external_potential
    block = 'coil_transition'
  []
[]

[AuxKernels]
  [update_coil_induced_e_field]
    type = MFEMGradAux
    variable = coil_induced_e_field
    source = coil_induced_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_transition_external_e_field]
    type = MFEMGradAux
    variable = transition_external_e_field
    source = transition_external_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_coil_e_field]
    type = MFEMSumAux
    variable = coil_e_field
    source_variables = 'coil_induced_e_field coil_external_e_field'
    execute_on = TIMESTEP_END
  []
  [project_coil_anticlockwise]
    type = MFEMVectorProjectionAux
    variable = coil_anticlockwise
    vector_coefficient = anticlockwise
    execute_on = TIMESTEP_END
  []
  # Reads coil_e_field and winding_orientation through coefficients, which do not
  # order it after the objects that compute them, so it is held back to FINAL.
  [project_coil_applied_current]
    type = MFEMVectorProjectionAux
    variable = coil_applied_current
    vector_coefficient = coil_e_field
    coefficient = winding_scale
    execute_on = FINAL
  []
[]

[Postprocessors]
  # Positive when -grad(V) runs anticlockwise.
  [winding_orientation]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = coil_e_field
    dual_variable = coil_anticlockwise
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [to_transition]
    type = MFEMSubMeshTransfer
    from_variable = coil_external_potential
    to_variable = transition_external_potential
    execute_on = TIMESTEP_END
  []
  [from_transition]
    type = MFEMSubMeshTransfer
    from_variable = transition_external_e_field
    to_variable = coil_external_e_field
    execute_on = TIMESTEP_END
  []
  [from_coil]
    type = MFEMSubMeshTransfer
    from_variable = coil_applied_current
    to_variable = applied_current
    execute_on = FINAL
  []
[]

[Solvers]
  # V_ind is only fixed up to a constant, the coil having no Dirichlet boundary, but
  # the right hand side is compatible and PCG copes with the null space.
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 300
    print_level = 0
  []
[]

[Executioner]
  type = MFEMSteady
[]

# No Outputs block: the applied current is transferred to the parent app.
