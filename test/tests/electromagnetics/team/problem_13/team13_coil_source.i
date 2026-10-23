# TEAM Problem 13: the coil current
# =================================
#
# Sub-app of team13.i, which is the input to run. It is launched automatically.
#
# The benchmark's coil is a winding carrying a uniform ampere-turn density. A
# sub-app hands back that current, J_applied: team13_coil_winding_analytic.i writes
# it down, and team13_coil_winding_laplace.i finds its direction from the mesh,
# chosen by coil_winding. It is not good enough for the main solve as it stands:
# the curl-curl operator there has every gradient in its null space, so its right
# hand side has to be orthogonal to every discrete gradient. J_applied is not,
# because the meshed coil boundary is faceted, so a little current leaks across it.
#
# This app removes that leak. It solves, over the whole mesh,
#
#   (grad(phi), grad(psi)) = (J_applied, grad(psi))_coil    for all psi,
#
# so that the current the main app actually uses,
#
#   J = J_applied restricted to the coil, minus grad(phi),
#
# satisfies (J, grad(psi)) = 0 for every psi by construction. That is a Helmholtz
# projection: the smallest correction that makes the source discretely divergence
# free, so the current stays as close to the specified uniform density as it can.
# Solving a conduction problem in the coil instead, the way an open coil source
# does, would drive the current to the inside of the bends and stop representing a
# winding.
#
# The two pieces are handed to the main app separately, as applied_current and
# correction_potential, because the first is integrated over the coil alone and the
# second over the whole domain. Sending a single combined field instead would need
# it to carry the coil's sharp edge, which an H(curl) field cannot do.
#
# phi is pinned on far_field. That is not merely to make the solve non-singular: it
# matches the main problem, whose tangential condition on the same boundary leaves
# exactly the gradients of potentials constant there in its null space.

ampere_turns = 1000 # AT; overridden by the main app
fe_order = FIRST # likewise; it has to match the main app's H(curl) space
coil_winding = analytic # likewise; analytic or laplace

[Mesh]
  type = MFEMFileMesh
  file = ./team13.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
[]

[Variables]
  [correction_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [applied_current]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[BCs]
  [pin_far_field]
    type = MFEMScalarDirichletBC
    variable = correction_potential
    boundary = 'far_field'
    coefficient = 0.0
  []
[]

[Kernels]
  [correction_diffusion]
    type = MFEMDiffusionKernel
    variable = correction_potential
  []
  # (-J_applied, grad(psi)) over the coil, which is where the leak comes from.
  [applied_current_divergence]
    type = MFEMMixedVectorWeakDivergenceKernel
    variable = correction_potential
    trial_variable = applied_current
    block = 'Coil'
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [gmres]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 300
  []
[]

[MultiApps]
  # The applied current has to be in place before the solve, so it is fetched at
  # INITIAL.
  [coil_winding]
    type = FullSolveMultiApp
    input_files = 'team13_coil_winding_${coil_winding}.i'
    cli_args = 'ampere_turns=${ampere_turns};fe_order=${fe_order}'
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_coil_winding]
    type = MultiAppMFEMCopyTransfer
    source_variables = 'applied_current'
    variables = 'applied_current'
    from_multi_app = coil_winding
  []
[]

[Executioner]
  type = MFEMSteady
[]

# No Outputs block: both fields are transferred to the main app and written out
# from there.
