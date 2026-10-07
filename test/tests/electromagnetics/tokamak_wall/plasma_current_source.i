# Tokamak wall induction: the plasma current
# ==========================================
#
# Sub-app of tokamak_wall.i, which is the input to run. It is launched
# automatically, and steps in time with it.
#
# The plasma current density J(x, y, z, t) is any vector function, defined in
# plasma_current_density.i; it is taken over the source_blocks only. The vector potential
# equation in the main app has every gradient in its null space outside the wall,
# so its right hand side has to be orthogonal to every discrete gradient there.
# The current density projected onto the mesh is not, even when J itself is
# divergence free: the meshed edge of the source region is faceted, and a little
# current leaks across it.
#
# Every time step, at the time the main app is stepping to, this app projects J
# onto H(curl) and removes that leak. It solves, over the whole mesh,
#
#   (grad(phi), grad(psi)) = (J_applied, grad(psi))_source    for all psi,
#
# so that the current the main app actually uses,
#
#   J = J_applied restricted to the source blocks, minus grad(phi),
#
# satisfies (J, grad(psi)) = 0 for every psi by construction. That is a Helmholtz
# projection: the smallest correction that makes the source discretely divergence
# free. J_applied and phi are handed to the main app separately, because the first
# is integrated over the source blocks alone and the second over the whole mesh.
#
# phi is pinned on far_field, which matches the main problem: its tangential
# condition there leaves exactly the gradients of potentials that are constant on
# that boundary in its null space.
#
# Gradients of phi lie in the H(curl) space of the main app only if the scalar space
# here has the same order, so fe_order is passed down from there.

fe_order = FIRST # overridden by the main app
source_blocks = 'plasma' # likewise

[Mesh]
  type = MFEMFileMesh
  file = tokamak_wall.e
[]

[Problem]
  type = MFEMProblem
[]

# J, as the vector function plasma_current_density, and whatever it is built from.
!include plasma_current_density.i

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
  [source_correction_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [applied_current_density]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  # The applied current is the right hand side of the solve, so it has to be in
  # place before the solve runs rather than at the end of the step. At TIMESTEP_BEGIN
  # the time is already the one the step is being taken to.
  [project_applied_current_density]
    type = MFEMVectorProjectionAux
    variable = applied_current_density
    vector_coefficient = plasma_current_density
    execute_on = TIMESTEP_BEGIN
  []
[]

[BCs]
  [pin_far_field]
    type = MFEMScalarDirichletBC
    variable = source_correction_potential
    boundary = far_field
    coefficient = 0.0
  []
[]

[Kernels]
  [correction_diffusion]
    type = MFEMDiffusionKernel
    variable = source_correction_potential
  []
  # (-J_applied, grad(psi)) over the source, which is where the leak comes from.
  [applied_current_divergence]
    type = MFEMMixedVectorWeakDivergenceKernel
    variable = source_correction_potential
    trial_variable = applied_current_density
    block = ${source_blocks}
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    # phi only has to make the source compatible to well within the tolerance the
    # main app solves to.
    l_tol = 1e-12
    l_max_its = 500
    print_level = 0
  []
[]

[Executioner]
  # Transient only so that it is stepped in time with the main app; there is no time
  # derivative, and each step is an independent solve at that step's time.
  type = MFEMTransient
  device = cpu
[]

# No Outputs block: both fields are transferred to the main app and written out
# from there.
