# TEAM Problem 4 (FELIX brick): the applied background field
# ==========================================================
#
# Sub-app of team4_induced_field.i, which is the input you actually run. The
# benchmark is described in doc/content/problems/team_problem_4.md.
#
# The benchmark applies a uniform axial field B_ext(t) = B0 * exp(-t/tau) * z_hat.
# This app builds its spatial part, B0 * z_hat, as an H(curl) field.
#
# That field could be written down directly, but projecting an analytic function
# onto H(curl) leaves a small discrete curl behind, which the main app would read as
# a spurious source current. Solving for a scalar potential and taking its gradient
# instead gives a field that is a discrete gradient by construction, so its discrete
# curl is zero to machine precision. The potential problem is
#
#   div(mu grad(phi)) = 0,    phi = z on the boundaries where B is normal,
#
# whose solution is phi = z. The x = 0 and y = 0 symmetry planes are left with the
# natural condition mu * dphi/dn = 0, which phi = z satisfies exactly, so the
# discrete solution recovers the uniform field to round-off. Finally
#
#   H_ext = (B0 / mu0) grad(phi) = (B0 / mu0) z_hat.

applied_b_field = 0.1 # T, B0 in B_ext(t) = B0 * exp(-t/tau) * z_hat
vacuum_permeability = '${fparse 4 * pi * 1e-7}' # H/m

# B is normal to the z = 0 symmetry plane and, since the brick perturbs it only
# locally, effectively normal to the far-field faces too, so the potential is pinned
# on all four. The x = 0 and y = 0 symmetry planes are left to the natural condition.
normal_b_boundaries = 'z_low z_high x_high y_high'

[Mesh]
  type = MFEMFileMesh
  file = ./team4_symmetrized.e
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  # Both the Dirichlet value and, by construction, the solution everywhere.
  [boundary_magnetic_potential]
    type = ParsedFunction
    expression = z
  []
[]

[FunctorMaterials]
  [Domain]
    type = MFEMGenericFunctorMaterial
    prop_names = permeability
    prop_values = ${vacuum_permeability}
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
[]

[Variables]
  [magnetic_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # Copied to the main app by MultiAppMFEMCopyTransfer, which needs the variable
  # name and finite element space to match on both sides.
  [external_h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = magnetic_potential
    coefficient = permeability
  []
[]

[BCs]
  [applied_potential]
    type = MFEMScalarDirichletBC
    variable = magnetic_potential
    boundary = ${normal_b_boundaries}
    coefficient = boundary_magnetic_potential
  []
[]

[AuxKernels]
  [update_external_h_field]
    type = MFEMGradAux
    variable = external_h_field
    source = magnetic_potential
    scale_factor = '${fparse applied_b_field / vacuum_permeability}'
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [gmres]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-10
    l_max_its = 100
  []
[]

[Executioner]
  type = MFEMSteady
[]

# No Outputs block: external_h_field is transferred to the main app and written out
# from there, so writing it twice would only duplicate files.
