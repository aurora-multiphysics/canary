# Cahn-Hilliard spinodal decomposition with MFEM
# ==============================================
#
# The MFEM counterpart of the phase field module's split Cahn-Hilliard example,
# modules/phase_field/examples/cahn-hilliard/Parsed_SplitCH.i: a conserved order
# parameter c, starting from small noise about c = 0, separates into regions of
# c = -1 and c = +1 which then coarsen. The same free energy, mobility, gradient
# energy coefficient, domain and time step are used. The problem and the comparison
# with the phase field module are documented in
# doc/content/problems/cahn_hilliard.md.
#
# Formulation
# -----------
# The Cahn-Hilliard equation for a free energy
#
#   F = int f(c) + kappa / 2 |grad c|^2 dV,    f(c) = W (c^2 - 1)^2,
#
# is
#
#   dc/dt = div(M grad mu),    mu = f'(c) - kappa lap(c).
#
# The phase field module splits it into two second order equations, one for c and
# one for the whole of mu, and puts f'(c) in the equation for mu. That term is a
# non-linear function of c tested against the mu test functions, and nothing in
# MFEM-MOOSE's base framework can supply it: kernels whose coefficients depend on
# the solution are only re-evaluated during the Newton iterations when they are
# non-linear form integrators, and those may only act on their own variable.
#
# Only the gradient energy part of the chemical potential, w = -kappa lap(c), is
# split off here instead. The bulk part goes into the equation for c, using
# grad f'(c) = f''(c) grad c, as a non-linear diffusion of c:
#
#   (dc/dt, v) + (M f''(c) grad c, grad v) + (M grad w, grad v) = 0
#   (w, q) - (kappa grad c, grad q)                               = 0
#
# The first is MFEMNLDiffusionKernel with k(c) = M f''(c), which is given dk/dc as
# well and so has an exact Jacobian. k is negative where f is concave, which is the
# uphill diffusion that drives the separation. The full chemical potential,
# mu = f'(c) + w, is recovered afterwards as an auxiliary variable.
#
# Every boundary is left natural, which is zero flux of c, M grad(mu).n = 0, and
# zero contact angle, kappa grad(c).n = 0. The phase field example is periodic
# instead; MFEM-MOOSE cannot yet make a periodic mesh.

# The free energy and kinetics of Parsed_SplitCH.i.
barrier_height = 0.25 # W, the constant_expressions 1.0/2^2 there
mobility = 1.0 # M
kappa = 0.5 # gradient energy coefficient

# The domain is 60 by 60 in Parsed_SplitCH.i, with 150 elements along each side.
domain_size = 60
elements_per_side = 150
fe_order = FIRST

# Amplitude of the initial noise about c = 0.
noise_amplitude = 0.1

dt = 2.0
end_time = 20.0

[Mesh]
  type = MFEMMeshGeneratorMesh
  [square]
    type = MFEMGeneratedMeshGenerator
    dim = 2
    nx = ${elements_per_side}
    ny = ${elements_per_side}
    xmax = ${domain_size}
    ymax = ${domain_size}
  []
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  # Noise uniformly distributed in [-noise_amplitude, noise_amplitude] from node to
  # node, like the phase field module's RandomIC, but generated from the coordinates
  # by a hash, fract(43758.5453 sin(12.9898 x + 78.233 y)), so that it is the same in
  # serial and in parallel and can be given to the phase field module too.
  [initial_noise]
    type = ParsedFunction
    expression = 'a * (2 * (43758.5453 * sin(12.9898 * x + 78.233 * y)
                            - floor(43758.5453 * sin(12.9898 * x + 78.233 * y))) - 1)'
    symbol_names = 'a'
    symbol_values = '${noise_amplitude}'
  []

  # M f''(c) and its derivative, M f'''(c), for the non-linear diffusion of c.
  [bulk_diffusivity]
    type = MFEMParsedFunction
    expression = 'M * 4 * W * (3 * c^2 - 1)'
    symbol_names = 'M W c'
    symbol_values = '${mobility} ${barrier_height} c'
  []
  [bulk_diffusivity_derivative]
    type = MFEMParsedFunction
    expression = 'M * 24 * W * c'
    symbol_names = 'M W c'
    symbol_values = '${mobility} ${barrier_height} c'
  []

  # The full chemical potential, f'(c) + w, which is the w of the phase field module.
  [chemical_potential_function]
    type = MFEMParsedFunction
    expression = '4 * W * c * (c^2 - 1) + w'
    symbol_names = 'W c w'
    symbol_values = '${barrier_height} c w'
  []
  # The two parts of the free energy density.
  [bulk_energy_density]
    type = MFEMParsedFunction
    expression = 'W * (c^2 - 1)^2'
    symbol_names = 'W c'
    symbol_values = '${barrier_height} c'
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
  []
  # The gradient of c lies in the H(curl) space of the same order, so it is held
  # exactly, and the energy integrals below are taken over it.
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
[]

[Variables]
  [c]
    type = MFEMVariable
    fespace = H1FESpace
  []
  # The gradient energy part of the chemical potential, -kappa lap(c).
  [w]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [chemical_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [c_gradient]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  # A unit vector, which turns the vector inner product postprocessor into a plain
  # integral of its coefficient.
  [x_direction]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[ICs]
  [c_noise]
    type = MFEMScalarIC
    variable = c
    coefficient = initial_noise
  []
[]

[Kernels]
  # (dc/dt, v) + (M f''(c) grad c, grad v) + (M grad w, grad v) = 0
  [c_dot]
    type = MFEMTimeDerivativeMassKernel
    variable = c
  []
  [bulk_flux]
    type = MFEMNLDiffusionKernel
    variable = c
    k_coefficient = bulk_diffusivity
    dk_du_coefficient = bulk_diffusivity_derivative
  []
  [gradient_energy_flux]
    type = MFEMMixedGradGradKernel
    variable = c
    trial_variable = w
    coefficient = ${mobility}
  []

  # (w, q) - (kappa grad c, grad q) = 0
  [w]
    type = MFEMMassKernel
    variable = w
  []
  [laplacian_c]
    type = MFEMMixedGradGradKernel
    variable = w
    trial_variable = c
    coefficient = '${fparse -kappa}'
  []
[]

[AuxKernels]
  [update_chemical_potential]
    type = MFEMScalarProjectionAux
    variable = chemical_potential
    coefficient = chemical_potential_function
  []
  [update_c_gradient]
    type = MFEMGradAux
    variable = c_gradient
    source = c
  []
  [project_x_direction]
    type = MFEMVectorProjectionAux
    variable = x_direction
    vector_coefficient = '1 0'
  []
[]

[Solvers]
  # The Jacobian couples c and w and is not symmetric, so it is factorised directly.
  [lu]
    type = MFEMMUMPS
    print_level = 0
  []
  [newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 20
    abs_tol = 1e-12
    rel_tol = 1e-9
    print_level = 1
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  assembly_level = legacy
  dt = ${dt}
  start_time = 0.0
  end_time = ${end_time}
[]

[Postprocessors]
  # The total free energy, which the flow can only decrease, in its two parts.
  [bulk_energy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = x_direction
    dual_variable = x_direction
    coefficient = bulk_energy_density
    outputs = none
  []
  [gradient_energy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = c_gradient
    dual_variable = c_gradient
    coefficient = '${fparse kappa / 2}'
    outputs = none
  []
  [total_free_energy]
    type = ParsedPostprocessor
    expression = 'bulk_energy + gradient_energy'
    pp_names = 'bulk_energy gradient_energy'
  []
  # The mean of c, which the flow conserves.
  [average_c]
    type = MFEMElementAverageValue
    variable = c
  []
[]

[VectorPostprocessors]
  # c and the chemical potential along the horizontal centre line at the end.
  [c_line]
    type = MFEMVariableLineValueSampler
    variable = c
    start_point = '0 ${fparse domain_size / 2} 0'
    end_point = '${domain_size} ${fparse domain_size / 2} 0'
    num_points = ${fparse elements_per_side + 1}
    execute_on = FINAL
  []
  [chemical_potential_line]
    type = MFEMVariableLineValueSampler
    variable = chemical_potential
    start_point = '0 ${fparse domain_size / 2} 0'
    end_point = '${domain_size} ${fparse domain_size / 2} 0'
    num_points = ${fparse elements_per_side + 1}
    execute_on = FINAL
  []
[]

[Outputs]
  # Everything here is computed at the end of each step. At the start, the initial
  # condition on c is not ordered before the objects that read c through a parsed
  # coefficient, so they would see c = 0; there is no row for t = 0.
  [CSV]
    type = CSV
    file_base = OutputData/CahnHilliard
    execute_on = 'TIMESTEP_END FINAL'
    execute_vector_postprocessors_on = FINAL
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/CahnHilliard
    vtk_format = BINARY
  []
[]
