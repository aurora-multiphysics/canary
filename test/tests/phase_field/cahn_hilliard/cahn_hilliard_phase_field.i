# Cahn-Hilliard spinodal decomposition with the phase field module
# ================================================================
#
# The reference solution cahn_hilliard.i is compared with in
# doc/content/problems/cahn_hilliard.md. It is the phase field module's example
# modules/phase_field/examples/cahn-hilliard/Parsed_SplitCH.i, set up the way
# cahn_hilliard.i has to be:
#
#   - the boundaries are left natural, zero flux, rather than periodic;
#   - the initial noise is the same hash of the coordinates rather than RandomIC,
#     so that both start from identical nodal values;
#   - the time integrator is implicit Euler, which is all MFEMTransient offers,
#     rather than BDF2.
#
# It needs an application with the phase field module, such as phase_field-opt,
# which moose_test is not, and so it is not run by ./run_tests. To reproduce the
# comparison,
#
#   phase_field-opt -i cahn_hilliard_phase_field.i
#   moose_test-opt -i cahn_hilliard.i Outputs/ParaViewDataCollection/vtk_format=ASCII
#   ./plot_cahn_hilliard.py

barrier_height = 0.25
mobility = 1.0
kappa = 0.5

domain_size = 60
elements_per_side = 150

noise_amplitude = 0.1

dt = 2.0
end_time = 20.0

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${elements_per_side}
  ny = ${elements_per_side}
  xmax = ${domain_size}
  ymax = ${domain_size}
[]

[Modules]
  [PhaseField]
    [Conserved]
      [c]
        free_energy = fbulk
        mobility = M
        kappa = kappa_c
        solve_type = REVERSE_SPLIT
      []
    []
  []
[]

[Functions]
  [initial_noise]
    type = ParsedFunction
    expression = 'a * (2 * (43758.5453 * sin(12.9898 * x + 78.233 * y)
                            - floor(43758.5453 * sin(12.9898 * x + 78.233 * y))) - 1)'
    symbol_names = 'a'
    symbol_values = '${noise_amplitude}'
  []
[]

[ICs]
  [c_noise]
    type = FunctionIC
    variable = c
    function = initial_noise
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '${mobility} ${kappa}'
  []
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = fbulk
    coupled_variables = c
    constant_names = W
    constant_expressions = ${barrier_height}
    expression = 'W * (1 - c)^2 * (1 + c)^2'
  []
[]

[AuxVariables]
  [local_energy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [local_energy]
    type = TotalFreeEnergy
    variable = local_energy
    f_name = fbulk
    interfacial_vars = c
    kappa_names = kappa_c
  []
[]

[Preconditioning]
  [full]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = implicit-euler
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  nl_max_its = 20
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-12
  dt = ${dt}
  end_time = ${end_time}
[]

[Postprocessors]
  [average_c]
    type = ElementAverageValue
    variable = c
  []
  [total_free_energy]
    type = ElementIntegralVariablePostprocessor
    variable = local_energy
  []
[]

[VectorPostprocessors]
  # c at every node at the end, to compare with cahn_hilliard.i point by point.
  [c_nodal]
    type = NodalValueSampler
    variable = c
    sort_by = id
    execute_on = FINAL
  []
  [c_line]
    type = LineValueSampler
    variable = c
    start_point = '0 ${fparse domain_size / 2} 0'
    end_point = '${domain_size} ${fparse domain_size / 2} 0'
    num_points = ${fparse elements_per_side + 1}
    sort_by = x
    execute_on = FINAL
  []
  [chemical_potential_line]
    type = LineValueSampler
    variable = chem_pot_c
    start_point = '0 ${fparse domain_size / 2} 0'
    end_point = '${domain_size} ${fparse domain_size / 2} 0'
    num_points = ${fparse elements_per_side + 1}
    sort_by = x
    execute_on = FINAL
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/CahnHilliardPhaseField
    execute_vector_postprocessors_on = FINAL
  []
[]
