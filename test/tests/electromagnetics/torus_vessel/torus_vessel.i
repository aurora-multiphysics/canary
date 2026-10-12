# Eddy currents and electromagnetic loads in a toroidal vessel shell driven by an
# axisymmetric source field written by gridded_fields.py. Here the source is
# the uniform vertical field ramp of generate_ramp_source.py; a FREEGSnke export
# from vde/freegsnke_source.py is read the same way.
#
# The shell is the TorusSheath block: major radius 5 m, minor radii 0.5 m to
# 1 m, symmetric about the z axis. TorusCore and Exterior are non-conducting.
#
# Reduced vector potential formulation. The source field enters only through
# its electric field E_s = -(1/R) dpsi_s/dt phi_hat = dpsi_s/dt (y, -x, 0) / R^2,
# and MFEM solves for the vector potential A of the eddy currents alone:
#
#   sigma dA/dt + curl(nu curl A) = sigma E_s   in the shell,
#             curl(nu curl A) = 0               elsewhere,
#
# with A x n = 0 on the outer boundary. The eddy current density is then
# J = sigma (E_s - dA/dt), the total field B = B_s + curl A, and the load is the
# Lorentz force density f = J x B.
#
# Everything the loads need is evaluated on a submesh of the shell, so the
# source data only has to cover the conductors.

# 316L stainless steel at room temperature [S/m].
conductivity = 1.35e6
# 1 / mu0 [m/H].
reluctivity = 795774.715459
# Mass term making the curl-curl operator definite in the non-conducting
# blocks, where it otherwise has the gradients in its kernel. It is 1e-6 of
# the curl-curl term for fields varying over 1 m, so it does not perturb A.
regularisation = 0.795774715459

[Mesh]
  type = MFEMFileMesh
  file = embedded_concentric_torus.e
[]

[Problem]
  type = MFEMProblem
[]

[SubMeshes]
  [vessel]
    type = MFEMDomainSubMesh
    block = TorusSheath
  []
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
  [VesselHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = vessel
  []
  [VesselHDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
    submesh = vessel
  []
  [VesselL2VectorFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = CONSTANT
    submesh = vessel
  []
[]

[Variables]
  # Vector potential of the eddy currents. Its time derivative is da_dt.
  [a]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [b_eddy]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [vessel_da_dt]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
  [vessel_b_eddy]
    type = MFEMVariable
    fespace = VesselHDivFESpace
  []
  [e_source]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
  # E_s - dA/dt, the total electric field in the shell.
  [e_total]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
  [b_source]
    type = MFEMVariable
    fespace = VesselHDivFESpace
  []
  [b_total]
    type = MFEMVariable
    fespace = VesselHDivFESpace
  []
  [force_density]
    type = MFEMVariable
    fespace = VesselL2VectorFESpace
  []
  # Weights turning volume integrals over the shell into toroidally averaged
  # currents and into force components.
  [phi_hat_over_2pi_r]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
  [r_hat]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
  [z_hat]
    type = MFEMVariable
    fespace = VesselHCurlFESpace
  []
[]

[UserObjects]
  [source]
    type = SolutionUserObject
    mesh = ramp_source_out.e
    system_variables = 'dpsi_dt b_r b_z b_phi'
  []
[]

[Functions]
  # Source fields mapped from the (R, Z) half-plane onto the 3D mesh, whose
  # symmetry axis is z.
  [dpsi_dt]
    type = Axisymmetric2D3DSolutionFunction
    solution = source
    from_variables = dpsi_dt
    3d_axis_point2 = '0 0 1'
  []
  [b_phi]
    type = Axisymmetric2D3DSolutionFunction
    solution = source
    from_variables = b_phi
    3d_axis_point2 = '0 0 1'
  []
  [b_poloidal_x]
    type = Axisymmetric2D3DSolutionFunction
    solution = source
    from_variables = 'b_r b_z'
    component = 0
    3d_axis_point2 = '0 0 1'
  []
  [b_poloidal_y]
    type = Axisymmetric2D3DSolutionFunction
    solution = source
    from_variables = 'b_r b_z'
    component = 1
    3d_axis_point2 = '0 0 1'
  []
  [b_poloidal_z]
    type = Axisymmetric2D3DSolutionFunction
    solution = source
    from_variables = 'b_r b_z'
    component = 2
    3d_axis_point2 = '0 0 1'
  []
  [e_source_function]
    type = ParsedVectorFunction
    expression_x = 'y * dpsi_dt / (x^2 + y^2)'
    expression_y = '-x * dpsi_dt / (x^2 + y^2)'
    expression_z = '0'
    symbol_names = dpsi_dt
    symbol_values = dpsi_dt
  []
  [sigma_e_source]
    type = ParsedVectorFunction
    expression_x = 'sigma * y * dpsi_dt / (x^2 + y^2)'
    expression_y = '-sigma * x * dpsi_dt / (x^2 + y^2)'
    expression_z = '0'
    symbol_names = 'dpsi_dt sigma'
    symbol_values = 'dpsi_dt ${conductivity}'
  []
  [b_source_function]
    type = ParsedVectorFunction
    expression_x = 'bx - y * bphi / sqrt(x^2 + y^2)'
    expression_y = 'by + x * bphi / sqrt(x^2 + y^2)'
    expression_z = 'bz'
    symbol_names = 'bx by bz bphi'
    symbol_values = 'b_poloidal_x b_poloidal_y b_poloidal_z b_phi'
  []
  [phi_hat_over_2pi_r_function]
    type = ParsedVectorFunction
    expression_x = '-y / (2 * pi * (x^2 + y^2))'
    expression_y = 'x / (2 * pi * (x^2 + y^2))'
    expression_z = '0'
  []
  [r_hat_function]
    type = ParsedVectorFunction
    expression_x = 'x / sqrt(x^2 + y^2)'
    expression_y = 'y / sqrt(x^2 + y^2)'
    expression_z = '0'
  []
  [z_hat_function]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = '1'
  []
[]

[Kernels]
  [eddy_current]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = a
    coefficient = ${conductivity}
    block = TorusSheath
  []
  [magnetic_energy]
    type = MFEMCurlCurlKernel
    variable = a
    coefficient = ${reluctivity}
  []
  [regularisation]
    type = MFEMVectorFEMassKernel
    variable = a
    coefficient = ${regularisation}
    block = 'TorusCore Exterior'
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = a
    vector_coefficient = sigma_e_source
    block = TorusSheath
  []
[]

[BCs]
  [far_field]
    type = MFEMVectorTangentialDirichletBC
    variable = a
    vector_coefficient = '0 0 0'
    boundary = Exterior
  []
[]

[AuxKernels]
  [b_eddy]
    type = MFEMCurlAux
    variable = b_eddy
    source = a
    execute_on = TIMESTEP_END
  []
  [e_source]
    type = MFEMVectorProjectionAux
    variable = e_source
    vector_coefficient = e_source_function
    execute_on = TIMESTEP_END
  []
  [e_total]
    type = MFEMSumAux
    variable = e_total
    source_variables = 'e_source vessel_da_dt'
    scale_factors = '1 -1'
    execute_on = TIMESTEP_END
  []
  [b_source]
    type = MFEMVectorProjectionAux
    variable = b_source
    vector_coefficient = b_source_function
    execute_on = TIMESTEP_END
  []
  [b_total]
    type = MFEMSumAux
    variable = b_total
    source_variables = 'b_source vessel_b_eddy'
    execute_on = TIMESTEP_END
  []
  [force_density]
    type = MFEMCrossProductAux
    variable = force_density
    first_source_vec = e_total
    second_source_vec = b_total
    coefficient = ${conductivity}
    execute_on = TIMESTEP_END
  []
  [phi_hat_over_2pi_r]
    type = MFEMVectorProjectionAux
    variable = phi_hat_over_2pi_r
    vector_coefficient = phi_hat_over_2pi_r_function
    execute_on = INITIAL
  []
  [r_hat]
    type = MFEMVectorProjectionAux
    variable = r_hat
    vector_coefficient = r_hat_function
    execute_on = INITIAL
  []
  [z_hat]
    type = MFEMVectorProjectionAux
    variable = z_hat
    vector_coefficient = z_hat_function
    execute_on = INITIAL
  []
[]

[Transfers]
  [da_dt_to_vessel]
    type = MFEMSubMeshTransfer
    from_variable = da_dt
    to_variable = vessel_da_dt
    execute_on = TIMESTEP_END
  []
  [b_eddy_to_vessel]
    type = MFEMSubMeshTransfer
    from_variable = b_eddy
    to_variable = vessel_b_eddy
    execute_on = TIMESTEP_END
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-12
    l_max_its = 500
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  # The shell's L/R time is about 1 s, so ten steps leave about 1e-3 of the
  # initial transient in the eddy current.
  dt = 1
  end_time = 10
[]

[Postprocessors]
  # Toroidally averaged current in the shell, the integral of J . phi_hat / (2 pi R)
  # over its volume. For an axisymmetric current this is the current through a
  # poloidal cross-section.
  [toroidal_current]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = e_total
    dual_variable = phi_hat_over_2pi_r
    coefficient = ${conductivity}
    execute_on = TIMESTEP_END
  []
  # The same for sigma E_s, the current the shell approaches once the eddy
  # currents stop changing under a steady ramp.
  [saturated_toroidal_current]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = e_source
    dual_variable = phi_hat_over_2pi_r
    coefficient = ${conductivity}
    execute_on = TIMESTEP_END
  []
  # L2 norm of dA/dt in the shell. Under a steady ramp it settles at the
  # gradient that keeps J . n = 0 on the faceted shell surface, where the
  # discrete E_s has a small normal component, rather than at zero.
  [da_dt_norm]
    type = MFEMVectorL2Error
    variable = vessel_da_dt
    function = '0 0 0'
    execute_on = TIMESTEP_END
  []
  [e_source_norm]
    type = MFEMVectorL2Error
    variable = e_source
    function = '0 0 0'
    execute_on = TIMESTEP_END
  []
  [radial_force]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = force_density
    dual_variable = r_hat
    execute_on = TIMESTEP_END
  []
  [vertical_force]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = force_density
    dual_variable = z_hat
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  [csv]
    type = CSV
  []
  [vessel]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TorusVessel
    submesh = vessel
  []
[]
