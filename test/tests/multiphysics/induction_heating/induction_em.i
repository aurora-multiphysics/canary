# Induction heating: the electromagnetic sub-app (MFEM)
# =====================================================
#
# Time-harmonic eddy currents in a cylindrical billet driven by a solenoid coil,
# solved with MFEM on the mesh built by mesh.i. This is the sub-app of
# induction_heating.i, which solves for the billet temperature; see
# doc/content/problems/induction_heating.md.
#
# Formulation
# -----------
# Writing the magnetic vector potential as A(x, t) = Re{(A_r + j A_i) e^{j omega t}},
# the eddy current A-formulation
#
#   curl(nu curl A) + j omega sigma A = J_s
#
# splits into a pair of real equations for the in-phase and quadrature parts,
#
#   curl(nu curl A_r) - omega sigma A_i = J_s
#   curl(nu curl A_i) + omega sigma A_r = 0
#
# with nu = 1/mu the magnetic reluctivity and J_s the source current density
# prescribed in the coil. Each is discretised in first order Nedelec elements;
# the off-diagonal loss terms couple them into one block system.
#
# MFEM can express this directly with MFEMComplexVariable instead, but a complex
# MFEMProblem may not hold real variables and the MFEM<->libMesh transfers carry
# only real ones, so the split form is what lets this app be coupled to a libMesh
# parent. It solves the same equations either way.
#
# The electric field is E = -j omega A, so the current density induced in the
# billet is J = sigma E and the cycle-averaged Ohmic dissipation is
#
#   q = |J|^2 / (2 sigma) = sigma |E|^2 / 2 = sigma omega^2 (|A_r|^2 + |A_i|^2) / 2
#
# which is the quantity handed to the thermal parent app.
#
# The curl-curl operator has a large gradient nullspace wherever sigma = 0, which
# would leave the system singular in the air. A small artificial conductivity
# (conductivity_air below, 1 S/m against 1.4e6 S/m in the billet) regularises it;
# its skin depth is several metres, so it does not perturb the solution.

frequency = 3000 # Hz
angular_frequency = ${fparse 2 * pi * frequency}

mu0 = ${fparse 4e-7 * pi}
reluctivity = ${fparse 1 / mu0}

# Stainless steel near room temperature, and its temperature coefficient of
# resistivity. The conductivity falls as the billet heats, which is the feedback
# from the thermal solve back into this one.
conductivity_ref = 1.4e6 # S/m
conductivity_alpha = 3e-3 # 1/K
reference_temperature = 300 # K

# Gauge regularisation for the non-conducting regions; see above.
conductivity_air = 1.0 # S/m

# 15 kA-turns spread over the 0.015 x 0.1 m coil cross section.
coil_current = 15000 # A-turns
coil_area = ${fparse (0.045 - 0.03) * 0.1}
coil_current_density = ${fparse coil_current / coil_area}

[Mesh]
  type = MFEMFileMesh
  file = induction_em.e
[]

[Problem]
  type = MFEMProblem
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
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [A_real]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [A_imag]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  # Billet temperature, transferred in from the thermal parent app. Its initial
  # value matters: the first electromagnetic solve happens before any transfer.
  [temperature]
    type = MFEMVariable
    fespace = H1FESpace
  []
  # Cycle-averaged Ohmic dissipation, transferred out to the thermal parent app.
  [joule_heating]
    type = MFEMVariable
    fespace = L2FESpace
  []

  # The fields themselves, for viewing. E = -j omega A gives E_r = omega A_i and
  # E_i = -omega A_r; B = curl A componentwise.
  [e_field_real]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [e_field_imag]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [b_field_real]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [b_field_imag]
    type = MFEMVariable
    fespace = HDivFESpace
  []

  # Cycle root-mean-square magnitudes of the two fields, which are what the
  # postprocessors below average. A scalar L2 space so that they can be block
  # restricted to the billet.
  [e_field_rms]
    type = MFEMVariable
    fespace = L2FESpace
  []
  [b_field_rms]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[ICs]
  [temperature_ic]
    type = MFEMScalarIC
    variable = temperature
    coefficient = ${reference_temperature}
  []
[]

[Functions]
  # The azimuthal unit vector, which the solenoid current follows. The coil is a
  # true annulus, so this source is divergence free on it.
  [coil_current_density]
    type = ParsedVectorFunction
    expression_x = '-J * y / hypot(x, y)'
    expression_y = 'J * x / hypot(x, y)'
    expression_z = '0'
    symbol_names = 'J'
    symbol_values = '${coil_current_density}'
  []

  # sigma(T), linearised about the reference temperature. The denominator is
  # floored well above zero so that a runaway temperature cannot flip its sign.
  [billet_conductivity]
    type = MFEMParsedFunction
    expression = 'sigma_ref / max(1 + alpha * (T - T_ref), 0.1)'
    symbol_names = 'T sigma_ref alpha T_ref'
    symbol_values = 'temperature ${conductivity_ref} ${conductivity_alpha} ${reference_temperature}'
  []

  # The loss coefficient omega*sigma multiplying the off-diagonal mass terms, and
  # the cycle-averaged dissipation sigma*omega^2*|A|^2/2.
  [loss_coefficient]
    type = MFEMParsedFunction
    expression = 'omega * sigma'
    symbol_names = 'sigma omega'
    symbol_values = 'electrical_conductivity ${angular_frequency}'
  []
  [negative_loss_coefficient]
    type = MFEMParsedFunction
    expression = '-omega * sigma'
    symbol_names = 'sigma omega'
    symbol_values = 'electrical_conductivity ${angular_frequency}'
  []
  # Cycle root-mean-square magnitudes. For a phasor field
  # F(t) = Re{(F_r + j F_i) exp(j omega t)}, the mean of |F|^2 over a cycle is
  # (|F_r|^2 + |F_i|^2) / 2. MFEM declares the magnitude of a vector variable as
  # <name>_mag, and the magnitude of its curl as <name>_curl_mag, so both of
  # these read straight off the solution rather than off an auxvariable computed
  # by another AuxKernel, which would need the two to be ordered.
  [e_field_rms_density]
    type = MFEMParsedFunction
    expression = 'omega * sqrt(0.5 * (a_r * a_r + a_i * a_i))'
    symbol_names = 'omega a_r a_i'
    symbol_values = '${angular_frequency} A_real_mag A_imag_mag'
  []
  [b_field_rms_density]
    type = MFEMParsedFunction
    expression = 'sqrt(0.5 * (b_r * b_r + b_i * b_i))'
    symbol_names = 'b_r b_i'
    symbol_values = 'A_real_curl_mag A_imag_curl_mag'
  []

  [ohmic_dissipation]
    type = MFEMParsedFunction
    expression = '0.5 * sigma * omega * omega * (a_r * a_r + a_i * a_i)'
    symbol_names = 'sigma omega a_r a_i'
    symbol_values = 'electrical_conductivity ${angular_frequency} A_real_mag A_imag_mag'
  []
[]

[FunctorMaterials]
  [billet]
    type = MFEMGenericFunctorMaterial
    prop_names = 'electrical_conductivity magnetic_reluctivity'
    prop_values = 'billet_conductivity ${reluctivity}'
    block = billet
  []
  [coil_and_air]
    type = MFEMGenericFunctorMaterial
    prop_names = 'electrical_conductivity magnetic_reluctivity'
    prop_values = '${conductivity_air} ${reluctivity}'
    block = 'coil air'
  []
[]

[BCs]
  # A x n = 0 on the outer surface of the air box: the flux is contained, which is
  # a reasonable truncation at five billet radii.
  [tangential_A_real]
    type = MFEMVectorTangentialDirichletBC
    variable = A_real
    boundary = 'outer domain_bottom domain_top'
  []
  [tangential_A_imag]
    type = MFEMVectorTangentialDirichletBC
    variable = A_imag
    boundary = 'outer domain_bottom domain_top'
  []
[]

[Kernels]
  # Diagonal blocks: curl(nu curl A_r) and curl(nu curl A_i).
  [curl_curl_real]
    type = MFEMCurlCurlKernel
    variable = A_real
    coefficient = magnetic_reluctivity
  []
  [curl_curl_imag]
    type = MFEMCurlCurlKernel
    variable = A_imag
    coefficient = magnetic_reluctivity
  []

  # Off-diagonal blocks: the eddy current loss term -omega sigma A_i in the real
  # equation and +omega sigma A_r in the imaginary one.
  [loss_real]
    type = MFEMMixedVectorMassKernel
    variable = A_real
    trial_variable = A_imag
    coefficient = negative_loss_coefficient
  []
  [loss_imag]
    type = MFEMMixedVectorMassKernel
    variable = A_imag
    trial_variable = A_real
    coefficient = loss_coefficient
  []

  # The coil drives the in-phase equation, which fixes the phase reference.
  [coil_drive]
    type = MFEMVectorFEDomainLFKernel
    variable = A_real
    vector_coefficient = coil_current_density
    block = coil
  []
[]

[AuxKernels]
  [ohmic_heating]
    type = MFEMScalarProjectionAux
    variable = joule_heating
    coefficient = ohmic_dissipation
    execute_on = TIMESTEP_END
  []

  # E = -j omega A, so the in-phase part of E comes from the quadrature part of A
  # and vice versa, with a sign change.
  [electric_field_real]
    type = MFEMSumAux
    variable = e_field_real
    source_variables = A_imag
    scale_factors = ${angular_frequency}
    execute_on = TIMESTEP_END
  []
  [electric_field_imag]
    type = MFEMSumAux
    variable = e_field_imag
    source_variables = A_real
    scale_factors = ${fparse -angular_frequency}
    execute_on = TIMESTEP_END
  []

  # B = curl A, which takes the H(curl) potential into the H(div) space.
  [magnetic_flux_density_real]
    type = MFEMCurlAux
    variable = b_field_real
    source = A_real
    execute_on = TIMESTEP_END
  []
  [magnetic_flux_density_imag]
    type = MFEMCurlAux
    variable = b_field_imag
    source = A_imag
    execute_on = TIMESTEP_END
  []

  [electric_field_rms]
    type = MFEMScalarProjectionAux
    variable = e_field_rms
    coefficient = e_field_rms_density
    execute_on = TIMESTEP_END
  []
  [magnetic_flux_density_rms]
    type = MFEMScalarProjectionAux
    variable = b_field_rms
    coefficient = b_field_rms_density
    execute_on = TIMESTEP_END
  []
[]

[Solvers]
  # The two Nedelec fields assemble into one 2x2 block operator, for which the
  # H(curl) preconditioners (AMS) are not set up, so it is solved directly.
  #
  # MUMPS rather than SuperLU deliberately. The two give bit-identical answers
  # here, but where PETSc is built --with-cuda its SuperLU_DIST puts the
  # factorisation on the GPU and takes around 700 MB of device memory per MPI
  # rank, which exhausts an 8 GB card at about ten ranks and fails with
  # 'cudaErrorMemoryAllocation'. MUMPS stays on the host and uses none.
  [main]
    type = MFEMMUMPS
    print_level = 0
  []
[]

[Executioner]
  type = MFEMTransient
  assembly_level = legacy
  device = cpu
  # The eddy current problem has no memory of its own: at every step it is solved
  # afresh for the conductivity implied by the temperature most recently received,
  # and the step size comes from the parent app.
  #
  # Both of these are set large deliberately. A TransientMultiApp that is not
  # sub-cycling constrains its parent to the smallest timestep among itself and
  # its sub-apps, so a small dt here would silently override whatever the parent
  # asked for. This app has no timescale of its own to impose, so it declares one
  # far larger than any parent would use and lets the parent decide.
  dt = 1e6
  end_time = 1e6
[]

[Postprocessors]
  [mean_power_density]
    type = MFEMElementAverageValue
    variable = joule_heating
    block = billet
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Cycle RMS field magnitudes averaged over the billet, which is where they set
  # the heating. Drop the block restriction, or point it at 'coil' or 'air', to
  # look at the applied field instead.
  [mean_e_field]
    type = MFEMElementAverageValue
    variable = e_field_rms
    block = billet
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mean_b_field]
    type = MFEMElementAverageValue
    variable = b_field_rms
    block = billet
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [Values]
    type = CSV
    file_base = OutputData/InductionEM
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/InductionEM
    execute_on = FINAL
  []
[]
