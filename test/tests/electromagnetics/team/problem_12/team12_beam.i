# TEAM Problem 12 (FELIX cantilevered beam): elastic motion of the beam
# ======================================================================
#
# Sub-app of team12.i, which is the input you actually run. The benchmark is described
# in doc/content/problems/team_problem_12.md.
#
# The beam is linear elastic and is driven by the Lorentz force density f = J x B that
# team12.i computes. Its root is held fixed over the 76 mm inside the clamp. The
# elastodynamic equation rho d^2u/dt^2 = div(sigma(u)) + f is second order in time, and
# MFEMTransient integrates first order systems, so it is solved for the displacement u
# and the velocity v together,
#
#   (rho dv/dt, w) + (sigma(u), grad(w)) = (f, w),
#   (du/dt, z) - (v, z) = 0.
#
# The momentum balance is tested with the displacement test functions so that the
# elasticity operator sits on the diagonal block of the system.
#
# The mechanics is solved on a submesh holding only the beam. The fields exchanged with
# team12.i are defined on the whole mesh, since that is where team12.i holds them, and are
# moved between the two with MFEMSubMeshTransfer.

beam_density = 2713 # kg/m^3, Table 5 of the LSTC test case document
youngs_modulus = 6.89e10 # Pa, Table 5 of the LSTC test case document
# Not given in the test case document; 0.33 is typical of aluminium alloys and is the
# value the LS-DYNA input deck distributed with it uses.
poissons_ratio = 0.33

[Mesh]
  type = MFEMFileMesh
  file = ./team12.e
[]

[Problem]
  type = MFEMProblem
[]

[SubMeshes]
  [beam]
    type = MFEMDomainSubMesh
    block = Beam
  []
[]

[FunctorMaterials]
  [Beam]
    type = MFEMGenericFunctorMaterial
    prop_names = 'lambda mu density'
    prop_values = '${fparse youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))}
                   ${fparse youngs_modulus / (2 * (1 + poissons_ratio))}
                   ${beam_density}'
  []
[]

[FESpaces]
  # Second order, since linear hexahedra lock in bending with a single element through
  # the thickness of the beam.
  [BeamH1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = SECOND
    range_dim = 3
    submesh = beam
  []
  [BeamL2FESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = FIRST
    range_dim = 3
    submesh = beam
  []
  # The spaces of the fields exchanged with team12.i, which must match the ones they are
  # defined on there.
  [VelocityFESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = SECOND
    range_dim = 3
  []
  [ForceFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = FIRST
    range_dim = 3
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = BeamH1FESpace
  []
  [velocity]
    type = MFEMVariable
    fespace = BeamH1FESpace
  []
[]

[AuxVariables]
  # Received from team12.i.
  [lorentz_force_density]
    type = MFEMVariable
    fespace = ForceFESpace
  []
  [beam_force_density]
    type = MFEMVariable
    fespace = BeamL2FESpace
  []
  # Sent to team12.i.
  [beam_velocity]
    type = MFEMVariable
    fespace = VelocityFESpace
  []
[]

[Kernels]
  # Momentum balance, tested with the displacement test functions.
  [elasticity]
    type = MFEMLinearElasticityKernel
    variable = displacement
    lambda = lambda
    mu = mu
  []
  [inertia]
    type = MFEMTimeDerivativeVectorMassKernel
    variable = displacement
    trial_variable = velocity
    coefficient = density
  []
  [lorentz_force]
    type = MFEMVectorDomainLFKernel
    variable = displacement
    vector_coefficient = beam_force_density
  []
  # Kinematic relation v = du/dt, tested with the velocity test functions.
  [displacement_rate]
    type = MFEMTimeDerivativeVectorMassKernel
    variable = velocity
    trial_variable = displacement
  []
  [velocity]
    type = MFEMVectorMassKernel
    variable = velocity
    coefficient = -1.0
  []
[]

[BCs]
  [clamped_displacement]
    type = MFEMVectorDirichletBC
    variable = displacement
    boundary = clamp
  []
  [clamped_velocity]
    type = MFEMVectorDirichletBC
    variable = velocity
    boundary = clamp
  []
[]

[Transfers]
  # The force for this step, ahead of the solve.
  [force_to_beam]
    type = MFEMSubMeshTransfer
    from_variable = lorentz_force_density
    to_variable = beam_force_density
    execute_on = TIMESTEP_BEGIN
  []
  [velocity_from_beam]
    type = MFEMSubMeshTransfer
    from_variable = velocity
    to_variable = beam_velocity
  []
[]

[Solvers]
  # The block system coupling displacement and velocity is indefinite.
  [main]
    type = MFEMMUMPS
  []
[]

[Executioner]
  type = MFEMTransient
  # The time step is set by team12.i.
  dt = 1e-4
[]

[VectorPostprocessors]
  # On the centreline of the beam, at the free end and 167 mm from it: the positions of
  # the displacement measurements.
  [deflection_points]
    type = MFEMVariablePointValueSampler
    variable = displacement
    points = '0.487 0 0
              0.320 0 0'
  []
[]

[Postprocessors]
  [TipDeflection]
    type = VectorPostprocessorComponent
    vectorpostprocessor = deflection_points
    vector_name = displacement_1
    index = 0
  []
  [StationDeflection]
    type = VectorPostprocessorComponent
    vectorpostprocessor = deflection_points
    vector_name = displacement_1
    index = 1
  []
[]

[Outputs]
  # The deflections are reported by team12.i.
  console = false
[]
