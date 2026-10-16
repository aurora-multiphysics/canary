# MFEM Example 34 (examples/ex34p.cpp)
#
# Solves the magnetostatic problem curl(curl(A)) + delta A = J on the fichera
# mesh, with A x n = 0 on all boundaries except the symmetry planes (boundaries
# 9 to 16), where the natural condition holds. The current density J flows
# through a conductor made of five of the mesh's elements, given attribute 2 in
# fichera-mixed-conductor.mesh as ex34p does in code. J is computed on the
# conductor, a submesh, by the sub-application current_density.i, which itself
# runs electric_potential.i for the potential that drives the current. As in
# ex34p, the mesh is refined once in serial and once in parallel, the vector
# potential is discretized with first order Nedelec elements, delta = 1e-6, and
# the system is solved with PCG preconditioned by AMS. The L2 norms of A and of
# the magnetic flux B = curl(A) are compared against the values computed by
# ex34p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/fichera-mixed-conductor.mesh
  serial_refine = 1
  parallel_refine = 1
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
[]

[Variables]
  [A]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [J]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [B]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[AuxKernels]
  [curl]
    type = MFEMCurlAux
    variable = B
    source = A
  []
[]

[BCs]
  [tangential_A]
    type = MFEMVectorTangentialDirichletBC
    variable = A
    boundary = '1 2 3 4 5 6 7 8 17 18 19 20 21 22 23 24'
    vector_coefficient = '0 0 0'
  []
[]

[Kernels]
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = A
  []
  [mass]
    type = MFEMVectorFEMassKernel
    variable = A
    coefficient = 1e-6
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = A
    vector_coefficient = J
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    print_level = 0
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 500
    print_level = 2
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [current]
    type = FullSolveMultiApp
    input_files = current_density.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_current]
    type = MultiAppMFEMCopyTransfer
    source_variables = J
    variables = J
    from_multi_app = current
  []
[]

[Postprocessors]
  [vector_potential_l2_norm]
    type = MFEMVectorL2Error
    variable = A
    function = '0 0 0'
  []
  [magnetic_flux_l2_norm]
    type = MFEMVectorL2Error
    variable = B
    function = '0 0 0'
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex34p
  []
[]
