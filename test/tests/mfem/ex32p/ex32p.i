# MFEM Example 32 (examples/ex32p.cpp)
#
# Computes the five lowest nonzero eigenmodes of the Maxwell eigenproblem,
# curl(curl(E)) = lambda epsilon E, on the inline-quad mesh with E x n = 0 on
# all boundaries. The mesh is two dimensional but E has three components: it
# is discretized with first order Nedelec elements for 3D vectors in 2D (ND_R2D
# in MFEM, range_dim = 3 here). The dielectric tensor epsilon is the constant,
# anisotropic 3x3 matrix of ex32p, with 2 on the diagonal and 1/sqrt(2) next to
# it. As in ex32p, the mesh is refined twice in serial and once in parallel, and
# the eigenmodes are found with AME preconditioned by AMS. The zero tangential
# condition is given as a two component vector, since MFEM rejects a three
# component one on the boundary of this space; zero values project to zero
# either way. The eigenvalues are compared against the values computed by
# ex32p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/inline-quad.mesh
  serial_refine = 2
  parallel_refine = 1
[]

[Problem]
  type = MFEMEigenproblem
  num_modes = 5
  rhs_matrix_coefficient = '2 0.7071067811865476 0; 0.7071067811865476 2 0.7071067811865476; 0 0.7071067811865476 2'
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    range_dim = 3
  []
[]

[Variables]
  [E]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[BCs]
  [boundary]
    type = MFEMVectorTangentialDirichletBC
    variable = E
    vector_coefficient = '0 0'
  []
[]

[Kernels]
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = E
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    singular = true
    print_level = 0
  []
  [ame]
    type = MFEMHypreAME
    preconditioner = ams
    l_tol = 1e-8
    l_max_its = 100
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[VectorPostprocessors]
  [eigenvalues]
    type = MFEMEigenvaluesPostprocessor
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex32p
  []
[]
