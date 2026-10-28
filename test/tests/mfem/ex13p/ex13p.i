# MFEM Example 13 (examples/ex13p.cpp)
#
# Computes the five lowest nonzero eigenmodes of the Maxwell eigenproblem,
# curl(curl(E)) = lambda E, on the beam-tet mesh with E x n = 0 on all
# boundaries. As in ex13p, the mesh is refined twice in serial and once in
# parallel, and the problem is discretized with first order Nedelec elements.
# The eigenmodes are found with AME preconditioned by AMS, with the settings of
# ex13p. The eigenvalues are compared against the values computed by ex13p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/beam-tet.mesh
  serial_refine = 2
  parallel_refine = 1
[]

[Problem]
  type = MFEMEigenproblem
  num_modes = 5
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
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
    vector_coefficient = '0 0 0'
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
    file_base = OutputData/ex13p
  []
[]
