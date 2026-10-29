# MFEM Example 11 (examples/ex11p.cpp)
#
# Computes the five lowest eigenmodes of the Laplacian, -div(grad(u)) = lambda u,
# on the star mesh with homogeneous Dirichlet conditions on all boundaries. As in
# ex11p, the mesh is refined twice in serial and once in parallel, and the
# problem is discretized with first order H1 elements. The eigenmodes are found
# with LOBPCG preconditioned by BoomerAMG, with the settings of ex11p. The
# eigenvalues are compared against the values computed by ex11p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/star.mesh
  serial_refine = 2
  parallel_refine = 1
[]

[Problem]
  type = MFEMEigenproblem
  num_modes = 5
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [boundary]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = 0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [lobpcg]
    type = MFEMHypreLOBPCG
    preconditioner = boomeramg
    l_tol = 1e-8
    l_max_its = 200
    print_level = 1
    random_seed = 75
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
    file_base = OutputData/ex11p
  []
[]
