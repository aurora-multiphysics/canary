# MFEM Example 39 (examples/ex39p.cpp)
#
# Solves the Poisson problem -div(k grad(u)) = f on a compass rose read from a
# Gmsh file whose regions are named. MOOSE looks names up in the attribute sets
# MFEM reads from the file:
#
#   k = 1e-6 everywhere, plus 1 on "Base" and 2 on the "Even" half of each
#       compass point,
#   f = 1 on the "Even" halves and 0 elsewhere,
#   u = 0 on every boundary.
#
# ex39p combines the four "Even" sets into a new set named "Rose Even", and the
# eight boundary sets into one named "Boundary". MOOSE lists the sets instead.
# The "Even" sets are listed by attribute number ("N Even" = 10, "W Even" = 12,
# "S Even" = 14, "E Even" = 16), because MOOSE splits lists on spaces and so
# cannot take a name that contains one.
# As in ex39p, the mesh is refined four times in serial and twice in parallel.
# The system is solved with first order Lagrange elements and MFEM's CG
# preconditioned by BoomerAMG, and the L2 norm of the solution is compared
# against the value computed by ex39p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/compass.msh
  serial_refine = 4
  parallel_refine = 2
[]

[Problem]
  type = MFEMProblem
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
    boundary = 'ENE NNE NNW WNW WSW SSW SSE ESE'
    coefficient = 0
  []
[]

[Kernels]
  [default_diffusion]
    type = MFEMDiffusionKernel
    variable = u
    coefficient = 1e-6
  []
  [base_diffusion]
    type = MFEMDiffusionKernel
    variable = u
    block = 'Base'
    coefficient = 1
  []
  [rose_diffusion]
    type = MFEMDiffusionKernel
    variable = u
    block = '10 12 14 16'
    coefficient = 2
  []
  [source]
    type = MFEMDomainLFKernel
    variable = u
    block = '10 12 14 16'
    coefficient = 1
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 1
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = boomeramg
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 2000
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [solution_l2_norm]
    type = MFEML2Error
    variable = u
    function = 0
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex39p
  []
[]
