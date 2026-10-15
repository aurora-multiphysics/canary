# MFEM Example 26 (examples/ex26p.cpp)
#
# Solves the Poisson problem -div(grad(u)) = 1 on the star mesh with u = 0 on all
# boundaries, as in ex1p, using CG preconditioned by p-multigrid. As in ex26p,
# the mesh is refined serially as many times as keeps it at or under 1,000
# elements (twice for the 20 elements of the star mesh), then twice more in
# parallel. The multigrid hierarchy has first, second and fourth order H1
# spaces on this mesh. The fourth order problem is solved with partial
# assembly. Each level above the coarsest is smoothed by a second order
# Chebyshev smoother, and the coarsest level is solved approximately by CG
# preconditioned by BoomerAMG. ex26p also partially assembles the second order
# level, but MFEMGeometricMultigridSolver currently crashes when a level other
# than the finest is partially assembled, so that level is fully assembled
# here. Each level has its own smoother object, since one object shared by
# several levels also crashes. The L2 norm of the solution is compared against
# the value computed by ex26p.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/star.mesh
  serial_refine = 2
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

[FESpaceHierarchies]
  [h1_hierarchy]
    type = MFEMFESpaceHierarchy
    fespace = H1FESpace
    refinements = '2 4'
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace_hierarchy = h1_hierarchy
  []
[]

[BCs]
  [dirichlet]
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
  [source]
    type = MFEMDomainLFKernel
    variable = u
    coefficient = 1
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    l_tol = 0
    l_max_its = 1
    print_level = -1
  []
  [coarse]
    type = MFEMCGSolver
    preconditioner = boomeramg
    l_tol = 1e-2
    l_abs_tol = 0
    l_max_its = 10
    print_level = -1
  []
  [chebyshev_second_order_level]
    type = MFEMOperatorChebyshevSmoother
    order = 2
  []
  [chebyshev_fourth_order_level]
    type = MFEMOperatorChebyshevSmoother
    order = 2
  []
  [pmg]
    type = MFEMGeometricMultigridSolver
    variable = u
    smoothers = 'chebyshev_second_order_level chebyshev_fourth_order_level'
    coarse_solver = coarse
    assembly_levels = 'legacy legacy partial'
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = pmg
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 2000
    print_level = 1
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
  assembly_level = partial
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
    file_base = OutputData/ex26p
  []
[]
