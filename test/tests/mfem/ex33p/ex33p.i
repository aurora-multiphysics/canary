# MFEM Example 33 (examples/ex33p.cpp)
#
# Solves the fractional Poisson problem (-div(grad))^alpha u = 1 on the star
# mesh with u = 0 on all boundaries, for alpha = 0.5. As in ex33p, the inverse
# fractional power is replaced by a rational approximation,
#
#   (-div(grad))^(-alpha) ~ sum_i c_i (-div(grad) + d_i)^(-1),
#
# so u is the sum of c_i x_i over 12 independent reaction-diffusion problems
#
#   -div(grad(x_i)) + d_i x_i = 1,   x_i = 0 on the boundary.
#
# The coefficients c_i and shifts d_i are those ex33p uses when MFEM is built
# without LAPACK (examples/ex33.hpp); with LAPACK, ex33p computes them at run
# time with the AAA algorithm. Each problem is its own weak form, solved by its
# own problem composer with CG preconditioned by BoomerAMG, as in ex33p. The
# mesh is refined three times in serial, as in ex33p. The L2 norm of u is
# compared against the value computed by ex33p.

c0 = 2.290262e+02
c1 = 2.641819e+01
c2 = 1.005566e+01
c3 = 5.390411e+00
c4 = 3.340725e+00
c5 = 2.211205e+00
c6 = 1.508883e+00
c7 = 1.049474e+00
c8 = 7.462709e-01
c9 = 5.482686e-01
c10 = 4.232510e-01
c11 = 3.578967e-01

d0 = 3.168211e+04
d1 = 3.236077e+03
d2 = 9.868287e+02
d3 = 3.945597e+02
d4 = 1.738889e+02
d5 = 7.925178e+01
d6 = 3.624992e+01
d7 = 1.629196e+01
d8 = 6.982956e+00
d9 = 2.679984e+00
d10 = 7.782607e-01
d11 = 7.649166e-02

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/star.mesh
  serial_refine = 3
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
  [x0]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x1]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x2]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x3]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x4]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x5]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x6]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x7]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x8]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x9]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x10]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [x11]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxKernels]
  [sum]
    type = MFEMSumAux
    variable = u
    source_variables = 'x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11'
    scale_factors = '${c0} ${c1} ${c2} ${c3} ${c4} ${c5} ${c6} ${c7} ${c8} ${c9} ${c10} ${c11}'
  []
[]

[BCs]
  [dirichlet_0]
    type = MFEMScalarDirichletBC
    variable = x0
    coefficient = 0
  []
  [dirichlet_1]
    type = MFEMScalarDirichletBC
    variable = x1
    coefficient = 0
  []
  [dirichlet_2]
    type = MFEMScalarDirichletBC
    variable = x2
    coefficient = 0
  []
  [dirichlet_3]
    type = MFEMScalarDirichletBC
    variable = x3
    coefficient = 0
  []
  [dirichlet_4]
    type = MFEMScalarDirichletBC
    variable = x4
    coefficient = 0
  []
  [dirichlet_5]
    type = MFEMScalarDirichletBC
    variable = x5
    coefficient = 0
  []
  [dirichlet_6]
    type = MFEMScalarDirichletBC
    variable = x6
    coefficient = 0
  []
  [dirichlet_7]
    type = MFEMScalarDirichletBC
    variable = x7
    coefficient = 0
  []
  [dirichlet_8]
    type = MFEMScalarDirichletBC
    variable = x8
    coefficient = 0
  []
  [dirichlet_9]
    type = MFEMScalarDirichletBC
    variable = x9
    coefficient = 0
  []
  [dirichlet_10]
    type = MFEMScalarDirichletBC
    variable = x10
    coefficient = 0
  []
  [dirichlet_11]
    type = MFEMScalarDirichletBC
    variable = x11
    coefficient = 0
  []
[]

[Kernels]
  [diffusion_0]
    type = MFEMDiffusionKernel
    variable = x0
  []
  [mass_0]
    type = MFEMMassKernel
    variable = x0
    coefficient = ${d0}
  []
  [source_0]
    type = MFEMDomainLFKernel
    variable = x0
    coefficient = 1
  []
  [diffusion_1]
    type = MFEMDiffusionKernel
    variable = x1
  []
  [mass_1]
    type = MFEMMassKernel
    variable = x1
    coefficient = ${d1}
  []
  [source_1]
    type = MFEMDomainLFKernel
    variable = x1
    coefficient = 1
  []
  [diffusion_2]
    type = MFEMDiffusionKernel
    variable = x2
  []
  [mass_2]
    type = MFEMMassKernel
    variable = x2
    coefficient = ${d2}
  []
  [source_2]
    type = MFEMDomainLFKernel
    variable = x2
    coefficient = 1
  []
  [diffusion_3]
    type = MFEMDiffusionKernel
    variable = x3
  []
  [mass_3]
    type = MFEMMassKernel
    variable = x3
    coefficient = ${d3}
  []
  [source_3]
    type = MFEMDomainLFKernel
    variable = x3
    coefficient = 1
  []
  [diffusion_4]
    type = MFEMDiffusionKernel
    variable = x4
  []
  [mass_4]
    type = MFEMMassKernel
    variable = x4
    coefficient = ${d4}
  []
  [source_4]
    type = MFEMDomainLFKernel
    variable = x4
    coefficient = 1
  []
  [diffusion_5]
    type = MFEMDiffusionKernel
    variable = x5
  []
  [mass_5]
    type = MFEMMassKernel
    variable = x5
    coefficient = ${d5}
  []
  [source_5]
    type = MFEMDomainLFKernel
    variable = x5
    coefficient = 1
  []
  [diffusion_6]
    type = MFEMDiffusionKernel
    variable = x6
  []
  [mass_6]
    type = MFEMMassKernel
    variable = x6
    coefficient = ${d6}
  []
  [source_6]
    type = MFEMDomainLFKernel
    variable = x6
    coefficient = 1
  []
  [diffusion_7]
    type = MFEMDiffusionKernel
    variable = x7
  []
  [mass_7]
    type = MFEMMassKernel
    variable = x7
    coefficient = ${d7}
  []
  [source_7]
    type = MFEMDomainLFKernel
    variable = x7
    coefficient = 1
  []
  [diffusion_8]
    type = MFEMDiffusionKernel
    variable = x8
  []
  [mass_8]
    type = MFEMMassKernel
    variable = x8
    coefficient = ${d8}
  []
  [source_8]
    type = MFEMDomainLFKernel
    variable = x8
    coefficient = 1
  []
  [diffusion_9]
    type = MFEMDiffusionKernel
    variable = x9
  []
  [mass_9]
    type = MFEMMassKernel
    variable = x9
    coefficient = ${d9}
  []
  [source_9]
    type = MFEMDomainLFKernel
    variable = x9
    coefficient = 1
  []
  [diffusion_10]
    type = MFEMDiffusionKernel
    variable = x10
  []
  [mass_10]
    type = MFEMMassKernel
    variable = x10
    coefficient = ${d10}
  []
  [source_10]
    type = MFEMDomainLFKernel
    variable = x10
    coefficient = 1
  []
  [diffusion_11]
    type = MFEMDiffusionKernel
    variable = x11
  []
  [mass_11]
    type = MFEMMassKernel
    variable = x11
    coefficient = ${d11}
  []
  [source_11]
    type = MFEMDomainLFKernel
    variable = x11
    coefficient = 1
  []
[]

[WeakForms]
  [term_0]
    type = MFEMWeakForm
    kernels = 'diffusion_0 mass_0 source_0'
    bcs = 'dirichlet_0'
  []
  [term_1]
    type = MFEMWeakForm
    kernels = 'diffusion_1 mass_1 source_1'
    bcs = 'dirichlet_1'
  []
  [term_2]
    type = MFEMWeakForm
    kernels = 'diffusion_2 mass_2 source_2'
    bcs = 'dirichlet_2'
  []
  [term_3]
    type = MFEMWeakForm
    kernels = 'diffusion_3 mass_3 source_3'
    bcs = 'dirichlet_3'
  []
  [term_4]
    type = MFEMWeakForm
    kernels = 'diffusion_4 mass_4 source_4'
    bcs = 'dirichlet_4'
  []
  [term_5]
    type = MFEMWeakForm
    kernels = 'diffusion_5 mass_5 source_5'
    bcs = 'dirichlet_5'
  []
  [term_6]
    type = MFEMWeakForm
    kernels = 'diffusion_6 mass_6 source_6'
    bcs = 'dirichlet_6'
  []
  [term_7]
    type = MFEMWeakForm
    kernels = 'diffusion_7 mass_7 source_7'
    bcs = 'dirichlet_7'
  []
  [term_8]
    type = MFEMWeakForm
    kernels = 'diffusion_8 mass_8 source_8'
    bcs = 'dirichlet_8'
  []
  [term_9]
    type = MFEMWeakForm
    kernels = 'diffusion_9 mass_9 source_9'
    bcs = 'dirichlet_9'
  []
  [term_10]
    type = MFEMWeakForm
    kernels = 'diffusion_10 mass_10 source_10'
    bcs = 'dirichlet_10'
  []
  [term_11]
    type = MFEMWeakForm
    kernels = 'diffusion_11 mass_11 source_11'
    bcs = 'dirichlet_11'
  []
[]

[ProblemComposers]
  [term_0]
    type = MFEMWeakFormProblemComposer
    weak_form = term_0
  []
  [term_1]
    type = MFEMWeakFormProblemComposer
    weak_form = term_1
  []
  [term_2]
    type = MFEMWeakFormProblemComposer
    weak_form = term_2
  []
  [term_3]
    type = MFEMWeakFormProblemComposer
    weak_form = term_3
  []
  [term_4]
    type = MFEMWeakFormProblemComposer
    weak_form = term_4
  []
  [term_5]
    type = MFEMWeakFormProblemComposer
    weak_form = term_5
  []
  [term_6]
    type = MFEMWeakFormProblemComposer
    weak_form = term_6
  []
  [term_7]
    type = MFEMWeakFormProblemComposer
    weak_form = term_7
  []
  [term_8]
    type = MFEMWeakFormProblemComposer
    weak_form = term_8
  []
  [term_9]
    type = MFEMWeakFormProblemComposer
    weak_form = term_9
  []
  [term_10]
    type = MFEMWeakFormProblemComposer
    weak_form = term_10
  []
  [term_11]
    type = MFEMWeakFormProblemComposer
    weak_form = term_11
  []
[]

[Solvers]
  # One solver serves all 12 problem composers: each solve sets the solver's operator to that
  # problem's matrix, which also sets up BoomerAMG for it. With several weak forms a solver must name
  # one of them for context, which these solvers do not otherwise use.
  [boomeramg]
    type = MFEMHypreBoomerAMG
    weak_form = term_0
    print_level = -1
  []
  [cg]
    type = MFEMCGSolver
    weak_form = term_0
    preconditioner = boomeramg
    use_initial_guess = false
    l_tol = 1e-12
    l_abs_tol = 0
    l_max_its = 2000
    print_level = 3
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
    file_base = OutputData/ex33p
  []
[]
