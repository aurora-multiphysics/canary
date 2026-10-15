# MFEM Example 22 (examples/ex22p.cpp)
#
# Solves the complex-valued damped harmonic oscillator
# -div(a grad(u)) - omega^2 b u + i omega c u = 0 on the inline-quad mesh, with
# a = 1/mu, b = epsilon and c = sigma, driven by the exact solution
# u = exp(-i kappa y), kappa = sqrt(mu omega (epsilon omega - i sigma)), on all
# boundaries. As in ex22p (problem 0), the mesh is refined once in serial and
# once in parallel, and the problem is discretized with first order H1
# elements. ex22p solves with FGMRES and a block diagonal preconditioner built
# from a separate real operator, which MOOSE cannot form, so the system is
# solved directly instead. The L2 error of the complex solution is compared
# against the value computed by ex22p.

mu = 1.0
epsilon = 1.0
sigma = 20.0
omega = 10.0
kappa_i = '${fparse -sqrt(mu*omega*sqrt(epsilon^2*omega^2+sigma^2)-epsilon*mu*omega^2)/sqrt(2)}'
kappa_r = '${fparse -mu*omega*sigma/kappa_i/2}'

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/inline-quad.mesh
  serial_refine = 1
  parallel_refine = 1
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
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
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [u_exact_real]
    type = ParsedFunction
    expression = 'exp(y * ${kappa_i}) * cos(y * ${kappa_r})'
  []
  [u_exact_imag]
    type = ParsedFunction
    expression = '-exp(y * ${kappa_i}) * sin(y * ${kappa_r})'
  []
[]

[BCs]
  [dirichlet]
    type = MFEMComplexScalarDirichletBC
    variable = u
    coefficient_real = u_exact_real
    coefficient_imag = u_exact_imag
  []
[]

[Kernels]
  [stiffness]
    type = MFEMComplexKernel
    variable = u
    [RealComponent]
      type = MFEMDiffusionKernel
      coefficient = '${fparse 1 / mu}'
    []
    [ImagComponent]
      type = MFEMDiffusionKernel
      coefficient = 0
    []
  []
  [mass_and_loss]
    type = MFEMComplexKernel
    variable = u
    [RealComponent]
      type = MFEMMassKernel
      coefficient = '${fparse -omega * omega * epsilon}'
    []
    [ImagComponent]
      type = MFEMMassKernel
      coefficient = '${fparse omega * sigma}'
    []
  []
[]

[Solvers]
  [direct]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [l2_error]
    type = MFEMComplexL2Error
    variable = u
    function_real = u_exact_real
    function_imag = u_exact_imag
  []
[]

[Outputs]
  [CSV]
    type = CSV
    file_base = OutputData/ex22p
  []
[]
