# Source current density for a homogenised stranded coil carrying a specified total current,
# for use by clf_Aformsolve.i and msqlf_Aformsolve.i.
#
# The coil is homogenised over its turns, so the current density has constant magnitude I/S
# across the cross-section, where I is the total current and S the cross-sectional area. The
# direction the current follows around the coil is obtained from the coil geometry alone: the
# electric field solved for below is everywhere tangential to the turns, so e/|e| is the unit
# vector following them, and the current density is
#
#   J = (I / S) e / |e|.
#
# S is measured as the flux of e/|e| through the cut surface, so no calibration constant relating
# the loop voltage to the total current is needed, and the loop voltage below only sets the sense
# in which the current circulates. Reversing its sign reverses e and S together, leaving J
# unchanged; the circulation direction is fixed by the orientation of the cut surface.

initial_coil_domains = 'First_half Second_half'
coil_cut_surface = 'Cut'
total_coil_current = 1260.0 # A
coil_loop_voltage = 1.0
coil_conductivity = 1.0

[Mesh]
  type = MFEMFileMesh
  file = ./team_coil_two_vols_plate_exterior_tet_m_offset_box.msh
[]

[Problem]
  type = MFEMProblem
[]

[FunctorMaterials]
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${coil_conductivity}
  []
  # Block restricted copy of the scaling below. A piecewise coefficient evaluates to zero on the
  # subdomains it has not been assigned to, so scaling the electric field by this gives a current
  # density that is zero outside the coil, rather than one of full magnitude in the layer of
  # exterior elements sharing edges, and so nonzero tangential degrees of freedom, with it.
  [CoilCurrentDensity]
    type = MFEMGenericFunctorMaterial
    prop_names = coil_current_density_scale
    prop_values = current_density_scale
    block = ${initial_coil_domains}
  []
[]

[Functions]
  # Reciprocal of |e| on the transition region, so that scaling e by it gives the unit vector
  # following the turns of the coil.
  [inverse_e_field_magnitude]
    type = MFEMParsedFunction
    expression = '1 / e_magnitude'
    symbol_names = 'e_magnitude'
    symbol_values = 'transition_e_field_mag'
  []
  # Scaling that turns the electric field into the homogenised stranded current density
  # J = (I / S) e / |e|, of constant magnitude I/S over the coil cross-section.
  [current_density_scale]
    type = MFEMParsedFunction
    expression = '${total_coil_current} / (area * e_magnitude)'
    symbol_names = 'area e_magnitude'
    symbol_values = 'CoilCrossSectionalArea total_e_field_mag'
  []
  # The same scaling, built from the copy of e restricted to the transition region, for
  # integrating the total current over the cut surface.
  [transition_current_density_scale]
    type = MFEMParsedFunction
    expression = '${total_coil_current} / (area * e_magnitude)'
    symbol_names = 'area e_magnitude'
    symbol_values = 'CoilCrossSectionalArea transition_e_field_mag'
  []
[]

[ICs]
  [coil_external_potential_ic]
    type = MFEMScalarBoundaryIC
    variable = coil_external_potential
    boundary = ${coil_cut_surface}
    coefficient = ${coil_loop_voltage}
  []
[]

[SubMeshes]
  [cut]
    type = MFEMCutTransitionSubMesh
    cut_boundary = ${coil_cut_surface}
    block = ${initial_coil_domains}
    transition_subdomain = transition_dom
    transition_subdomain_boundary = transition_bdr
    closed_subdomain = coil_dom
  []
  [coil]
    type = MFEMDomainSubMesh
    block = coil_dom
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [CoilH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = coil
  []
  [TransitionH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = cut
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
  [TransitionHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = cut
  []
[]

[Variables]
  [coil_induced_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
[]

[AuxVariables]
  [coil_external_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
  [transition_external_potential]
    type = MFEMVariable
    fespace = TransitionH1FESpace
  []
  [transition_external_e_field]
    type = MFEMVariable
    fespace = TransitionHCurlFESpace
  []
  # Total electric field restricted to the transition region, used to measure the coil
  # cross-sectional area on the cut surface bounding that region.
  [transition_e_field]
    type = MFEMVariable
    fespace = TransitionHCurlFESpace
  []
  [induced_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [induced_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [external_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [total_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  # Homogenised stranded current density, transferred out to the A-form solves.
  [source_j_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[AuxKernels]
  [update_induced_e_field]
    type = MFEMGradAux
    variable = induced_e_field
    source = induced_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_external_e_field]
    type = MFEMGradAux
    variable = transition_external_e_field
    source = transition_external_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_total_e_field]
    type = MFEMSumAux
    variable = total_e_field
    source_variables = 'induced_e_field external_e_field'
    execute_on = TIMESTEP_END
  []
  # Scale the direction field up to the homogenised stranded current density. The result is
  # projected onto an H(div) space: the current density is tangential to the surface of the coil,
  # so it is normal-continuous there and across the coil interior, and an H(div) projection
  # preserves the flux, and so the total current, through every face. An H(curl) space is not
  # suitable, as its tangential continuity cannot represent the jump at the coil surface.
  [update_source_j_field]
    type = MFEMScaledVectorAux
    variable = source_j_field
    vector_coefficient = total_e_field
    coefficient = coil_current_density_scale
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [CoilPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = conductivity
    dual_variable = total_e_field
    primal_variable = total_e_field
    block = 'First_half Second_half'
  []
  # The flux of the unit vector e/|e| through the cut surface is the coil cross-sectional area.
  [CoilCrossSectionalArea]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    coefficient = inverse_e_field_magnitude
    variable = transition_e_field
    boundary = ${coil_cut_surface}
  []
  # Total current carried by the source current density, which should recover the value
  # requested above.
  [TotalCoilCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    coefficient = transition_current_density_scale
    variable = transition_e_field
    boundary = ${coil_cut_surface}
    execute_on = TIMESTEP_END
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = coil_induced_potential
    coefficient = conductivity
  []
  [source]
    type = MFEMMixedGradGradKernel
    trial_variable = coil_external_potential
    variable = coil_induced_potential
    coefficient = conductivity
    block = 'transition_dom'
  []
[]

[Solvers]
  [superlu]
    type = MFEMMUMPS
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Transfers]
  [submesh_transfer_from_coil]
    type = MFEMSubMeshTransfer
    from_variable = coil_induced_potential
    to_variable = induced_potential
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_to_transition]
    type = MFEMSubMeshTransfer
    from_variable = coil_external_potential
    to_variable = transition_external_potential
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_from_transition]
    type = MFEMSubMeshTransfer
    from_variable = transition_external_e_field
    to_variable = external_e_field
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_total_to_transition]
    type = MFEMSubMeshTransfer
    from_variable = total_e_field
    to_variable = transition_e_field
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ClosedCoilSourceSubMesh
    vtk_format = ASCII
   submesh = coil
  []
  [GlobalParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/WholePotentialCoil
    vtk_format = ASCII
  []
[]
