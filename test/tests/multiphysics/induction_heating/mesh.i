# Mesh for the induction heating example
# ======================================
#
# Builds the two meshes the example needs, from one common cross section so
# that the billet occupies exactly the same space in both:
#
#   induction_em.e      the full electromagnetic domain (billet, coil, air)
#   induction_billet.e  the billet on its own, for the thermal solve
#
# Regenerate both with
#
#   <exe> -i mesh.i --mesh-only induction_em.e
#   <exe> -i mesh.i Mesh/final_generator=billet --mesh-only induction_billet.e
#
# Geometry (axis along z, all dimensions in metres):
#
#            r=0.10   outer air boundary
#     .--------------------------------.
#     |                                |
#     |     +--+              +--+     |  coil:   0.030 < r < 0.045
#     |     |##|  .--------.  |##|     |          -0.05 < z < 0.05
#     |     |##|  | billet |  |##|     |
#     |     |##|  |  r<    |  |##|     |  billet: r < 0.020
#     |     |##|  | 0.020  |  |##|     |          -0.10 < z < 0.10
#     |     +--+  '--------'  +--+     |
#     |                                |
#     '--------------------------------'

billet_radius = 0.02
coil_inner_radius = 0.03
coil_outer_radius = 0.045
domain_radius = 0.10

billet_half_length = 0.10
coil_half_length = 0.05
domain_half_length = 0.15

[Mesh]
  # The cross section: four concentric rings, which the extrusion below turns
  # into four nested cylinders. 'preserve_volumes = false' puts the nodes exactly
  # on the nominal radii, so the billet boundary is identical in both meshes.
  [cross_section]
    type = ConcentricCircleMeshGenerator
    num_sectors = 6 # per quadrant, so 24 azimuthally
    radii = '${billet_radius} ${coil_inner_radius} ${coil_outer_radius} ${domain_radius}'
    rings = '4 1 2 4'
    has_outer_square = false
    preserve_volumes = false
    smoothing_max_it = 3
  []

  # Extrude along the axis. The layer boundaries are placed at +/- the billet and
  # coil half lengths so that both can be cut out by bounding box below.
  [extrude]
    type = AdvancedExtruderGenerator
    input = cross_section
    direction = '0 0 1'
    bottom_boundary = 1000
    top_boundary = 2000
    heights = '0.05 0.05 0.10 0.05 0.05'
    num_layers = '1 2 4 2 1'
  []
  [shift]
    type = TransformGenerator
    input = extrude
    transform = TRANSLATE
    vector_value = '0 0 -${domain_half_length}'
  []

  # Carve the billet and the coil out of their rings, over the axial extent each
  # occupies; every ring element left over is air.
  [billet_block]
    type = SubdomainBoundingBoxGenerator
    input = shift
    block_id = 10
    block_name = billet
    restricted_subdomains = '1'
    bottom_left = '-${billet_radius} -${billet_radius} -${billet_half_length}'
    top_right = '${billet_radius} ${billet_radius} ${billet_half_length}'
  []
  [coil_block]
    type = SubdomainBoundingBoxGenerator
    input = billet_block
    block_id = 20
    block_name = coil
    restricted_subdomains = '3'
    bottom_left = '-${coil_outer_radius} -${coil_outer_radius} -${coil_half_length}'
    top_right = '${coil_outer_radius} ${coil_outer_radius} ${coil_half_length}'
  []
  [air_block]
    type = RenameBlockGenerator
    input = coil_block
    old_block = '1 2 3 4'
    new_block = 'air air air air'
  []

  # The electromagnetic mesh is the whole thing. Only its outer surface carries a
  # boundary condition, so no interior sidesets are added to it.
  [em]
    type = RenameBoundaryGenerator
    input = air_block
    old_boundary = '1000 2000'
    new_boundary = 'domain_bottom domain_top'
  []

  # The thermal mesh is the billet alone. Its curved side is where the coolant
  # channel takes heat away; the two flat ends are left adiabatic.
  [billet_only]
    type = BlockDeletionGenerator
    input = air_block
    block = 'air coil'
  []
  [billet_surface]
    type = SideSetsAroundSubdomainGenerator
    input = billet_only
    block = billet
    new_boundary = billet_surface
  []
  [billet_side]
    type = ParsedGenerateSideset
    input = billet_surface
    included_boundaries = billet_surface
    combinatorial_geometry = 'abs(z) < ${billet_half_length} - 1e-8'
    new_sideset_name = billet_side
  []
  [billet]
    type = ParsedGenerateSideset
    input = billet_side
    included_boundaries = billet_surface
    combinatorial_geometry = 'abs(z) > ${billet_half_length} - 1e-8'
    new_sideset_name = billet_ends
  []

  final_generator = em
[]
