# TEAM Problem 12 (FELIX cantilevered beam): mesh
# ================================================
#
# Generates team12.e with
#
#   moose_test-opt -i team12_mesh.i --mesh-only team12.e
#
# The aluminium beam is a 487 x 3.175 x 100 mm plate lying along x, with its thickness
# along y and its width along z. Its root, 0 <= x <= 76 mm, is held in the clamp, so the
# free length is 411 mm. Every feature of the problem is a box aligned with the axes, so
# the whole domain is a tensor-product grid of hexahedra: the beam is meshed uniformly
# and the surrounding free space grows away from it, one element per interval.
#
# Node planes are placed at the edge of the clamp (x = 76 mm), at the measurement
# stations 167 mm from the free end (x = 320 mm) and at the free end (x = 487 mm), and at
# mid-length (x = 243.5 mm), about which the circulating current is evaluated.
#
# Dimensions are in metres.

[Mesh]
  [grid]
    type = CartesianMeshGenerator
    dim = 3
    # x: 400 mm of free space behind the root, the beam (clamp | to mid-length | to the
    # 320 mm station | to the tip, at about 10 mm per element) and 400 mm beyond the tip.
    dx = '0.15 0.10 0.07 0.05 0.03 0.076 0.1675 0.0765 0.167 0.03 0.05 0.07 0.10 0.15'
    ix = '1    1    1    1    1    8     17     8      17    1    1    1    1    1'
    # y: 300 mm of free space either side of the 3.175 mm beam, which has one element
    # through its thickness. That is enough for both physics: the beam is far thinner
    # than the skin depth at the rate the field decays (about 20 mm), and the second
    # order displacement represents the linear variation of the bending strain through
    # the thickness. A second layer changes the static deflection by less than 0.2%.
    dy = '0.15 0.08 0.04 0.02 0.01 0.0044125 0.003175 0.0044125 0.01 0.02 0.04 0.08 0.15'
    iy = '1    1    1    1    1    1         1        1         1    1    1    1    1'
    # z: 400 mm of free space either side of the 100 mm wide beam.
    dz = '0.15 0.10 0.07 0.05 0.03 0.1 0.03 0.05 0.07 0.10 0.15'
    iz = '1    1    1    1    1    10  1    1    1    1    1'
  []
  [shift]
    type = TransformGenerator
    input = grid
    transform = TRANSLATE
    vector_value = '-0.4 -0.306 -0.45'
  []
  [beam]
    type = SubdomainBoundingBoxGenerator
    input = shift
    block_id = 1
    bottom_left = '0 -0.0015875 -0.05'
    top_right = '0.487 0.0015875 0.05'
  []
  # MFEM reserves ID 0, so block and sideset IDs are renumbered from 1 before naming.
  [block_ids]
    type = RenameBlockGenerator
    input = beam
    old_block = '0 1'
    new_block = '1 2'
  []
  [block_names]
    type = RenameBlockGenerator
    input = block_ids
    old_block = '1 2'
    new_block = 'Vacuum Beam'
  []
  [boundary_ids]
    type = RenameBoundaryGenerator
    input = block_names
    old_boundary = 'back bottom right top left front'
    new_boundary = '1 2 3 4 5 6'
  []
  [boundary_names]
    type = RenameBoundaryGenerator
    input = boundary_ids
    old_boundary = '1 2 3 4 5 6'
    new_boundary = 'back bottom right top left front'
  []
  # The faces of the beam inside the clamp, on which its displacement is held at zero.
  [clamp]
    type = ParsedGenerateSideset
    input = boundary_names
    combinatorial_geometry = 'x < 0.0761'
    included_subdomains = Beam
    included_neighbors = Vacuum
    new_sideset_name = clamp
  []
[]
