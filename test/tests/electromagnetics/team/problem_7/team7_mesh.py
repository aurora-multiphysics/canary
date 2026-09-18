#!/usr/bin/env python3
"""Build the TEAM Problem 7 mesh with gmsh.

TEAM Problem 7 is a 3-D time-harmonic eddy current benchmark: a racetrack coil
carrying a sinusoidal current sits above a thick aluminium plate with an
off-centre square hole.  The benchmark and the geometry are documented in
doc/content/problems/team_problem_7.md.

Every part of the model is a prism standing on the x-y plane - the plate, the
hole through it and the coil - so the mesh is built the same way: the plan is
triangulated once and then extruded through a stack of slabs, giving wedge
elements throughout.  That puts layer boundaries exactly on the faces of the
plate and the coil and on the z = 34 mm plane the flux density is measured on,
and lets the plate be resolved finely through its 19 mm thickness, where the
skin depth sets the scale, without paying for the same resolution in plan.

The problem has no symmetry, so the whole of it is meshed.  Dimensions come
from Fig. 1 and Fig. 2 of the problem statement and are in millimetres here;
the mesh is written in metres.  Run with --help for the meshing controls.
"""

import argparse
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# Geometry, in millimetres, straight from Figs. 1 and 2 of the problem statement.
# ---------------------------------------------------------------------------
PLATE = 294.0  # the plate is a 294 x 294 mm square with a corner at the origin
PLATE_Z = 19.0  # and is 19 mm thick
HOLE_X0, HOLE_Y0, HOLE = 18.0, 18.0, 108.0  # 108 mm square hole, 18 mm in

COIL_CX, COIL_CY = 194.0, 100.0  # centre of the coil in plan
COIL_OUT_HALF, COIL_OUT_R = 100.0, 50.0  # outer outline: 200 mm square, R50 corners
COIL_IN_HALF, COIL_IN_R = 75.0, 25.0  # inner outline: 150 mm square, R25 corners
COIL_Z0, COIL_Z1 = 49.0, 149.0  # 30 mm above the plate, 100 mm tall

MEASURE_Z = 34.0  # the height of the lines A1-B1 and A2-B2 the field is measured on

# The outermost boundary of the recommended mesh (Fig. 2): a box centred in plan
# on the middle of the plate and running from 300 mm below it to 300 mm above the
# coil.  The benchmark asks for B.n = 0 there.
AIR_X0, AIR_X1 = -1353.0, 1647.0
AIR_Z0, AIR_Z1 = -300.0, 449.0

MM = 1.0e-3

# Comparison tolerance in millimetres, used for locating entities by position.
TOL = 1.0e-4


def graded(n, ratio, fine_at_top):
    """Cumulative layer heights, normalised to 1, growing geometrically by ratio.

    gmsh takes the heights of the layers of an extrusion as a cumulative list, in
    the direction of the extrusion.  The fine end is the one next to the device.
    """
    widths = [ratio**i for i in range(n)]
    if fine_at_top:
        widths.reverse()
    total = sum(widths)
    heights, acc = [], 0.0
    for w in widths:
        acc += w / total
        heights.append(acc)
    heights[-1] = 1.0
    return heights


def build(args):
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("team7")
    occ = gmsh.model.occ

    # -- the plan -----------------------------------------------------------
    # Everything is laid out on the bottom face of the free space box and then
    # extruded upwards.  fragment cuts the outlines into the pieces of the plan
    # that each lie wholly inside or outside every outline.
    z = AIR_Z0
    air = occ.addRectangle(AIR_X0, AIR_X0, z, AIR_X1 - AIR_X0, AIR_X1 - AIR_X0)
    plate = occ.addRectangle(0.0, 0.0, z, PLATE, PLATE)
    hole = occ.addRectangle(HOLE_X0, HOLE_Y0, z, HOLE, HOLE)
    coil_out = occ.addRectangle(COIL_CX - COIL_OUT_HALF, COIL_CY - COIL_OUT_HALF, z,
                                2 * COIL_OUT_HALF, 2 * COIL_OUT_HALF,
                                roundedRadius=COIL_OUT_R)
    coil_in = occ.addRectangle(COIL_CX - COIL_IN_HALF, COIL_CY - COIL_IN_HALF, z,
                               2 * COIL_IN_HALF, 2 * COIL_IN_HALF,
                               roundedRadius=COIL_IN_R)
    outlines = [air, plate, hole, coil_out, coil_in]
    _, fmap = occ.fragment([(2, t) for t in outlines], [])
    occ.synchronize()

    # fragment reports, for each input outline, the pieces it ended up as, so
    # whether a piece is inside an outline is read straight off that rather than
    # guessed at from its position.
    inside = [set(t for _, t in ents) for ents in fmap]
    pieces = sorted(inside[0])
    in_plate = inside[1] - inside[2]  # the plate, less its hole
    in_coil = inside[3] - inside[4]  # the winding, between the two outlines

    # -- the slabs ----------------------------------------------------------
    # (top of slab, number of layers, layer heights).  The plate carries the
    # finest layers: the skin depth in the aluminium is 12 mm at 50 Hz and 6 mm at
    # 200 Hz, against a plate 19 mm thick.  The free space above and below is
    # graded away from the device.
    slabs = [
        (0.0, args.far_layers, graded(args.far_layers, args.grading, fine_at_top=True)),
        (PLATE_Z, args.plate_layers, None),
        (MEASURE_Z, args.gap_layers, None),
        (COIL_Z0, args.gap_layers, None),
        (COIL_Z1, args.coil_layers, None),
        (AIR_Z1, args.far_layers, graded(args.far_layers, args.grading, fine_at_top=False)),
    ]

    plate_vols, coil_vols, air_vols = [], [], []
    bottom = [(2, t) for t in pieces]
    z0 = AIR_Z0
    for z1, n, heights in slabs:
        # Uniform slabs are n layers in one go; graded ones are n single layers at
        # the cumulative heights given. recombine makes the layers wedges rather
        # than splitting them into tetrahedra.
        if heights is None:
            num_elements, heights = [n], [1.0]
        else:
            num_elements = [1] * n
        out = occ.extrude(bottom, 0.0, 0.0, z1 - z0, numElements=num_elements,
                          heights=heights, recombine=True)
        # The extrusion returns, for each surface in order, its top face, then the
        # volume swept out, then the lateral faces, whose number varies. The volume
        # is the only 3-D entity in each group, so it marks where a group starts.
        tops, vols = [], []
        for i, (dim, tag) in enumerate(out):
            if dim == 3:
                tops.append(out[i - 1])
                vols.append(tag)
        if len(vols) != len(pieces):
            sys.exit(f"expected {len(pieces)} volumes in the slab below z = {z1}, "
                     f"got {len(vols)}")
        for piece, vol in zip(pieces, vols):
            if z1 == PLATE_Z and piece in in_plate:
                plate_vols.append(vol)
            elif z1 == COIL_Z1 and piece in in_coil:
                coil_vols.append(vol)
            else:
                air_vols.append(vol)
        bottom = tops
        z0 = z1
    occ.synchronize()

    named = [("Plate", plate_vols), ("Coil", coil_vols), ("Air", air_vols)]
    for name, tags in named:
        if not tags:
            sys.exit(f"no volumes were assigned to block {name}")
        gmsh.model.setPhysicalName(3, gmsh.model.addPhysicalGroup(3, tags), name)
    print("volumes: " + ", ".join(f"{n}={len(t)}" for n, t in named))

    # -- the truncation boundary --------------------------------------------
    # The only named surface: the outside of the free space box, which carries the
    # tangential condition on A, i.e. B.n = 0.
    far = []
    for _, tag in gmsh.model.getEntities(2):
        x0, y0, bz0, x1, y1, bz1 = gmsh.model.getBoundingBox(2, tag)
        if (abs(bz1 - bz0) < TOL and (abs(bz0 - AIR_Z0) < TOL or abs(bz0 - AIR_Z1) < TOL)) or \
           (abs(x1 - x0) < TOL and (abs(x0 - AIR_X0) < TOL or abs(x0 - AIR_X1) < TOL)) or \
           (abs(y1 - y0) < TOL and (abs(y0 - AIR_X0) < TOL or abs(y0 - AIR_X1) < TOL)):
            far.append(tag)
    if not far:
        sys.exit("found no faces on the outer boundary")
    gmsh.model.setPhysicalName(2, gmsh.model.addPhysicalGroup(2, far), "far_field")
    print(f"surfaces: far_field={len(far)}")

    # -- element sizing in plan ---------------------------------------------
    # A box field over the plate and the coil, which covers both measurement lines,
    # grading out to the far field.  Only the plan is sized; the extrusion sets the
    # layers.
    box = gmsh.model.mesh.field.add("Box")
    for k, v in (
        ("XMin", -args.margin), ("XMax", PLATE + args.margin),
        ("YMin", -args.margin), ("YMax", PLATE + args.margin),
        ("ZMin", AIR_Z0 - 1.0), ("ZMax", AIR_Z0 + 1.0),
        ("VIn", args.device_size), ("VOut", args.air_size),
        ("Thickness", args.transition),
    ):
        gmsh.model.mesh.field.setNumber(box, k, v)
    fields = [box]

    # The eddy currents crowd into a band about a skin depth wide along the edges
    # of the plate and of the hole, where they turn to run parallel to them, and
    # three of the measurement points sit on those edges. So the plan is refined
    # towards them, grading back to the device size over edge_band.
    if args.edge_size is not None and args.edge_size < args.device_size:
        def on_edge(bb):
            x0, y0, _, x1, y1, _ = bb
            for lo, hi in ((0.0, PLATE), (HOLE_X0, HOLE_X0 + HOLE)):
                if abs(x1 - x0) < TOL and min(abs(x0 - lo), abs(x0 - hi)) < TOL and \
                   lo - TOL < y0 and y1 < hi + TOL:
                    return True
                if abs(y1 - y0) < TOL and min(abs(y0 - lo), abs(y0 - hi)) < TOL and \
                   lo - TOL < x0 and x1 < hi + TOL:
                    return True
            return False

        edges = [tag for _, tag in gmsh.model.getEntities(1)
                 if abs(gmsh.model.getBoundingBox(1, tag)[2] - AIR_Z0) < TOL
                 and abs(gmsh.model.getBoundingBox(1, tag)[5] - AIR_Z0) < TOL
                 and on_edge(gmsh.model.getBoundingBox(1, tag))]
        if not edges:
            sys.exit("found no plate or hole edges to refine towards")
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "CurvesList", edges)
        gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
        thr = gmsh.model.mesh.field.add("Threshold")
        for k, v in (("InField", dist), ("SizeMin", args.edge_size),
                     ("SizeMax", args.device_size), ("DistMin", 0.0),
                     ("DistMax", args.edge_band), ("StopAtDistMax", 1)):
            gmsh.model.mesh.field.setNumber(thr, k, v)
        fields.append(thr)

    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    for opt, val in (
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeMax", args.air_size),
        ("Mesh.Algorithm", 6),  # Frontal-Delaunay: well shaped triangles
    ):
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ScalingFactor", MM)  # millimetres -> metres
    gmsh.write(args.msh)
    n_prism = len(gmsh.model.mesh.getElementsByType(6)[0])
    gmsh.finalize()
    print(f"wrote {args.msh}: {n_prism} wedges")
    return n_prism


def to_exodus(args):
    """Convert to Exodus so the block and sideset names reach MFEM.

    MFEM picks up named attribute sets from Exodus only; its gmsh reader gives
    bare integer attributes.  MOOSE reads the physical names out of the .msh and
    writes them as element block and sideset names, so a --mesh-only pass is all
    that is needed.
    """
    exe = args.exe or shutil.which("moose_test-opt")
    if not exe:
        moose_dir = os.environ.get("MOOSE_DIR")
        if moose_dir:
            cand = os.path.join(moose_dir, "test", "moose_test-opt")
            exe = cand if os.path.exists(cand) else None
    if not exe:
        sys.exit("no MOOSE executable found; pass --exe or set MOOSE_DIR")

    here = os.path.dirname(os.path.abspath(args.msh)) or "."
    conv = os.path.join(here, ".team7_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by team7_mesh.py; converts the gmsh mesh to Exodus so that the\n"
            "# block and sideset names reach MFEM, which reads named attribute sets from\n"
            "# Exodus only.\n"
            "[Mesh]\n  [fmg]\n    type = FileMeshGenerator\n"
            f"    file = {os.path.basename(args.msh)}\n  []\n[]\n"
            "[Problem]\n  solve = false\n[]\n[Executioner]\n  type = Steady\n[]\n"
        )
    try:
        subprocess.run(
            [exe, "-i", os.path.basename(conv), "--mesh-only", os.path.basename(args.exodus)],
            cwd=here, check=True,
            stdout=None if args.verbose else subprocess.DEVNULL,
        )
    finally:
        os.path.exists(conv) and os.remove(conv)
    print(f"wrote {args.exodus}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--device-size", type=float, default=12.0,
                   help="element size in plan over the plate and the coil (mm); the "
                        "measurements are 18 mm apart")
    p.add_argument("--edge-size", type=float, default=None,
                   help="refine the plan towards the edges of the plate and of the hole, to "
                        "this element size (mm). Off by default: refining them to 6 mm "
                        "moves the values sampled right above or on those edges by up to 8%% "
                        "of the peak, but leaves the overall agreement with the measurements "
                        "unchanged, and nearly doubles the cost")
    p.add_argument("--edge-band", type=float, default=30.0,
                   help="distance from those edges over which the size grades back to the "
                        "device size (mm)")
    p.add_argument("--margin", type=float, default=20.0,
                   help="how far the fine region extends beyond the plate in plan (mm)")
    p.add_argument("--transition", type=float, default=150.0,
                   help="distance in plan over which the size grades out to the far field (mm)")
    p.add_argument("--air-size", type=float, default=500.0,
                   help="element size in plan at the outer boundary (mm)")
    p.add_argument("--plate-layers", type=int, default=6,
                   help="layers through the 19 mm plate")
    p.add_argument("--gap-layers", type=int, default=2,
                   help="layers in each of the two 15 mm slabs between the plate and the coil")
    p.add_argument("--coil-layers", type=int, default=4,
                   help="layers through the 100 mm height of the coil")
    p.add_argument("--far-layers", type=int, default=4,
                   help="layers in each of the free space slabs above and below the device")
    p.add_argument("--grading", type=float, default=1.6,
                   help="growth ratio of the free space layers away from the device")
    p.add_argument("--msh", default="team7.msh")
    p.add_argument("--exodus", default="team7.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
