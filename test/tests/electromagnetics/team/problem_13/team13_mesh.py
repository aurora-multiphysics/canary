#!/usr/bin/env python3
"""Build the TEAM Problem 13 mesh with gmsh.

TEAM Problem 13 is a 3-D non-linear magnetostatic benchmark: a racetrack coil
drives flux up a central steel plate, across a 0.5 mm air gap, and back down two
steel channels whose openings face the plate.  The benchmark and the geometry are
documented in doc/content/problems/team_problem_13.md.

Only the half of the model above z = 0 is meshed.  The geometry and the coil
current are symmetric under reflection in that plane, which makes B normal to it,
so the plane carries the natural boundary condition of the A-formulation and no
constraint has to be applied there.

Dimensions come from Fig. 1 of the problem statement and are in millimetres here;
the mesh is written in metres.  Run with --help for the meshing controls.
"""

import argparse
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# Geometry, in millimetres, straight from Fig. 1 of the problem statement.
# ---------------------------------------------------------------------------
PLATE_HX = 1.6  # half of the 3.2 mm plate thickness
PLATE_HY = 25.0  # half of the 50 mm plate width
PLATE_Z = 63.2  # plate top, level with the top of the channels

GAP = 0.5  # plate face to channel flange end

CH_X0 = PLATE_HX + GAP  # 2.1, inner end of the flanges
CH_X1 = 125.3  # outer face of the web
CH_WEB_X0 = CH_X1 - 3.2  # 122.1, inner face of the web
CH_Y0, CH_Y1 = 15.0, 65.0  # the channels are offset in y
CH_Z0 = 60.0  # underside of the top flange
CH_Z1 = PLATE_Z  # 63.2

COIL_A = 50.0  # centre of the corner arcs, at (+-50, +-50)
COIL_RI, COIL_RO = 25.0, 50.0  # inner and outer corner radii
COIL_HALF = 75.0  # inner opening half-width (150 mm across)
COIL_Z = 50.0  # coil half-height (100 mm tall)

# Positions 1 to 7 of Table 1 of the problem statement: the heights up the centre
# plate at which the average flux density in the steel is reported.
MEASURE_Z = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]

# The plate is built as a stack of slabs, one block per measurement height, so
# that the average over each can be reported.  These are the cuts between them,
# placed so that each slab is centred on its measurement height.  The lowest
# slab is only the top half of its band, the rest of it being below the symmetry
# plane, and the highest runs to the top of the plate rather than to 65 mm.
PLATE_CUTS = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0]

# Positions 8 to 18 of Table 1: stations along the top flange of the +x channel,
# each reported over the full width 15 <= y <= 65 and thickness 60 <= z <= 63.2.
FLANGE_X = [2.1, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 110.0, 122.1]

# Positions 19 to 25: stations down the web, over 122.1 <= x <= 125.3 and the same
# width in y.
WEB_Z = [60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 0.0]

# Those are sections rather than volumes, so each is reported as the average over
# a band of steel centred on it, carved out of the channel.  The band is narrow
# because the field changes quickly along the flange near the plate.  Where a
# station sits at the end of its part the band is clipped and so one-sided: that
# is x = 2.1 and x = 122.1 on the flange and z = 60 on the web.  At z = 0 the
# clipped half-band is still centred, the rest of it lying below the symmetry
# plane.  The other channel is the mirror of this one and carries no blocks.
MEASURE_BAND = 5.0

MM = 1.0e-3

# Comparison tolerance in millimetres, used for locating entities by position.
TOL = 1.0e-4


def build(args):
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("team13")
    occ = gmsh.model.occ

    air = args.air_xy
    air_z = args.air_z

    # -- the steel ----------------------------------------------------------
    # The plate is built as a stack of slabs, one per measurement height.
    zs = [0.0] + PLATE_CUTS + [PLATE_Z]
    plate = [
        occ.addBox(-PLATE_HX, -PLATE_HY, z0, 2 * PLATE_HX, 2 * PLATE_HY, z1 - z0)
        for z0, z1 in zip(zs[:-1], zs[1:])
    ]

    # Each channel is a C-section: a top flange reaching in towards the plate and
    # a web closing it off at the outer end.  The second channel is the first
    # rotated by 180 degrees about the z axis, which is what makes the two
    # channels sit on opposite sides in y.
    channel = [
        occ.addBox(CH_X0, CH_Y0, CH_Z0, CH_X1 - CH_X0, CH_Y1 - CH_Y0, CH_Z1 - CH_Z0),
        occ.addBox(CH_WEB_X0, CH_Y0, 0.0, CH_X1 - CH_WEB_X0, CH_Y1 - CH_Y0, CH_Z1),
        occ.addBox(-CH_X1, -CH_Y1, CH_Z0, CH_X1 - CH_X0, CH_Y1 - CH_Y0, CH_Z1 - CH_Z0),
        occ.addBox(-CH_X1, -CH_Y1, 0.0, CH_X1 - CH_WEB_X0, CH_Y1 - CH_Y0, CH_Z1),
    ]
    # Carved back out of the channel below; each lies wholly inside the steel.
    def band(lo, hi, centre):
        return max(lo, centre - MEASURE_BAND / 2), min(hi, centre + MEASURE_BAND / 2)

    measure, measure_names = [], []
    for x in FLANGE_X:
        x0, x1 = band(CH_X0, CH_WEB_X0, x)
        measure.append(
            occ.addBox(x0, CH_Y0, CH_Z0, x1 - x0, CH_Y1 - CH_Y0, CH_Z1 - CH_Z0))
        measure_names.append("FlangeX%03d" % round(x))
    for z in WEB_Z:
        z0, z1 = band(0.0, CH_Z0, z)
        measure.append(
            occ.addBox(CH_WEB_X0, CH_Y0, z0, CH_X1 - CH_WEB_X0, CH_Y1 - CH_Y0, z1 - z0))
        measure_names.append("WebZ%03d" % round(z))

    # -- the coil -----------------------------------------------------------
    # A rounded rectangular ring: a filleted box with a smaller filleted box cut
    # out of it.  Filleting the four vertical edges of a box of half-width a with
    # radius r gives exactly the rounded rectangle of the problem statement,
    # whose corner arcs are centred at (+-(a - r), +-(a - r)).
    def rounded_box(half, radius, z0, dz):
        tag = occ.addBox(-half, -half, z0, 2 * half, 2 * half, dz)
        occ.synchronize()
        curves = set()
        for face in gmsh.model.getBoundary([(3, tag)], combined=False, oriented=False):
            for c in gmsh.model.getBoundary([face], combined=False, oriented=False):
                curves.add(abs(c[1]))
        vertical = []
        for c in curves:
            x0, y0, zc0, x1, y1, zc1 = gmsh.model.getBoundingBox(1, c)
            # OCC pads bounding boxes by its own gap tolerance, so this test has
            # to be loose compared with 1e-9 but is still far tighter than any
            # feature in the model (the smallest is the 0.5 mm gap).
            if max(x1 - x0, y1 - y0) < TOL and abs(zc1 - zc0 - dz) < TOL:
                vertical.append(c)
        return occ.fillet([tag], vertical, [radius])[0][1]

    outer = rounded_box(COIL_A + COIL_RO, COIL_RO, 0.0, COIL_Z)
    inner = rounded_box(COIL_A + COIL_RI, COIL_RI, -5.0, COIL_Z + 10.0)
    coil = occ.cut([(3, outer)], [(3, inner)])[0][0][1]

    # The coil is a closed loop.  The alternative to the analytic winding direction
    # finds it by solving for a potential that jumps by a fixed amount across a cut
    # through the loop (see team13_coil_winding_laplace.i); the analytic one does
    # not need the cut.  The cut is the coil's cross-section in the plane x = 0 on
    # the +y side, where the coil is straight.  It is meshed as an interior face,
    # and named as a sideset in to_exodus().
    section = [(COIL_HALF, 0.0), (COIL_A + COIL_RO, 0.0),
               (COIL_A + COIL_RO, COIL_Z), (COIL_HALF, COIL_Z)]
    points = [occ.addPoint(0.0, y, z) for y, z in section]
    lines = [occ.addLine(points[i], points[(i + 1) % 4]) for i in range(4)]
    cut = occ.addPlaneSurface([occ.addCurveLoop(lines)])

    # -- free space and the conformal assembly ------------------------------
    box = occ.addBox(-air, -air, 0.0, 2 * air, 2 * air, air_z)
    solids = [(3, t) for t in plate + channel + measure + [coil]]
    _, fmap = occ.fragment([(3, box)], solids + [(2, cut)])
    occ.synchronize()

    # -- name the volumes ---------------------------------------------------
    # fragment reports, for each input solid, the pieces it ended up as, so the
    # blocks are read straight off that rather than guessed at from positions.
    # The first entry is the free space box, which fragment splits into every
    # piece of the model; whatever is not claimed by a solid is the air. The last
    # entry is the cut, which is a surface rather than a volume.
    if not fmap[-1]:
        sys.exit("the cut through the coil was lost in the fragment")
    pieces = [set(t for _, t in ents) for ents in fmap[:-1]]
    n_p, n_c, n_m = len(plate), len(channel), len(measure)
    plate_tags = set().union(*pieces[1:1 + n_p])
    channel_tags = set().union(*pieces[1 + n_p:1 + n_p + n_c])
    measure_tags = [pieces[1 + n_p + n_c + i] for i in range(n_m)]
    coil_tags = pieces[1 + n_p + n_c + n_m]
    air_tags = pieces[0] - plate_tags - channel_tags - coil_tags
    # The measurement boxes overlap the channel boxes, so take them back out.
    channel_tags -= set().union(*measure_tags)

    slab_tags = [sorted(pieces[1 + i]) for i in range(len(plate))]
    named = [(f"PlateZ{int(MEASURE_Z[i]):03d}", t) for i, t in enumerate(slab_tags)]
    named += [(n, sorted(t)) for n, t in zip(measure_names, measure_tags)]
    named += [("Channel", sorted(channel_tags)), ("Coil", sorted(coil_tags)),
              ("Air", sorted(air_tags))]
    for name, tags in named:
        if not tags:
            sys.exit(f"no volumes were assigned to block {name}")
        gmsh.model.setPhysicalName(3, gmsh.model.addPhysicalGroup(3, tags), name)
    print("volumes: " + ", ".join(f"{n}={len(t)}" for n, t in named))

    # -- name the surfaces --------------------------------------------------
    # Only the two faces the solve needs: the symmetry plane and the truncation
    # boundary that carries the tangential condition on A.  The cut through the
    # coil is named later, in to_exodus().  The measured sections
    # are reported as averages over the blocks above rather than as fluxes over
    # sidesets, so no section surfaces are named.  The tests are on bounding boxes
    # rather than centres of mass so that a large face cannot be mistaken for a
    # small one it happens to be centred on.
    surf = {}

    def add(name, tag):
        surf.setdefault(name, []).append(tag)

    def within(bb, lo, hi):
        return all(lo[i] - TOL < bb[i] and bb[i + 3] < hi[i] + TOL for i in range(3))

    BIG = 1e9
    for _, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.getBoundingBox(2, tag)
        flat_z = abs(bb[5] - bb[2]) < TOL
        if flat_z and abs(bb[2]) < TOL:
            add("symmetry", tag)
        elif (abs(abs(bb[0]) - air) < TOL and abs(bb[3] - bb[0]) < TOL) or \
             (abs(abs(bb[1]) - air) < TOL and abs(bb[4] - bb[1]) < TOL) or \
             (flat_z and abs(bb[2] - air_z) < TOL):
            add("far_field", tag)

    for name, tags in sorted(surf.items()):
        gmsh.model.setPhysicalName(2, gmsh.model.addPhysicalGroup(2, tags), name)
    print("surfaces: " + ", ".join(f"{k}={len(v)}" for k, v in sorted(surf.items())))
    missing = [n for n in ("far_field", "symmetry") if n not in surf]
    if missing:
        sys.exit("missing expected surface groups: " + ", ".join(missing))

    # -- element sizing -----------------------------------------------------
    # Boxes rather than a distance field: the features that set the resolution
    # are all axis aligned, and a box field with a transition thickness grades
    # out of them without needing the geometry to be interrogated.
    def box_field(bounds, size, thickness):
        f = gmsh.model.mesh.field.add("Box")
        (x0, y0, z0), (x1, y1, z1) = bounds
        for k, v in (
            ("XMin", x0), ("XMax", x1), ("YMin", y0), ("YMax", y1),
            ("ZMin", z0), ("ZMax", z1), ("VIn", size), ("VOut", 1e22),
            ("Thickness", thickness),
        ):
            gmsh.model.mesh.field.setNumber(f, k, v)
        return f

    fields = [
        # the two 0.5 mm gaps, which are the tightest feature in the model
        box_field(((1.2, 12.0, 56.0), (2.7, 28.0, 64.0)), args.gap_size, 8.0),
        box_field(((-2.7, -28.0, 56.0), (-1.2, -12.0, 64.0)), args.gap_size, 8.0),
        # the steel, which has to resolve the 3.2 mm thicknesses
        box_field(((-CH_X1 - 1, -CH_Y1 - 1, -1.0), (CH_X1 + 1, CH_Y1 + 1, PLATE_Z + 1)),
                  args.steel_size, 30.0),
        # the coil and the air immediately around the device
        box_field(((-110.0, -110.0, -1.0), (110.0, 110.0, 60.0)), args.coil_size, 60.0),
    ]
    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    for opt, val in (
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeMax", args.air_size),
        ("Mesh.Algorithm3D", 10),  # HXT: fast, and robust on this geometry
        ("Mesh.Optimize", 1),
        ("Mesh.OptimizeNetgen", 0),
    ):
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ScalingFactor", MM)  # millimetres -> metres
    gmsh.write(args.msh)
    n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
    gmsh.finalize()
    print(f"wrote {args.msh}: {n_tet} tetrahedra")
    return n_tet


def to_exodus(args):
    """Convert to Exodus so the block and sideset names reach MFEM.

    MFEM picks up named attribute sets from Exodus only; its gmsh reader gives
    bare integer attributes.  MOOSE reads the physical names out of the .msh and
    writes them as element block and sideset names, so a --mesh-only pass is all
    that is needed.

    The cut through the coil is named here rather than in gmsh.  A physical
    surface on an interior face reaches the Exodus file once from each element
    either side of it, and MFEM rejects a boundary with every face doubled.
    Taking only the sides whose outward normal is +x keeps one of each.
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
    conv = os.path.join(here, ".team13_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by team13_mesh.py; converts the gmsh mesh to Exodus so that the\n"
            "# block and sideset names reach MFEM, which reads named attribute sets from\n"
            "# Exodus only.\n"
            "[Mesh]\n  [fmg]\n    type = FileMeshGenerator\n"
            f"    file = {os.path.basename(args.msh)}\n  []\n"
            "  [coil_cut]\n    type = ParsedGenerateSideset\n    input = fmg\n"
            "    combinatorial_geometry = 'abs(x) < 1e-7 & y > 0'\n"
            "    included_subdomains = Coil\n    included_neighbors = Coil\n"
            "    normal = '1 0 0'\n    new_sideset_name = coil_cut\n  []\n[]\n"
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
    p.add_argument("--gap-size", type=float, default=0.25,
                   help="element size in the 0.5 mm plate-to-channel gaps (mm); the "
                        "setting the answer is by far the most sensitive to")
    p.add_argument("--steel-size", type=float, default=4.0,
                   help="element size in the steel (mm); the parts are 3.2 mm thick")
    p.add_argument("--coil-size", type=float, default=14.0,
                   help="element size in the coil and the air around the device (mm)")
    p.add_argument("--air-size", type=float, default=40.0,
                   help="element size at the outer boundary (mm)")
    p.add_argument("--air-xy", type=float, default=250.0,
                   help="half-width of the free space box in x and y (mm)")
    p.add_argument("--air-z", type=float, default=200.0,
                   help="height of the free space box above the symmetry plane (mm)")
    p.add_argument("--msh", default="team13.msh")
    p.add_argument("--exodus", default="team13.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
