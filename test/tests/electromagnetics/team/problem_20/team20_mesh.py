#!/usr/bin/env python3
"""Build the TEAM Problem 20 mesh with gmsh.

TEAM Problem 20 is a 3-D non-linear magnetostatic force benchmark: a coil wound
round a steel centre pole drives flux up the pole, across into the top plates of
a steel yoke, down its two legs and back across a 1.5 mm gap into the foot of the
pole.  The quantity asked for is the force on the pole.  The benchmark and the
geometry are documented in doc/content/problems/team_problem_20.md.

Only the quarter of the model with x >= 0 and y >= 0 is meshed, which is the
region Fig. 2 of the problem statement allows.  The planes x = 0 and y = 0 both
contain the axis of the coil, and the field lies in them, so they carry the same
tangential condition on A as the outer boundary.

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
POLE_HX = 12.5  # half of the 25 mm pole width in x
POLE_HY = 5.0  # half of the 10 mm pole depth in y
POLE_Z0 = 26.5  # foot of the pole, 1.5 mm above the yoke
POLE_Z1 = 125.0  # top of the pole, level with the underside of the top plates

YOKE_HX = 63.5  # outer face of the legs: 50 + 27 / 2
YOKE_HY = 12.5  # half of the 25 mm yoke depth in y
YOKE_LEG_X0 = YOKE_HX - 25.0  # 38.5, inner face of the legs
YOKE_OPEN_HX = 13.5  # half of the 27 mm opening between the top plates
YOKE_Z0, YOKE_Z1 = 25.0, 125.0  # the window, between the bottom and top plates
YOKE_TOP = 150.0

COIL_A = 14.5  # corner arcs are centred at (+-14.5, +-14.5)
COIL_RI, COIL_RO = 5.0, 23.0  # inner and outer corner radii
COIL_Z0, COIL_Z1 = YOKE_Z0 + 1.7, YOKE_Z1 - 1.7  # 1.7 mm clear of the yoke

# Table 3 of the problem statement asks for the average Bz over the section of the
# pole and of the yoke leg at z = 75 mm.  Those are sections rather than volumes,
# so each is reported as the average over a band of steel centred on it, carved
# out of the part.  The field in both is close to uniform along z there, so the
# band stands in for the section closely.
MEASURE_Z = 75.0
MEASURE_BAND = 5.0

MM = 1.0e-3

# Comparison tolerance in millimetres, used for locating entities by position.
TOL = 1.0e-4


def build(args):
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("team20")
    occ = gmsh.model.occ

    air_xy = args.air_xy
    air_z0, air_z1 = -args.air_below, YOKE_TOP + args.air_above
    z_lo, z_hi = MEASURE_Z - MEASURE_BAND / 2, MEASURE_Z + MEASURE_BAND / 2

    # -- the steel ----------------------------------------------------------
    pole = occ.addBox(0.0, 0.0, POLE_Z0, POLE_HX, POLE_HY, POLE_Z1 - POLE_Z0)
    # The yoke is a rectangular frame in the x-z plane with a slot cut through
    # the top: the bottom plate, a leg, and a top plate stopping short of the pole.
    yoke = [
        occ.addBox(0.0, 0.0, 0.0, YOKE_HX, YOKE_HY, YOKE_Z0),
        occ.addBox(YOKE_LEG_X0, 0.0, YOKE_Z0, YOKE_HX - YOKE_LEG_X0, YOKE_HY, YOKE_Z1 - YOKE_Z0),
        occ.addBox(YOKE_OPEN_HX, 0.0, YOKE_Z1, YOKE_HX - YOKE_OPEN_HX, YOKE_HY, YOKE_TOP - YOKE_Z1),
    ]
    # Carved back out of the pole and the leg below; each lies wholly in the steel.
    pole_band = occ.addBox(0.0, 0.0, z_lo, POLE_HX, POLE_HY, z_hi - z_lo)
    yoke_band = occ.addBox(YOKE_LEG_X0, 0.0, z_lo, YOKE_HX - YOKE_LEG_X0, YOKE_HY, z_hi - z_lo)

    # -- the coil -----------------------------------------------------------
    # A rounded square ring: a filleted box with a smaller filleted box cut out
    # of it.  Filleting the four vertical edges of a box of half-width a with
    # radius r gives a rounded square whose corner arcs are centred at
    # (+-(a - r), +-(a - r)), so both boxes share the arc centres the problem
    # statement draws.  The ring is then cut down to the quarter being meshed.
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
            # feature in the model (the smallest is the 1 mm clearance at the top
            # of the pole).
            if max(x1 - x0, y1 - y0) < TOL and abs(zc1 - zc0 - dz) < TOL:
                vertical.append(c)
        return occ.fillet([tag], vertical, [radius])[0][1]

    coil_h = COIL_Z1 - COIL_Z0
    outer = rounded_box(COIL_A + COIL_RO, COIL_RO, COIL_Z0, coil_h)
    inner = rounded_box(COIL_A + COIL_RI, COIL_RI, COIL_Z0 - 5.0, coil_h + 10.0)
    ring = occ.cut([(3, outer)], [(3, inner)])[0][0][1]
    quarter = occ.addBox(0.0, 0.0, COIL_Z0 - 5.0, 100.0, 100.0, coil_h + 10.0)
    coil = occ.intersect([(3, ring)], [(3, quarter)])[0][0][1]

    # -- free space and the conformal assembly ------------------------------
    box = occ.addBox(0.0, 0.0, air_z0, air_xy, air_xy, air_z1 - air_z0)
    solids = [pole] + yoke + [pole_band, yoke_band, coil]
    _, fmap = occ.fragment([(3, box)], [(3, t) for t in solids])
    occ.synchronize()

    # -- name the volumes ---------------------------------------------------
    # fragment reports, for each input solid, the pieces it ended up as, so the
    # blocks are read straight off that rather than guessed at from positions.
    # The first entry is the free space box, which fragment splits into every
    # piece of the model; whatever is not claimed by a solid is the air.
    pieces = [set(t for _, t in ents) for ents in fmap]
    pole_tags = pieces[1]
    yoke_tags = set().union(*pieces[2:5])
    pole_band_tags, yoke_band_tags, coil_tags = pieces[5], pieces[6], pieces[7]
    air_tags = pieces[0] - pole_tags - yoke_tags - coil_tags
    # The measurement bands overlap the pole and the leg, so take them back out.
    pole_tags -= pole_band_tags
    yoke_tags -= yoke_band_tags

    named = [
        ("Pole", pole_tags),
        ("PoleZ075", pole_band_tags),
        ("Yoke", yoke_tags),
        ("YokeZ075", yoke_band_tags),
        ("Coil", coil_tags),
        ("Air", air_tags),
    ]
    for name, tags in named:
        if not tags:
            sys.exit(f"no volumes were assigned to block {name}")
        gmsh.model.setPhysicalName(3, gmsh.model.addPhysicalGroup(3, sorted(tags)), name)
    print("volumes: " + ", ".join(f"{n}={len(t)}" for n, t in named))

    # -- name the surfaces --------------------------------------------------
    # The two symmetry planes and the truncation boundary.  All three carry the
    # same condition in the solve, but they are named apart so that the
    # documentation and any later change of boundary condition can tell them
    # apart.  The tests are on bounding boxes rather than centres of mass so
    # that a large face cannot be mistaken for a small one it is centred on.
    surf = {}

    def add(name, tag):
        surf.setdefault(name, []).append(tag)

    for _, tag in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
        if abs(x1 - x0) < TOL and abs(x0) < TOL:
            add("symmetry_x", tag)
        elif abs(y1 - y0) < TOL and abs(y0) < TOL:
            add("symmetry_y", tag)
        elif (abs(x1 - x0) < TOL and abs(x0 - air_xy) < TOL) or \
             (abs(y1 - y0) < TOL and abs(y0 - air_xy) < TOL) or \
             (abs(z1 - z0) < TOL and (abs(z0 - air_z0) < TOL or abs(z0 - air_z1) < TOL)):
            add("far_field", tag)

    for name, tags in sorted(surf.items()):
        gmsh.model.setPhysicalName(2, gmsh.model.addPhysicalGroup(2, tags), name)
    print("surfaces: " + ", ".join(f"{k}={len(v)}" for k, v in sorted(surf.items())))
    missing = [n for n in ("far_field", "symmetry_x", "symmetry_y") if n not in surf]
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
        # the 1.5 mm gap under the foot of the pole
        box_field(((0.0, 0.0, YOKE_Z0 - 1.0), (POLE_HX + 2.0, POLE_HY + 2.0, POLE_Z0 + 1.0)),
                  args.gap_size, 6.0),
        # the top of the pole, where the flux crosses the 1 mm to the top plates
        box_field(((POLE_HX - 2.0, 0.0, POLE_Z1 - 3.0),
                   (YOKE_OPEN_HX + 2.0, POLE_HY + 2.0, POLE_Z1 + 3.0)),
                  args.gap_size, 6.0),
        # the pole and the air immediately round it, which the force is taken over
        box_field(((0.0, 0.0, POLE_Z0 - 1.0), (POLE_HX + 1.0, POLE_HY + 1.0, POLE_Z1 + 1.0)),
                  args.pole_size, 10.0),
        # the rest of the steel
        box_field(((0.0, 0.0, 0.0), (YOKE_HX, YOKE_HY, YOKE_TOP)), args.steel_size, 20.0),
        # the coil and the air immediately around the device
        box_field(((0.0, 0.0, 0.0), (COIL_A + COIL_RO + 5.0, COIL_A + COIL_RO + 5.0, YOKE_TOP)),
                  args.coil_size, 40.0),
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
    conv = os.path.join(here, ".team20_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by team20_mesh.py; converts the gmsh mesh to Exodus so that the\n"
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
    p.add_argument("--gap-size", type=float, default=0.3,
                   help="element size in the gap under the pole and at its top corner (mm)")
    p.add_argument("--pole-size", type=float, default=1.0,
                   help="element size in the pole and the air round it (mm)")
    p.add_argument("--steel-size", type=float, default=3.0,
                   help="element size in the yoke (mm)")
    p.add_argument("--coil-size", type=float, default=6.0,
                   help="element size in the coil and the air around the device (mm)")
    p.add_argument("--air-size", type=float, default=40.0,
                   help="element size at the outer boundary (mm)")
    p.add_argument("--air-xy", type=float, default=250.0,
                   help="extent of the free space box in x and y (mm)")
    p.add_argument("--air-below", type=float, default=200.0,
                   help="free space below the yoke (mm)")
    p.add_argument("--air-above", type=float, default=200.0,
                   help="free space above the yoke (mm)")
    p.add_argument("--msh", default="team20.msh")
    p.add_argument("--exodus", default="team20.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
