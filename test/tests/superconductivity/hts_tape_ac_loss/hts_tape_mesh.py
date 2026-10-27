#!/usr/bin/env python3
"""Build the mesh for the HTS tape AC loss benchmark with gmsh.

The benchmark is the first of the HTS Modelling Workgroup: the cross-section of a
single superconducting tape, 12 mm wide and 1 um thick, carrying an alternating
transport current.  It is documented in doc/content/problems/hts_tape_ac_loss.md.

Only the quarter of the cross-section with x >= 0 and y >= 0 is meshed: the tape
is symmetric about both axes, and so is the current it carries.  The free space
around it is a quarter disc, on whose arc the field of the transport current is
imposed.

The tape is a block of quadrilaterals, structured so that its 1 um thickness costs
no more elements than are asked for across it.  They are graded towards the edge
of the tape, where the current first penetrates and where the field is singular.
The free space is meshed with triangles, which grow away from the tip of the tape
at a fixed rate from its 0.5 um half-thickness to the size asked for at the arc.

Dimensions are in millimetres here; the mesh is written in metres.  Run with
--help for the meshing controls.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys

TAPE_HALF_WIDTH = 6.0  # half of the 12 mm width
TAPE_HALF_THICKNESS = 0.5e-3  # half of the 1 um thickness

MM = 1.0e-3


def build(args):
    import gmsh

    a, b, r = TAPE_HALF_WIDTH, TAPE_HALF_THICKNESS, args.air_radius

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("hts_tape")
    geo = gmsh.model.geo

    origin = geo.addPoint(0, 0, 0)
    tape_edge = geo.addPoint(a, 0, 0)
    tape_tip = geo.addPoint(a, b, 0)
    tape_top = geo.addPoint(0, b, 0)
    arc_x = geo.addPoint(r, 0, 0)
    arc_y = geo.addPoint(0, r, 0)

    tape_bottom = geo.addLine(origin, tape_edge)
    tape_end = geo.addLine(tape_edge, tape_tip)
    tape_face = geo.addLine(tape_tip, tape_top)
    tape_axis = geo.addLine(tape_top, origin)
    air_bottom = geo.addLine(tape_edge, arc_x)
    arc = geo.addCircleArc(arc_x, origin, arc_y)
    air_axis = geo.addLine(arc_y, tape_top)

    tape = geo.addPlaneSurface([geo.addCurveLoop([tape_bottom, tape_end, tape_face, tape_axis])])
    air = geo.addPlaneSurface(
        [geo.addCurveLoop([air_bottom, arc, air_axis, -tape_face, -tape_end])])

    # Element lengths across the width shrink geometrically towards the edge, by
    # edge_grading overall. gmsh's progression runs in the direction of the curve,
    # and the two long sides of the tape run in opposite directions.
    n = args.width_elements
    ratio = args.edge_grading ** (-1.0 / max(n - 1, 1))
    geo.mesh.setTransfiniteCurve(tape_bottom, n + 1, "Progression", ratio)
    geo.mesh.setTransfiniteCurve(tape_face, n + 1, "Progression", 1.0 / ratio)
    for c in (tape_end, tape_axis):
        geo.mesh.setTransfiniteCurve(c, args.thickness_elements + 1)
    # The transport current enters through the tangential field on the arc, and each
    # edge of it gets the field at its midpoint: the current comes out too high by a
    # fraction (edge angle)^2 / 12.
    geo.mesh.setTransfiniteCurve(arc, math.ceil(0.5 * math.pi * r / args.arc_size) + 1)
    geo.mesh.setTransfiniteSurface(tape)
    geo.mesh.setRecombine(2, tape)

    geo.synchronize()

    for dim, tags, name in (
        (2, [tape], "Tape"),
        (2, [air], "Air"),
        (1, [arc], "far_field"),
        (1, [tape_bottom, air_bottom], "x_axis"),
        (1, [tape_axis, air_axis], "y_axis"),
    ):
        gmsh.model.setPhysicalName(dim, gmsh.model.addPhysicalGroup(dim, tags), name)

    # The size grows linearly with distance from the tip of the tape, from the
    # spacing the tape's end is cut into there. The tape's own boundary nodes set
    # the size along its face, which the free space takes up from there.
    tip_size = b / args.thickness_elements
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "PointsList", [tape_tip])
    grow = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        grow, "F", f"Min({args.air_size}, {tip_size} + {args.growth} * F{dist})")
    gmsh.model.mesh.field.setAsBackgroundMesh(grow)
    for opt, val in (
        ("Mesh.MeshSizeExtendFromBoundary", 1),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeMax", args.air_size),
        ("Mesh.Algorithm", 6),  # frontal-Delaunay
    ):
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ScalingFactor", MM)  # millimetres -> metres
    gmsh.write(args.msh)
    n_quad = len(gmsh.model.mesh.getElementsByType(3)[0])
    n_tri = len(gmsh.model.mesh.getElementsByType(2)[0])
    gmsh.finalize()
    print(f"wrote {args.msh}: {n_quad} quadrilaterals in the tape, {n_tri} triangles")


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
    conv = os.path.join(here, ".hts_tape_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by hts_tape_mesh.py; converts the gmsh mesh to Exodus so that the\n"
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
    p.add_argument("--width-elements", type=int, default=100,
                   help="elements across the half-width of the tape")
    p.add_argument("--thickness-elements", type=int, default=1,
                   help="elements across the half-thickness of the tape")
    p.add_argument("--edge-grading", type=float, default=20.0,
                   help="ratio of the widest element in the tape, at its centre, to "
                        "the narrowest, at its edge")
    p.add_argument("--growth", type=float, default=0.3,
                   help="rate at which the element size in free space grows with "
                        "distance from the tip of the tape")
    p.add_argument("--air-size", type=float, default=4.0,
                   help="largest element size in free space (mm)")
    p.add_argument("--arc-size", type=float, default=1.0,
                   help="element size along the arc bounding free space (mm)")
    p.add_argument("--air-radius", type=float, default=60.0,
                   help="radius of the free space around the tape (mm)")
    p.add_argument("--msh", default="hts_tape.msh")
    p.add_argument("--exodus", default="hts_tape.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
