#!/usr/bin/env python3
"""Build the mesh for the electrical contact example with gmsh.

Two coaxial rods, copper above and aluminium below, are pressed together end to
end and touch over a circular spot on their axis. The example is documented in
doc/content/problems/electrical_contact.md. The problem is axisymmetric, so the
mesh is of the half-plane r = x >= 0, with z = y along the axis.

The spot is covered by a film, which is meshed as a layer of the given thickness:
a single row of quadrilaterals, graded towards the edge of the spot, where the
current crowds and the field is singular. Outside the spot the rods are separated
by a gap of the same thickness, which is not meshed: their faces either side of it
are insulated. The rods themselves are meshed with triangles, which grow away from
the edge of the spot at a fixed rate from the size of the film's first element to
the largest size asked for.

With --spot-radius equal to --rod-radius the rods touch over their whole faces
and the problem is one-dimensional.

Dimensions are in millimetres here; the mesh is written in metres. Run with --help
for the meshing controls.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys

MM = 1.0e-3


def build(args):
    import gmsh

    a, b, half_length = args.spot_radius, args.rod_radius, args.rod_half_length
    h = 0.5 * args.film_thickness
    if not 0 < a <= b:
        sys.exit("the spot radius must be positive and no larger than the rod radius")
    full_contact = math.isclose(a, b)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("electrical_contact")
    geo = gmsh.model.geo

    # Points on the film, and on the faces of the rods either side of the gap.
    film_axis = [geo.addPoint(0, z, 0) for z in (-h, h)]
    spot_edge = [geo.addPoint(a, z, 0) for z in (-h, h)]
    rim = spot_edge if full_contact else [geo.addPoint(b, z, 0) for z in (-h, h)]
    electrode_axis = [geo.addPoint(0, z, 0) for z in (-half_length, half_length)]
    electrode_rim = [geo.addPoint(b, z, 0) for z in (-half_length, half_length)]

    # The film. Its faces run from the axis to the edge of the spot.
    film_face = [geo.addLine(film_axis[i], spot_edge[i]) for i in (0, 1)]
    film_axis_line = geo.addLine(film_axis[0], film_axis[1])
    film_edge_line = geo.addLine(spot_edge[0], spot_edge[1])
    film = geo.addPlaneSurface([geo.addCurveLoop(
        [film_face[0], film_edge_line, -film_face[1], -film_axis_line])])

    # The rods: i = 0 is the aluminium below the film, i = 1 the copper above it.
    rods, gap_faces, electrodes, sides, axes = [], [], [], [], []
    for i, sign in ((0, -1), (1, 1)):
        gap = None if full_contact else geo.addLine(spot_edge[i], rim[i])
        side = geo.addLine(rim[i], electrode_rim[i])
        electrode = geo.addLine(electrode_rim[i], electrode_axis[i])
        axis = geo.addLine(electrode_axis[i], film_axis[i])
        loop = [film_face[i]] + ([gap] if gap else []) + [side, electrode, axis]
        if sign < 0:  # keep the loop anticlockwise below the film
            loop = [-c for c in reversed(loop)]
        rods.append(geo.addPlaneSurface([geo.addCurveLoop(loop)]))
        gap_faces.append(gap)
        sides.append(side)
        electrodes.append(electrode)
        axes.append(axis)

    # Element lengths along the film's faces shrink geometrically towards the edge
    # of the spot, from the largest at the axis to the first, film_size, at the edge,
    # growing by 1 + growth from one element to the next.
    edge_size = args.edge_size * args.film_thickness
    ratio = 1 + args.growth
    n = max(1, math.ceil(math.log(1 + args.growth * a / edge_size) / math.log(ratio)))
    for c in film_face:
        geo.mesh.setTransfiniteCurve(c, n + 1, "Progression", 1 / ratio)
    for c in (film_axis_line, film_edge_line):
        geo.mesh.setTransfiniteCurve(c, args.film_elements + 1)
    geo.mesh.setTransfiniteSurface(film)
    geo.mesh.setRecombine(2, film)
    geo.synchronize()

    groups = [
        (2, [rods[1]], "Copper"),
        (2, [rods[0]], "Aluminium"),
        (2, [film], "Film"),
        (1, [electrodes[1]], "copper_electrode"),
        (1, [electrodes[0]], "aluminium_electrode"),
        (1, axes + [film_axis_line], "axis"),
        (1, sides, "outer"),
    ]
    if not full_contact:
        groups.append((1, [film_edge_line] + gap_faces, "gap"))
    for dim, tags, name in groups:
        gmsh.model.setPhysicalName(dim, gmsh.model.addPhysicalGroup(dim, tags), name)

    # In the rods the size grows linearly with distance from the edge of the spot.
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "PointsList", spot_edge)
    grow = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        grow, "F", f"Min({args.max_size}, {edge_size} + {args.growth} * F{dist})")
    gmsh.model.mesh.field.setAsBackgroundMesh(grow)
    for opt, val in (
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeMax", args.max_size),
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
    print(f"wrote {args.msh}: {n_quad} quadrilaterals in the film, {n_tri} triangles")


def to_exodus(args):
    """Convert to Exodus so the block and sideset names reach MFEM.

    MFEM picks up named attribute sets from Exodus only; its gmsh reader gives
    bare integer attributes. MOOSE reads the physical names out of the .msh and
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
    conv = os.path.join(here, ".electrical_contact_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by electrical_contact_mesh.py; converts the gmsh mesh to Exodus so\n"
            "# that the block and sideset names reach MFEM, which reads named attribute\n"
            "# sets from Exodus only.\n"
            "[Mesh]\n  [fmg]\n    type = FileMeshGenerator\n"
            f"    file = {os.path.basename(args.msh)}\n  []\n[]\n"
            "[Problem]\n  solve = false\n[]\n[Executioner]\n  type = Steady\n[]\n"
        )
    try:
        subprocess.run(
            [exe, "-i", os.path.basename(conv), "--mesh-only",
             os.path.abspath(args.exodus)],
            cwd=here, check=True,
            stdout=None if args.verbose else subprocess.DEVNULL,
        )
    finally:
        os.path.exists(conv) and os.remove(conv)
    print(f"wrote {args.exodus}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spot-radius", type=float, default=1.0,
                   help="radius of the contact spot, a (mm)")
    p.add_argument("--rod-radius", type=float, default=5.0, help="radius of the rods, b (mm)")
    p.add_argument("--rod-half-length", type=float, default=20.0,
                   help="distance from the contact plane to each electrode (mm)")
    p.add_argument("--film-thickness", type=float, default=1e-3,
                   help="thickness of the film over the spot, and of the gap (mm)")
    p.add_argument("--film-elements", type=int, default=1,
                   help="elements across the thickness of the film")
    p.add_argument("--edge-size", type=float, default=0.25,
                   help="size of the elements at the edge of the spot, as a multiple "
                        "of the film thickness")
    p.add_argument("--growth", type=float, default=0.1,
                   help="rate at which the element size grows with distance from the "
                        "edge of the spot")
    p.add_argument("--max-size", type=float, default=1.0, help="largest element size (mm)")
    p.add_argument("--msh", default="electrical_contact.msh")
    p.add_argument("--exodus", default="electrical_contact.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
