#!/usr/bin/env python3
"""Build the mesh for the tokamak wall induction example with gmsh.

A torus of plasma current sits inside a toroidal conducting wall, the whole
surrounded by a cylinder of air. The example is documented in
doc/content/problems/tokamak_wall_induction.md.

The geometry is axisymmetric, so the mesh is built by revolving a mesh of the
poloidal (R, z) half-plane about the z axis. Revolving gives elements that are
long in the toroidal direction, along which the axisymmetric solution does not
vary, and short across the thin wall, which a tetrahedral mesh generator working
in 3-D would not produce: it would fill the wall with elements as small across
as they are through it. The prisms of the revolved layers are split into
tetrahedra so that they conform with the unstructured column of air around the
axis, R < axis_column_radius, which is not revolved because elements touching
the axis would be degenerate.

The wall annulus is meshed with a structured grid of triangles,
wall_layers elements through its thickness, so that the thickness is resolved
without the wall's element size leaking into the vacuum either side of it.

Blocks: plasma (r < a), vacuum (a < r < b), wall (b < r < b + d) and air
(everything else), with r the minor radius about the circle R = R0, z = 0.
Sideset: far_field, the outer surface of the air cylinder.

Run with --help for the dimensions and meshing controls. All lengths are in
metres. --refine scales every element size and the number of toroidal layers
together, for mesh convergence studies.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys


def build(args):
    import gmsh

    R0, a, b, d = args.major_radius, args.plasma_radius, args.wall_radius, args.wall_thickness
    rc, rout, zout = args.axis_column_radius, args.air_radius, args.air_half_height
    if not (a < b and R0 - b - d > rc and R0 + b + d < rout and b + d < zout):
        sys.exit("the plasma must lie inside the wall, and the wall clear of the axis column "
                 "and inside the air cylinder")

    f = args.refine
    n_toroidal = 4 * max(1, round(args.toroidal_layers * f / 4))
    wall_layers = max(1, round(args.wall_layers * f))
    # Elements around each quarter of the wall, sized to match the vacuum on its inside.
    wall_arc = max(2, round(0.5 * math.pi * b / (args.vacuum_size / f)))
    h_plasma, h_vacuum, h_air = args.plasma_size / f, args.vacuum_size / f, args.air_size / f

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("tokamak_wall")
    g = gmsh.model.geo

    # The poloidal half-plane is the xz plane, x >= 0, with x playing the part of R.
    centre = g.addPoint(R0, 0, 0)

    def circle(r, h):
        """Four quarter arcs of the circle of minor radius r, anticlockwise from the
        outboard midplane, and the points between them."""
        pts = [g.addPoint(R0 + r * math.cos(k * math.pi / 2), 0, r * math.sin(k * math.pi / 2), h)
               for k in range(4)]
        arcs = [g.addCircleArc(pts[k], centre, pts[(k + 1) % 4]) for k in range(4)]
        return pts, arcs

    _, plasma_arcs = circle(a, h_plasma)
    in_pts, in_arcs = circle(b, h_vacuum)
    out_pts, out_arcs = circle(b + d, h_vacuum)

    corners = [g.addPoint(rc, 0, -zout, h_air), g.addPoint(rout, 0, -zout, h_air),
               g.addPoint(rout, 0, zout, h_air), g.addPoint(rc, 0, zout, h_air)]
    # Bottom, outer side, top and the side facing the axis column, in that order.
    box = [g.addLine(corners[k], corners[(k + 1) % 4]) for k in range(4)]

    plasma = g.addPlaneSurface([g.addCurveLoop(plasma_arcs)])
    vacuum = g.addPlaneSurface([g.addCurveLoop(in_arcs), g.addCurveLoop(plasma_arcs)])
    air = g.addPlaneSurface([g.addCurveLoop(box), g.addCurveLoop(out_arcs)])

    # The wall annulus as four structured quarters.
    radial = [g.addLine(in_pts[k], out_pts[k]) for k in range(4)]
    wall = []
    for k in range(4):
        loop = g.addCurveLoop([radial[k], out_arcs[k], -radial[(k + 1) % 4], -in_arcs[k]])
        wall.append(g.addPlaneSurface([loop]))
    for c in radial:
        g.mesh.setTransfiniteCurve(c, wall_layers + 1)
    for c in in_arcs + out_arcs:
        g.mesh.setTransfiniteCurve(c, wall_arc + 1)
    for s in wall:
        g.mesh.setTransfiniteSurface(s)

    # Revolve in four quarter turns, since the built-in kernel cannot revolve by pi or
    # more at once. Each quarter is revolved from the end faces of the last, so that
    # the layers line up, and the faces left at the end of the fourth, which coincide
    # with the ones the first started from, are merged with them afterwards.
    sources = {"plasma": [plasma], "vacuum": [vacuum], "wall": wall, "air": [air]}
    volumes = {name: [] for name in sources}
    far_field, column_side = [], []
    for _ in range(4):
        for name, surfaces in sources.items():
            ends = []
            for s in surfaces:
                out = g.revolve([(2, s)], 0, 0, 0, 0, 0, 1, math.pi / 2,
                                numElements=[n_toroidal // 4])
                ends.append(out[0][1])
                volumes[name].append(out[1][1])
                if name == "air":
                    # The swept faces follow the source surface's boundary curves, the
                    # four sides of the box coming first.
                    bottom, outer, top, inner = (t for (_, t) in out[2:6])
                    far_field += [bottom, outer, top]
                    column_side.append(inner)
            sources[name] = ends
    g.synchronize()
    g.removeAllDuplicates()
    g.synchronize()

    # Close off the axis column with discs bounded by the arcs that the top and bottom
    # corners of the box swept out, and fill it with unstructured tetrahedra.
    top_arcs, bottom_arcs = [], []
    for s in column_side:
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            zmin, zmax = gmsh.model.getBoundingBox(1, c)[2::3]
            if abs(zmax - zmin) < 1e-9 * zout:
                (top_arcs if zmin > 0 else bottom_arcs).append(c)
    discs = [g.addPlaneSurface([g.addCurveLoop(arcs, reorient=True)])
             for arcs in (top_arcs, bottom_arcs)]
    column = g.addVolume([g.addSurfaceLoop(column_side + discs)])
    volumes["air"].append(column)
    far_field += discs
    g.synchronize()

    for name, vols in volumes.items():
        gmsh.model.addPhysicalGroup(3, vols, name=name)
    gmsh.model.addPhysicalGroup(2, far_field, name="far_field")

    # Grade the air from the size at the wall out to air_size at the far field.
    field = gmsh.model.mesh.field
    distance = field.add("Distance")
    field.setNumbers(distance, "CurvesList", out_arcs)
    field.setNumber(distance, "Sampling", 100)
    threshold = field.add("Threshold")
    field.setNumber(threshold, "InField", distance)
    field.setNumber(threshold, "SizeMin", h_vacuum)
    field.setNumber(threshold, "SizeMax", h_air)
    field.setNumber(threshold, "DistMin", 0.0)
    field.setNumber(threshold, "DistMax", args.grading_distance)
    # The field only applies in the air; the plasma and vacuum keep their own sizes.
    restrict = field.add("Restrict")
    field.setNumber(restrict, "InField", threshold)
    field.setNumbers(restrict, "SurfacesList", [air])
    field.setNumbers(restrict, "VolumesList", [column])
    field.setAsBackgroundMesh(restrict)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    counts = {}
    for name, vols in volumes.items():
        counts[name] = sum(len(gmsh.model.mesh.getElements(3, v)[1][0]) for v in vols)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(args.msh)
    gmsh.finalize()
    summary = ", ".join(f"{n} {c}" for n, c in counts.items())
    print(f"wrote {args.msh}: {sum(counts.values())} tetrahedra ({summary}), "
          f"{n_toroidal} toroidal layers")


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
    conv = os.path.join(here, ".tokamak_wall_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by tokamak_wall_mesh.py; converts the gmsh mesh to Exodus so that\n"
            "# the block and sideset names reach MFEM, which reads named attribute sets\n"
            "# from Exodus only.\n"
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
    geom = p.add_argument_group("geometry (m)")
    geom.add_argument("--major-radius", type=float, default=1.0,
                      help="R0, the major radius of the plasma and of the wall")
    geom.add_argument("--plasma-radius", type=float, default=0.3,
                      help="a, the minor radius of the plasma")
    geom.add_argument("--wall-radius", type=float, default=0.4,
                      help="b, the minor radius of the inside of the wall")
    geom.add_argument("--wall-thickness", type=float, default=0.03, help="d")
    geom.add_argument("--axis-column-radius", type=float, default=0.4,
                      help="radius of the unstructured column of air around the axis")
    geom.add_argument("--air-radius", type=float, default=6.0,
                      help="radius of the air cylinder")
    geom.add_argument("--air-half-height", type=float, default=6.0,
                      help="half-height of the air cylinder")
    mesh = p.add_argument_group("meshing")
    mesh.add_argument("--toroidal-layers", type=int, default=64,
                      help="elements around the torus; rounded to a multiple of four")
    mesh.add_argument("--wall-layers", type=int, default=2,
                      help="elements through the wall thickness")
    mesh.add_argument("--plasma-size", type=float, default=0.04,
                      help="element size in the plasma (m)")
    mesh.add_argument("--vacuum-size", type=float, default=0.04,
                      help="element size in the vacuum between the plasma and the wall, "
                           "and along the wall (m)")
    mesh.add_argument("--air-size", type=float, default=0.75,
                      help="element size at the far field (m)")
    mesh.add_argument("--grading-distance", type=float, default=3.0,
                      help="distance from the wall over which the air coarsens (m)")
    mesh.add_argument("--refine", type=float, default=1.0,
                      help="divide every element size by this, and multiply the number of "
                           "toroidal and wall layers by it")
    p.add_argument("--msh", default="tokamak_wall.msh")
    p.add_argument("--exodus", default="tokamak_wall.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
