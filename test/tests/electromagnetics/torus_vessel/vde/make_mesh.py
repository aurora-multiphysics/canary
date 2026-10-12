"""
Generates the torus of MOOSE's embedded_concentric_torus.e inside a larger air
region, so that the vessel's inductance is not reduced by the far field boundary.

- TorusCore: torus of major radius 5 m and minor radius 0.5 m, about the z axis.
- TorusSheath: the shell between minor radii 0.5 m and 1 m.
- Exterior: the air in a cube of half-width HALF_WIDTH, whose surface is the
  boundary named Exterior.

Elements are NEAR_SIZE throughout the torus and up to UNIFORM_DISTANCE from the
shell surfaces, then grow linearly with distance to FAR_SIZE at the cube, so the
torus is meshed alike whatever the size of the cube.

Usage: python make_mesh.py HALF_WIDTH OUTPUT.msh [NEAR_SIZE]
"""

import sys

import gmsh

R0 = 5.0
CORE_RADIUS = 0.5
SHEATH_RADIUS = 1.0
# Default element size in and near the torus [m], giving about as many shell
# elements as embedded_concentric_torus.e.
NEAR_SIZE = 0.25
# Element size at the far field boundary [m].
FAR_SIZE = 4.0
# Distance from the shell surfaces within which elements keep NEAR_SIZE [m]. It
# exceeds the largest distance from them of any point in the torus, the core's
# 0.5 m minor radius.
UNIFORM_DISTANCE = 0.6


def main():
    half_width = float(sys.argv[1])
    output = sys.argv[2]
    near_size = float(sys.argv[3]) if len(sys.argv) > 3 else NEAR_SIZE

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("torus_vessel")
    occ = gmsh.model.occ

    core = occ.addTorus(0, 0, 0, R0, CORE_RADIUS)
    outer = occ.addTorus(0, 0, 0, R0, SHEATH_RADIUS)
    box = occ.addBox(-half_width, -half_width, -half_width, 2 * half_width, 2 * half_width,
                     2 * half_width)
    # Fragmenting makes the three regions share conforming interfaces.
    _, fragment_map = occ.fragment([(3, box)], [(3, outer), (3, core)])
    occ.synchronize()

    # Classify the fragments by their distance from the torus centre line.
    volumes = {"TorusCore": [], "TorusSheath": [], "Exterior": []}
    for dim, tag in gmsh.model.getEntities(3):
        x, y, z = occ.getCenterOfMass(dim, tag)
        mass = occ.getMass(dim, tag)
        if mass > 8 * half_width**3 / 2:
            volumes["Exterior"].append(tag)
        else:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            minor = (zmax - zmin) / 2
            volumes["TorusCore" if minor < 0.75 else "TorusSheath"].append(tag)
    for tag, (name, tags) in enumerate(volumes.items(), start=1):
        if len(tags) != 1:
            raise RuntimeError(f"Expected one {name} volume, found {tags}.")
        gmsh.model.addPhysicalGroup(3, tags, tag, name)

    exterior_surfaces = [
        tag
        for dim, tag in gmsh.model.getBoundary([(3, volumes["Exterior"][0])], oriented=False)
        if abs(abs(gmsh.model.occ.getCenterOfMass(dim, tag)[0]) - half_width) < 1e-6
        or abs(abs(gmsh.model.occ.getCenterOfMass(dim, tag)[1]) - half_width) < 1e-6
        or abs(abs(gmsh.model.occ.getCenterOfMass(dim, tag)[2]) - half_width) < 1e-6
    ]
    if len(exterior_surfaces) != 6:
        raise RuntimeError(f"Expected six cube faces, found {exterior_surfaces}.")
    gmsh.model.addPhysicalGroup(2, exterior_surfaces, 1, "Exterior")

    # Size field: NEAR_SIZE in and near the torus, growing to FAR_SIZE at the cube.
    sheath_surfaces = [
        tag
        for dim, tag in gmsh.model.getBoundary([(3, volumes["TorusSheath"][0])], oriented=False)
    ]
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "SurfacesList", sheath_surfaces)
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", near_size)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", FAR_SIZE)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", UNIFORM_DISTANCE)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", half_width - R0 - SHEATH_RADIUS)
    gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    # MFEM reads version 2.2 of the Gmsh format.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(output)

    for name, tags in volumes.items():
        _, elements, _ = gmsh.model.mesh.getElements(3, tags[0])
        print(f"{name}: {sum(len(e) for e in elements)} elements")
    gmsh.finalize()


if __name__ == "__main__":
    main()
