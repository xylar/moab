from pymoab import core,types


# start a MOAB instance
mb = core.Core()

# load the file
mb.load_file("YOUR_MOAB_MESH_FILE_HERE")

# get the root set of the MOAB instance
root_set = mb.get_root_set()

# query the root set for all triangles
tris = mb.get_entities_by_type(root_set, types.MBTRI, recur = True)
print("Found " + str(tris.size()) + " triangles in this model.")

# similar query for vertices
verts = mb.get_entities_by_type(root_set, types.MBVERTEX, recur = True)
print("Found " + str(verts.size()) + " vertices in this model.")

# get the vertex coordinates
vert_coords = mb.get_coords(verts)
print("Coordinates of vertex " + str(verts[0]) + ":")
print(vert_coords[:3])


