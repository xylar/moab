# MOAB TODO LIST

## General

1. Document FileOptions class in UG
2. Add an example for using various tree decompositions.
3. AHF should be updated to handle mixed and polygonal/polyhedral meshes.
4. Make `mbpart` execute in parallel with support for *ParMetis* and *Zoltan*
5. Add tests for MCNP5 meshtal reader. (@pshriwise, @gonuke) Issue #32
6. Explore usage of `clang-tidy` for cleaning up MOAB sources
7. Unify formatting in all sources with `astyle`
8. Update and improve the re-order entities method to have better memory complexity when rearranging entities for maximizing contiguity
9. Add an example to use NearestNeighbor/KNN queries; Optionally adding NanoFlann as a dependency.
10. Add Tag support for `bit`, `bool`, `long`, `float` and `long double` data types
11. Introduce `MOABInitialize(argc,argv)/MOABFinalize()` methods so that we get access to command-line arguments in order to dynamically drive internal behavior exposed through an options database. Note that this routine will not be responsible for `MPI_Init/MPIFinalize` but can initialize global MOAB objects like timers, loggers, error handlers etc.

## Adaptivity

1. Investigate memory management ideas for dynamically resizing the Sequences when performing adaptivity

### Non-Conformal

1. Update adjacency queries to return the right list of elements, when the mesh has hanging nodes. For example, for a non-conformal element (say 1 quad element adjacent to 2 quads in 2-D), the up adjacency for that hanging vertex (to elem) should return 3 elements. The adjacency for that split edges (and refined elements that contain it) should return the same coarse element. And the coarser element adjacency query will return 2 finer elements that share the half-edges.
2. Support standard templates for splitting various types of elements (tri/quad/tet/hex/poly) in the reference frame
3. Ensure defining and manipulating tags on the hanging nodes work correctly

### Conformal

1. Explore automatic quad and triangle refinement schemes such that the meshes remain conformal. Use Mesquite interfaces to apply smoothing techniques in order to improve quality of these refined regions.
2. Extend the algorithms to 3-D (hexes, tetrahedra)

## PyMOAB: Python interface for MOAB

1. User's Guide Section
2. Variable length tag support via the tag_set/get_by_ptr methods
3. Better implementation of the PyMOAB Range class's __str__ and __repr__ methods
4. More examples and usage tutorials for PyMOAB

## Intersection

1. The intersection points are currently not shared correctly between tasks. When intersection point is on the boundary edge (edge that is shared between partitions), the settle_intersection_points method just sends from the owning processor the correct position for the 3d point. It should also send the entity handle, which should be used for shared tag handle. It can never be multi-shared, so this makes it easier. Also, the non-owning task should send the entity handle towards the owner (would be owner), because this intx point will now be shared. Without this addition, parallel file I/O of the intersection mesh will fail.
2. Verify the DoF numbering for FEM weight generation with `mbtempest` tool
3. Remove or minimize the MOAB overload for the mesh loops in TempestOfflineMap.cpp
4. Modify TempestRemap sources so that it is aware of parallel callers (MOAB)

## Local Discretization

1. Add basis function calculations to simplify FEM operator assembly over MOAB meshes for arbitrary element topologies
2. Consolidate ElemUtil and LD functions to create a more unified interface 

## iRel

Currently, iRel usage requires that we match the geometry entity with the mesh set based on **global id** for both the descriptions; i.e., the geometry entity from iGeom/CGM and mesh set from iMesh/MOAB. The mesh sets have to exist already in the mesh file, and the way we do this is through implicit assumptions.

1. Create the stp file and export from Cubit
2. Reset in Cubit, and import the stp file back
3. Mesh the model depending on user specification and save the *cub* file
4. `mbconvert` from *cub* file to *h5m* format
5. Then load the geometry with iGeom, load mesh with iMesh, and call "relate", which will match the geometry entity with the corresponding mesh set based on global IDs of the geometry entity and the mesh set. (they have to match exactly, by ID and dimension)

However, if these models were created through arbitrary workflows like say a geometry model through `Cubit` and a mesh corresponding to that geometry through `Gmsh`, iRel cannot work out the associations.

1. Create a stp file from either Cubit or Gmsh
2) Mesh with any one of the supported mesher like Gmsh or Netgen
3) If Gmsh has the same strategy as Cubit (boundary first, then interior), in principle we could identify all 0d, 1d, 2d, 3d cells in the mesh; also, depending on the format Gmsh has, some "parent" type information could be saved in the mesh itself.
4) Then do a geometric search, to find the associativity relations by "brute" force, if nothing else is possible.

   Brute force means that by looking at a mesh edge, for example, find the curve in the stp file, that is the closest to the mesh edge. It is a complex process, and we will need some iterations to compute the inverse, etc (For any point in space, we could get what is the closest point on a curve in iGeom. There are also some geometry trees that can be used (similar to Kd-tree, for geometry entities, not for mesh). This needs to be investigated further within iGeom and CGM.

   This search procedure might need to happen for every point in the mesh first. Then edges, etc. For a mesh face (a triangle) for example, we could associate it with a geometry face only if all 3 points are "close" to the same geometry face. All those mesh entities associated to a geometry face needs to be put in a mesh set, with the dimension 2, and global id the same as the geometry entity itself. (or just simple relate). For 3d mesh elements, they should be in the interior of a geometry volume. One condition would be that all points are in the interior of the geometric volume.

## Mesquite integration

1. FileTokenizer is one example, used for reading vtk files; Maybe it would be good to unify this, to "maintain" only one; On a glimpse, Mesquite VTK reader allows more datasets. For example structured mesh and rectilinear grids, which are not supported in MOAB for now.
```C
switch( datatype )
{
  case 1: vtk_read_structured_points( tokens, err ); break;
  case 2: vtk_read_structured_grid  ( tokens, err ); break;
  case 3: vtk_read_unstructured_grid( tokens, err ); break;
  case 4: vtk_read_polydata         ( tokens, err ); break;
  case 5: vtk_read_rectilinear_grid ( tokens, err ); break;
  case 6: vtk_read_field            ( tokens, err ); break;
}
```
2. Consolidate usage of MBMesquite::Vector3D and MBMesquite::Matrix3D classes are similar to moab::CartVect and moab::Matrix3
3. On a more general note, the MOAB database is array-based, while Mesquite is designed based on std::vector of elements/vertices. So in comparison, Mesquite is a more heterogeneous data structure, where each element has an explicit connectivity list, as indices in the vertex array, and there is a C++ object (struct) for each element and each vertex; i.e., each vertex stores its upward adjacency in a vector.
4. There are several memory leaks in Mesquite tests and examples. These should be fixed and cleaned up. Valgrind is NOT happy.
5. Add more examples showing Mesquite usage with MOAB
6. Verify mesh optimization algorithms in both serial and parallel and see if all of the indicators work as they should
7. Unify error propagation between MOAB and Mesquite
8. Unify EntityHandle definition between MOAB and Mesquite (void*)

## Hybrid computing

1. Profile first and investigate algorithms that are amenable for OpenMP parallelism
2. Make sure all Core interface functions are thread-safe
3. Add examples that demonstrate hybrid MPI-OpenMP parallelism
