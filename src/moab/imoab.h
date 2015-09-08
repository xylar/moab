#ifndef IMOAB_H
#define IMOAB_H
/** \file iMOAB.h
  iMOAB: a language-agnostic, lightweight interface to MOAB

  Supports usage from C/C++, Fortran (77/90/2003), Python
 
  \remark 1) All data in the interface are exposed via POD-types.
  \remark 2) Pass everything by reference, so we do not have to use %VAL()
  \remark 3) Arrays are allocated by the client code. No concerns about 
     de-allocation of the data will be taken up by the interface.
  \remark 4) Always pass the pointer to the start of array along with the
     total allocated size for the array.
  \remark 5) Return the filled array requested by client along with 
     optionally the actual length of the array that was filled.
     (for typical cases, should be the allocated length)
*/

/*
Comments from Mike and Emily:

1) Fortran MPI_Comm won't work. Take an integer argument and use MPI_F2C calls to get the C-Comm object
2) ReadHeaderInfo - Does it need the pid ? 
3) Reuse the comm object form the registration for both load and write ops. Do not take comm objects again to avoid confusion.
4) Decipher the global comm object and the subset partioning for each application based on the comm object
5) GetMeshInfo - return separately the owned and ghosted vertices/elements -- not together in visible_** but rather owned_** and ghosted_**. Make these arrays of size 2.
6) Should we sort the vertices after we do ghosting ? So that the locally owned is first and then the ghosted is appended.
7) RCM only for hte owned part of the mesh -- do not screw with the ghosted layers
8) GetBlockID - identical to GetVertexID -- return the global numberign for block
9) GetVertexID -- remember that the order of vertices returned have an implicit numbering embedded in it. DO NOT CHANGE THIS ORDERING...
10) GetBlockInfo takes global Block ID; Remove blockname unless there is a separate use case for it..
11) GetElementConnectivity - clarify whether we return global or local vertex numbering. Preferably local numbering else lot of deciphering for global.
*/

#define iMOAB_AppID    int*
#define iMOAB_String   char*
#define iMOAB_GlobalID int
#define iMOAB_LocalID  int
#define ErrCode        int

/** 
  \fn ErrCode iMOABInitialize( int argc, iMOAB_String* argv )
  \brief Initialize the iMOAB interface implementation and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective

  \param[in] argc (int)           Number of command line arguments 
  \param[in] argv (iMOAB_String*) Command line arguments
*/
ErrCode iMOABInitialize( int argc, iMOAB_String* argv );

/** 
  \fn ErrCode iMOABInitializeFortran( )
  \brief Initialize the iMOAB interface implementation from Fortran driver and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective
*/
ErrCode iMOABInitializeFortran( );

/**
  \fn ErrCode iMOABFinalize()
  \brief Finalize the iMOAB interface implementation and delete the internally reference counted MOAB instance.

  <B>Operations:</B> Collective
*/
ErrCode iMOABFinalize();

/**
  \fn ErrCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid )
  \brief Register application - Create a unique application ID and bootstrap interfaces for further queries.
  
  \note
  Internally, a mesh set will be associated with the application ID and all subsequent queries on the MOAB
  instance will be directed to this mesh/file set.

  <B>Operations:</B> Collective

  \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
  \param[in]  comm (MPI_Comm*)        MPI communicator to be used for all mesh-releated queries originating from this application
  \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
*/
ErrCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid );

/**
  \fn ErrCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length )
  \brief Register a Fortran-basedapplication - Create a unique application ID and bootstrap interfaces for further queries.
  
  \note
  Internally, the Comm object will be converted and stored as MPI_Comm. Additionally, a mesh set will be associated with the 
  application ID and all subsequent queries on the MOAB instance will be directed to this mesh/file set.

  <B>Operations:</B> Collective

  \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
  \param[in]  comm (int*)             MPI communicator to be used for all mesh-releated queries originating from this application
  \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
  \param[in]  app_name_length (int)   Length of application name string
*/
ErrCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length );

/**
  \fn ErrCode DeregisterApplication( iMOAB_AppID pid )
  \brief De-Register application: delete mesh (set) associated with the application ID

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID) The unique pointer to the application ID
*/
ErrCode DeregisterApplication( iMOAB_AppID pid );

/**
  \fn ErrCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
  \brief Get global information from the file

  <B>Operations:</B> Not collective

  \param[in]  filename (iMOAB_String)    The MOAB mesh file (H5M) to probe for header information
  \param[out] num_global_vertices (int*) The total number of vertices in the mesh file
  \param[out] num_global_elements (int*) The total number of elements (of highest dimension only) 
  \param[out] num_dimension (int*)       The highest dimension of elements in the mesh (Edge=1, Tri/Quad=2, Tet/Hex/Prism/Pyramid=3)
  \param[out] num_parts (int*)           The total number of partitions available in the mesh file, typically partitioned with mbpart during pre-processing
  \param[in]  filename_length (int)      Length of the file name string
*/
ErrCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length );

/**
  \fn ErrCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int num_ghost_layers, int filename_length, int read_options_length )
  \brief Load a MOAB mesh file in parallel and exchange ghost layers as requested

  \note
  This will exchange ghosts and the dense/sparse tags that are specified in the mesh.
  Do we need an interface to exchange tags explicitly that user specifies separately ?
  In which case, do we assume that implicit tags like GLOBAL_ID, MATERIAL_SET, NEUMANN_SET, DIRICHLET_SET are exchanged by default ?

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to load onto the internal application mesh set
  \param[in] read_options (iMOAB_String)  Additional options for reading the MOAB mesh file in parallel 
  \param[in] num_ghost_layers (int)       The total number of ghost layers to exchange during mesh loading
  \param[in] filename_length (int)        Length of the filename string
  \param[in] read_options_length (int)    Length of the read options string  
*/
ErrCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int * num_ghost_layers, int filename_length, int read_options_length );

/**
  \fn ErrCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
  \brief Write a MOAB mesh along with the solution tags to a file

  \note
  The interface will write one single file (H5M) and for serial files (VTK/Exodus), it will write one file per task

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to write all the entities contained in the internal application mesh set
  \param[in] write_options (iMOAB_String) Additional options for writing the MOAB mesh in parallel
  \param[in] filename_length (int*)       Length of the filename string
  \param[in] write_options_length (int*)  Length of the write options string
*/
ErrCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length );

/**
  \fn ErrCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
  \brief Obtain local mesh size information based on the loaded file

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[out] num_visible_vertices (int*)  The number of vertices in the current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_elements (int*)  The number of elements in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_blocks (int*)    The number of material sets in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_surfaceBC (int*) The number of surfaces that have a NEUMANN_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_vertexBC (int*)  The number of vertices that have a DIRICHLET_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
*/
ErrCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC );

/**
  \fn ErrCode GetVertexID( iMOAB_AppID pid, int vertices_length, iMOAB_GlobalID* global_vertex_ID, iMOAB_LocalID* local_vertex_ID )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertices_length (int)               The allocated size of array (typical <TT>size := num_visible_vertices</TT>)
  \param[out] global_vertex_ID (iMOAB_GlobalID*)  The global IDs for all locally visible vertices (array allocated by client)
  \param[out] local_vertex_ID (iMOAB_LocalID*)    (<I><TT>Optional</TT></I>) The local IDs for all locally visible vertices (array allocated by client)
*/
ErrCode GetVertexID( iMOAB_AppID pid, int vertices_length, iMOAB_GlobalID* global_vertex_ID, iMOAB_LocalID* local_vertex_ID );

/**
  \fn ErrCode GetVertexOwnership( iMOAB_AppID pid, int vertices_length, int* visible_global_rank_ID )
  \brief Get vertex ownership information i.e., for each vertex based on the local ID, return the process that owns the vertex (local, shared or ghost)

  \note
  Do we need to implement this for owned ? That doesn't make sense.
  If we query only for shared, how do we relate the ordering ?

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)             The unique pointer to the application ID
  \param[in]  vertices_length (int)         The allocated size of array (typically <TT>size := num_visible_vertices</TT>)
  \param[out] visible_global_rank_ID (int*) The processor rank owning each of the local vertices 
*/
ErrCode GetVertexOwnership( iMOAB_AppID pid, int vertices_length, int* visible_global_rank_ID );

/**
  \fn ErrCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int coords_length, double* coordinates )
  \brief Get vertex coordinates for all local (owned and ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)     The unique pointer to the application ID
  \param[in]  coords_length (int)   The size of the allocated coordinate array (array allocated by client, <TT>size := 3*num_visible_vertices</TT>)
  \param[out] coordinates (double*) The pointer to client allocated memory that will be filled with interleaved coordinates (do need an option for blocked coordinates ?)
*/
ErrCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int coords_length, double* coordinates );

/**
  \fn ErrCode GetBlockID( iMOAB_AppID pid, int block_length, iMOAB_GlobalID* global_block_IDs, iMOAB_LocalID* local_block_IDs )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID
  \param[in]  block_length (int)                 The allocated size of array (typical <TT>size := num_visible_blocks</TT>)
  \param[out] global_block_IDs (iMOAB_GlobalID*) The global IDs for all locally visible blocks (array allocated by client)
  \param[out] local_block_IDs (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible blocks (array allocated by client)
*/
ErrCode GetBlockID( iMOAB_AppID pid, int block_length, iMOAB_GlobalID* global_block_IDs, iMOAB_LocalID* local_block_IDs );

/**
  \fn ErrCode  GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int* vertices_per_element, int* num_elements_in_block)
  \brief Get the global block information and elements of certain type or belonging to MATERIAL_SET

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set to be queried
  \param[out] vertices_per_element (int*)       The number of vertices per element
  \param[out] num_elements_in_block (int*)      The number of elements in block
*/
ErrCode GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int* vertices_per_element, int* num_elements_in_block);

/** 
  \fn ErrCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int connectivity_length, int* element_connectivity)
  \brief Get the connectivity for elements within a certain block, ordered based on global element IDs

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  connectivity_length (int)         The allocated size of array (typical <TT>size := vertices_per_element*num_visible_elements</TT>)
  \param[out] element_connectivity (int*)       The connectivity array to store element ordering in MOAB canonical numbering scheme (array allocated by client); array contains vertex identifiers with global ID numbering
*/
ErrCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int connectivity_length, int* element_connectivity);

/**
  \fn ErrCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, int* element_ownership)
  \brief Get the element ownership within a certain block i.e., processor ID of the element owner

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)       The allocated size of ownership array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] element_ownership (int*)          The ownership array to store processor ID for all elements (array allocated by client) 
*/
ErrCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, int* element_ownership);

/**
  \fn ErrCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID)
  \brief Get the global IDs for all locally visible elements belonging to a particular block 

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)    The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)         The allocated size of global element ID array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] global_element_ID (iMOAB_GlobalID*) The global IDs for all locally visible elements (array allocated by client)
  \param[out] local_element_ID (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible elements (array allocated by client)
*/
ErrCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID);

/**
  \fn ErrCode GetPointerToSurfaceBC(iMOAB_AppID pid, int surface_BC_length, iMOAB_GlobalID* global_element_ID, int* reference_surface_ID, int* boundary_condition_value)
  \brief Get the surface boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  surface_BC_length (int)             The allocated size of surface boundary condition array, same as <TT>num_visible_surfaceBC</TT> returned by GetMeshInfo()
  \param[out] global_element_ID (iMOAB_GlobalID*) The global element IDs that contains the side with the surface BC
  \param[out] reference_surface_ID (int*)         The surface number with the BC in the canonical reference element (e.g., 1 to 6 for HEX, 1-4 for TET)
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the NeumannSet defined on the element)
*/
ErrCode GetPointerToSurfaceBC(iMOAB_AppID pid, int surface_BC_length, iMOAB_GlobalID* global_element_ID, int* reference_surface_ID, int* boundary_condition_value);

/**
  \fn ErrCode GetPointerToVertexBC(iMOAB_AppID pid, int vertex_BC_length, iMOAB_GlobalID* global_vertext_ID, int* num_vertex_BC, int* boundary_condition_value)
  \brief Get the vertex boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertex_BC_length (int)              The allocated size of vertex boundary condition array, same as <TT>num_visible_vertexBC</TT> returned by GetMeshInfo()
  \param[out] global_vertext_ID (iMOAB_GlobalID*) The global vertex ID that has Dirichlet BC defined
  \param[out] num_vertex_BC (int*)                The allocated size of vertex boundary condition array, same as num_visible_vertexBC
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the DirichletSet defined on the vertex)
*/
ErrCode GetPointerToVertexBC(iMOAB_AppID pid, int vertex_BC_length, iMOAB_GlobalID* global_vertext_ID, int* num_vertex_BC, int* boundary_condition_value);

/**
  \fn ErrCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int tag_storage_name_length)
  \brief Define a MOAB Tag corresponding to the application depending on requested types.

  \note
  In MOAB, for most solution vectors, we only need to create a "Dense", "Double" Tag

  \todo 1) Do we care about sparse/bit/integer/handle tags, and variable-length tags ?

  <B>Operations:</B> Collective

   \param[in] pid (iMOAB_AppID)               The unique pointer to the application ID
   \param[in] tag_storage_name (iMOAB_String) The tag name to store/retreive the data in MOAB
   \param[in] tag_type (int*)                 The type of MOAB tag (Dense/Sparse on Vertices/Elements, Double/Int/EntityHandle)
   \param[in] components_per_entity (int*)    The total size of vector dimension per entity for the tag (e.g., number of doubles per entity)
   \param[in] tag_storage_name_length (int)   The length of the tag_storage_name string
*/
ErrCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int tag_storage_name_length);

/**
   \fn ErrCode SetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to replace the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrCode SetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)  The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)     The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (double*)       The array data of type <I>double</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (int)    The length of the tag_storage_name string
*/
ErrCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrCode GetNeighborElements(iMOAB_AppID pid, iMOAB_GlobalID global_element_ID, int* num_adjacent_elements, iMOAB_GlobalID* adjacent_element_IDs)
   \brief Compute the adjacencies for the element entities

   <B>Operations:</B> Collective

   \bug 1) How do we decide the elements/vertices to query if we are not passing EntityHandles ?
   \bug 2) We need a dimension for up/down adjacencies. Or do we care only about 1-ring neighbors ?

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  global_element_ID (iMOAB_GlobalID)     The global element ID for which adjacency information is needed
   \param[out] num_adjacent_elements (int*)           The total number of adjacent elements
   \param[out] adjacent_element_IDs (iMOAB_GlobalID*) The global element IDs of all adjacent elements to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrCode GetNeighborElements(iMOAB_AppID pid, iMOAB_GlobalID global_element_ID, int* num_adjacent_elements, iMOAB_GlobalID* adjacent_element_IDs);

/**
   \fn ErrCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_GlobalID global_vertex_ID, int* num_adjacent_vertices, iMOAB_GlobalID* adjacent_vertex_IDs)
   \brief Compute the adjacencies for the vertex entities

   <B>Operations:</B> Collective

   \bug 1) How do we decide the elements/vertices to query if we are not passing EntityHandles ?
   \bug 2) We need a dimension for up/down adjacencies. Or do we care only about 1-ring neighbors ?

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  global_vertex_ID (iMOAB_GlobalID)      The global vertex ID for which adjacency information is needed
   \param[out] num_adjacent_vertices (int*)           The total number of adjacent vertices
   \param[out] adjacent_vertex_IDs (iMOAB_GlobalID*)  The global element IDs of all adjacent vertices to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_GlobalID global_vertex_ID, int* num_adjacent_vertices, iMOAB_GlobalID* adjacent_vertex_IDs);

#endif
