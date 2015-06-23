
/**
  imoab: simple interface to moab
  callable from c, fortran77, fortran90; fortran 2003 ?
 
  Notes:
  pass everything by reference, so we do not have to use %VAL()
  arrays are allocated by the client; 
  pass the pointer to the start of array, and the allocated length
  return the filled array, and the actual length (should be 
  most of the time allocated length)
*/

/**
  this will create the moab instance, if not created already
  pass pointer to the  number of command line arguments and 
   pointer to command line arguments
*/
ErrorCode InitializeMoab(int * pargc, char *** pargv);

/**
  deletes the moab instance 
*/
ErrorCode FinalizeMoab();

/**
  register application 
  (internally, a mesh set will be associated with this integer; all mesh 
   for this application will reside in this mesh/file set)
  whenever something is required about the mesh, this integer id will need to be passed
  
  \param app_name application name (PROTEUS, NEK5000, etc)
  \param (out) pid  application id pointer 
  \param (in) length of application name 
*/

ErrorCode RegisterApplication(char * app_name, int * pid, int len_name);



/**
  Get global information from the file
  \param (in) pid application id
  \param filename 
  \param (out) GlobalVertices
  \param (out) GlobalElements (highest dimension only?)
  \param (out) NumDimensions ( 2 or 3 ) 
  \param (out) NumPartitions 
  \param (in) file name length
*/
ErrorCode  ReadHeaderInfo (int *pid, char * filename, int * GlobalVertices, int * GlobalElements, int * NumDimensions, int * NumPartitions, int len_filename);

/**
  load mesh and ghost if needed
  \param (in) pid application id 
  \param (in) filename 
  \param (in) comm   MPI communicator
  \param (in) ghost_layers  
  \param (in) len_filename 

  (this will exchange ghosts and exchange all important tags, like 
   global id, material(block) tags, neumann tags and dirichlett tags
  or should the exchange happen explicitly for the tags user specifies? )
*/
ErrorCode LoadMesh(int * pid, char * filename, MPI_Comm * comm, int * ghost_layers, int len_filename);

/**
  obtain local mesh size information
  \param (in) pid  application id
  \param (out) visibleVertices
  \param (out) VisibleBlocks
  \param (out) VisibleSurfaceBC (is this the count of surface elem that have a bc?) 
  \param (out) VisibleVertexBCa (is this the count of vertices that have a BC?)
*/

ErrorCode GetMeshInfo(int *pid, int * visibleVertices, int *VisibleBlocks,
int * VisibleSurfaceBC, int * VisibleVertexBC);

/**
  get vertex coordinates

  \param (in) pid 
  \param (in/out) coords  pointer to memory that will be filled with 
      interleaved coordinates; client allocates this 
  \param (in/out) len; at input, usable memory (3*numv?); on output, actual  
*/
ErrorCode GetVisibleVerticesCoordinates(int *pid, double * coords, int * len);

/**
  get rank that owns each vertex 
  
  \param (in) pid
  \param (in/out) mesh rank for each vertex (array allocated by client, size visibleVertices)
        (should this be long for mesh  > 2B ?)
*/

ErrorCode GetVertexOwnership(int * pid, int * VisibleGlobalRankID);

/** 
  obtain block information
  \param (in) pid 
  \param (in) Block  block ID
  \param (out) VerticesPerElement  number of vertices per element
  \param (out) NumElements          number of elements in block
  \param (out) BlockName  return for the material set the name (if given as NAME in h5m file?)
  \param (out) lenBlockName  name length 
*/
  
ErrorCode  GetBlockInfo(int *pid, int * Block, int * VerticesPerElement,
    int * NumElements, char * BlockName, int lenBlockName);

/** 
  get element connectivity for block
  \param (in) pid
  \param (in) block ID
  \param (in/out) connectivity array (allocated by client, size 
                VerticesPerElement*NumElements 
*/

ErrorCode GetElementConnectivity(int *pid, int *Block, int * Connectivity);

/**
   get element ownership information 
  \param (in) pid
  \param (in) block ID
  \param (in/out) ownership array (allocated by client, size NumElements) 
                      (this will be global ID in moab terms)
*/
ErrorCode GetElementOwnership(int * pid, int * Block,int * ElementRankID);

/**
   surface boundary condition information
   (all arrays allocated by client, size VisibleSurfaceBC?)

   \param (in) pid
   \param (out) element global id (mesh_rank)
   \param (out) (from 1 to 6 for hex, 1-4 for tetras)  side number 
   \param (out) boundary condition type ( a number corresponding to NeumannSet ?)
*/
ErrorCode GetPointerToSurfaceBC(int *pid, int * ElementID,int * ReferenceSurfaceID,
  int* BoundaryConditionType);

/**
   vertex boundary condition info
   \param (in) pid
   (all arrays allocated by client, size VisibleVertexBC)

   \param (out) vertex global id (mesh rank?)
   \param (out) boundary condition type ( a number corresponding to Dirichlet Set ?)
*/
ErrorCode GetPointerToVertexBC(int *pid, int * VertexID, int * BoundaryConditionType);








/**
  vertex ownership for all visible vertices

  \param (in)  pid  application id
  \param (in/out) vertex_owner  array with mesh rank id  (will correspond to GLOBAL_ID tag in MOAB,
                      if there are no gaps in GLOBAL_ID tag space )
  \param (in/out) len; at input: usable memory (numv) ; on output, actual length
*/

ErrorCode  get_vertices_ownership(int * pid, int * vertex_ID, int * len);

/**
  \param (in) pid
  \param (in) block (1 to num visiblt blocks)  
  \param (out) blockname : integer corresponding to material set
  \param (out) vert_per_elem
  \param (out) num_elem 

*/
ErrorCode get_block_info(int *pid, int * block, int * blockname, int * vert_per_elem, int * num_elem);

/**
  get connectivity for the whole block
  (represented by indices in the vertex array, 0 or 1 based?)
  \param (in) pid
  \param (in) block 
  \param (in/out) connectivity  
  \param (in/out) len  

  Notes: block type (num of vertices known from before)
  A block corresponds to a material set in moab, and a block in cubit  
*/
ErrorCode get_connectivity (int * pid, int *block, int * connectivity, int *len);

/**
  similar to vertex ownership; is this GLOBAL_ID or processor rank?
*/
ErrorCode get_element_ownership(int *pid, int *block, int * elem_ID, int * len);


/**
  \param pid
  \param bcid boundary condition index
  \param (out) len  number of quads in the set
  output will be length for dimensioning the array for next call
*/
ErrorCode get_surface_bc_len(int *pid, int *bcid, int *len);

/**
 so it would be a set of faces in a neumann set
 output will be element mesh rank and face index (1 to 6) for each quad in the original set
  \param (out) elementID  element mesh rank (or element index in local array?)
  \param (out) referenceSurface face index from 1 to 6 for each 
 these arrays need to be dimensioned by the client
*/

ErrorCode get_surface_bc(int *pid, int *bcid, int * elementID, int* referenceSurface, int *len);

/**
  \param pid
  \param bcid boundary condition index
  \param (out) len  number of vertices in the set
  output will be length for dimensioning the array for next call
*/
ErrorCode get_vertex_bc_len(int *pid, int *bcid, int *len);

/**
 so it would be a set of faces in a dirichlett set
 output will be vertex rank in the original set
  \param (out) vertexID  vertex  mesh rank (or vertex index in local array?)
 these arrays need to be dimensioned by the client
*/

ErrorCode get_vertex_bc(int *pid, int *bcid, int * vertexID, int *len);
 


