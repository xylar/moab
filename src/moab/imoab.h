
/**
  imoab: simple interface to moab
  callable from c, fortran77, fortran90; fortran 2003 ?
 
  Notes:
  pass everything by reference, so we do not have to use %VAL()
  arrays are allocated by the client; 
  pass the pointer to the start of array, and the allocated length ?
  return the filled array, and the actual length (should be 
  most of the time allocated length)
  or should we assume that the user allocated everything fine?
*/


/**
  this will create the moab instance, if not created already
  \param (in) argc   number of command line arguments 
  \param (in) argv   command line arguments
*/

ErrorCode InitializeMoab(int argc, char **argv);

/**
  deletes the moab instance 
*/
ErrorCode FinalizeMoab();

/**
  register application 
  (internally, a mesh set will be associated with this integer; all mesh 
   for this application will reside in this mesh/file set)
  whenever something is required about the mesh, this integer id will need to be passed
  collective  

  \param (in) app_name application name (PROTEUS, NEK5000, etc)
  \param (in) comm  MPI communicator
  \param (out) pid  application id pointer 
  \param (in) length of application name  
*/

ErrorCode RegisterApplication(char * app_name, MPI_Comm * comm, int * pid, int len_name);

/**
  deregister application: delete mesh associated with it (collective)
  \param (in) pid  application id 
*/

ErrorCode DeregisterApplication( int * pid );

/**
  Get global information from the file
  \param (in) pid application id
  \param (in) filename  application mesh file 
  \param (out) GlobalVertices  number of global vertices 
  \param (out) GlobalElements  number of global elements (highest dimension only?) 
  \param (out) NumDimensions ( 2 or 3 ) 
  \param (out) NumPartitions  num partitions in the file
  \param (in) len_filename  file name length
*/

ErrorCode  ReadHeaderInfo (int *pid, char * filename, int * GlobalVertices, int * GlobalElements, int * NumDimensions, int * NumPartitions, int len_filename);

/**
  load mesh and ghost if needed (collective)
  \param (in) pid application id 
  \param (in) filename 
  \param (in) comm   MPI communicator (is this needed if already passed at registration?)
  \param (in) ghost_layers  number of layers 
  \param (in) len_filename  filename length

  (this will exchange ghosts and exchange all important tags, like 
   global id, material(block) tags, neumann tags and dirichlett tags
  or should the exchange happen explicitly for the tags user specifies? )
*/
ErrorCode LoadMesh(int * pid, char * filename, MPI_Comm * comm, int * ghost_layers, int len_filename);

/**
  obtain local mesh size information 
  \param (in) pid  application id
  \param (out) visibleVertices number of local vertices (including ghosts + shared)
  \param (out) VisibleBlocks number of visible material sets in local mesh (that includes ghosts)
  \param (out) VisibleSurfaceBC (is this the count of surface elem that have a bc?) 
                                    is this including ghosts or not?
  \param (out) VisibleVertexBCa (is this the count of vertices that have a BC?)
                                    is this including ghosts or not?
*/

ErrorCode GetMeshInfo(int *pid, int * visibleVertices, int *VisibleBlocks,
int * VisibleSurfaceBC, int * VisibleVertexBC);

/**
  get vertex coordinates

  \param (in) pid  application id
  \param (in/out) coords  pointer to memory that will be filled with 
      interleaved coordinates; client allocates this 
  \param (in/out) len; at input, usable memory (numdim*numv?); on output, actual  
    or is this parameter not needed at all? assume everything is allocated fine?
*/
ErrorCode GetVisibleVerticesCoordinates(int *pid, double * coords, int * len);

/**
  get mesh rank for each vertex  (local, shared or ghost)
  
  \param (in) pid  application id
  \param (in/out) mesh rank for each vertex (array allocated by client, size visibleVertices)
        (should this be long for mesh  > 2B ?)
*/

ErrorCode GetVertexOwnership(int * pid, int * VisibleGlobalRankID);

/** 
  obtain block information
  \param (in) pid  application id
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
  \param (in) pid  application id
  \param (in) block ID (index from 1 to VisibleBlocks? ) 
  \param (in/out) connectivity array (allocated by client, size 
                VerticesPerElement*NumElements 
*/

ErrorCode GetElementConnectivity(int *pid, int *Block, int * Connectivity);

/**
   get element ownership information 
  \param (in) pid  application id
  \param (in) block ID (num from 1 to VisibleBlocks) : we don't know how many global blocks
  \param (in/out) ownership array (allocated by client, size NumElements) 
                      (this will be global ID in moab terms)
*/
ErrorCode GetElementOwnership(int * pid, int * Block,int * ElementRankID);

/**
   surface boundary condition information
   (all arrays allocated by client, size VisibleSurfaceBC?)

  \param (in) pid  application id
   \param (in/out) element global id (mesh_rank)
   \param (in/out) (from 1 to 6 for hex, 1-4 for tetras)  side number 
   \param (in/out) boundary condition type ( a number corresponding to NeumannSet ?)
*/
ErrorCode GetPointerToSurfaceBC(int *pid, int * ElementID,int * ReferenceSurfaceID,
  int* BoundaryConditionType);

/**
   vertex boundary condition info
   (all arrays allocated by client, size VisibleVertexBC)
   \param (in) pid  application id

   \param (in/out) vertex global id (mesh rank?)
   \param (in/out) boundary condition type ( a number corresponding to Dirichlet Set ?)
*/
ErrorCode GetPointerToVertexBC(int *pid, int * VertexID, int * BoundaryConditionType);

 


