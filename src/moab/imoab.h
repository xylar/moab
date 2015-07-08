
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
  \param (in) readOptions additional options for reading 
  \param (in) comm   MPI communicator (is this needed if already passed at registration?)
  \param (in) ghost_layers  number of layers 
  \param (in) len_filename  filename length
  \param (in) len_options  read options length

  (this will exchange ghosts and exchange all important tags, like 
   global id, material(block) tags, neumann tags and dirichlett tags
  or should the exchange happen explicitly for the tags user specifies? )
*/
ErrorCode LoadMesh(int * pid, char * filename, char * readOptions,  MPI_Comm * comm, int * ghost_layers, int len_filename, int len_options);

/**
  write mesh (collective)
  \param (in) pid application id 
  \param (in) filename 
  \param (in) writeOptions additional options for writing 
  \param (in) comm   MPI communicator (is this needed if already passed at registration?)
  \param (in) len_filename  filename length
  \param (in) len_options  write options length
( we write one single file; in serial, it will write one file per task)
*/
ErrorCode WriteMesh(int * pid, char * filename, char * writeOptions,  MPI_Comm * comm, int len_filename, int len_options);

/**
  obtain local mesh size information 
  \param (in) pid  application id
  \param (out) VisibleVertices number of vertices on current process (including ghosts + shared)
  \param (out) VisibleElements number of elements on current process (including ghosts )
  \param (out) VisibleBlocks number of visible material sets in local mesh (that includes ghosts)
  \param (out) VisibleSurfaceBC (is this the count of surface elem that have a bc?) 
                                    is this including ghosts or not?
  \param (out) VisibleVertexBC  (is this the count of vertices that have a BC?)
                                    is this including ghosts or not?
*/

ErrorCode GetMeshInfo(int *pid, int * VisibleVertices, int * VisibleElements, int *VisibleBlocks,
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
  GetVertexOwnership
  get mesh rank (processor that owns it) for each vertex  (local, shared or ghost)
  
  \param (in) pid  application id
  \param (in/out) VisibleGlobalRankID processor rank for each vertex (array allocated by client, should be size VisibleVertices)
  \param (in/out) len  allocated size of array (on input); on output, it should be actual length, VisibleVertices. 
*/
ErrorCode GetVertexOwnership(int * pid, int * VisibleGlobalRankID, int * len);

/**
  GetVertexID
  get global ID for each vertex  (local, shared or ghost)
  
  \param (in) pid  application id
  \param (in/out) VisibleID global ID for each vertex (array allocated by client, should be size VisibleVertices)
  \param (in/out) len  allocated size of array (on input); on output, it should be actual length, VisibleVertices. 
*/
ErrorCode GetVertexID(int * pid, int * VisibleID, int * len);

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
  GetElementConnectivity
  get element connectivity for block
  \param (in) pid  application id
  \param (in) Block block ID (index from 1 to VisibleBlocks? ) 
  \param (in/out) Connectivity  connectivity array (allocated by client, size 
                VerticesPerElement*NumElements; it will be local index in coords array)
  \param (in/out) len  allocated size of array (on input); 
                on output, it should be actual length, VerticesPerElement*NumElements. 
*/
ErrorCode GetElementConnectivity(int *pid, int *Block, int * Connectivity, int * len);

/**
  GetElementOwnership
   get element ownership information 
  \param (in) pid  application id
  \param (in) block ID (num from 1 to VisibleBlocks) : we don't know how many global blocks
  \param (in/out) ownership array (allocated by client, size NumElements) 
  \param (in/out) len  allocated size of array (on input); on output, actual length (should be NumElements) 
*/
ErrorCode GetElementOwnership(int * pid, int * Block,int * ElementRankID, int *len);

/**
  GetElementID
   get element ownership information 
  \param (in) pid  application id
  \param (in) Block:  block ID (num from 1 to VisibleBlocks) : we don't know how many global blocks
  \param (in/out) ElementID  global ID for each element (allocated by client, size NumElements) 
                      (this will be global ID in moab terms)
  \param (in/out) len  allocated size of array (on input); 
*/
ErrorCode GetElementID(int * pid, int * Block, int * ElementID, int * len);

/**
  GetPointerToSurfaceBC
   surface boundary condition information
   (all arrays allocated by client, size VisibleSurfaceBC?)

  \param (in) pid  application id
  \param (in/out) ElementID element global id ( corresponding to moab global ID ) 
  \param (in/out) ReferenceSurfaceID, (from 1 to 6 for hex, 1-4 for tetras)  side number 
  \param (in/out) BoundaryConditionType boundary condition type ( a number corresponding to NeumannSet ?)
  \param (in/out) len  allocated size of both arrays (on input); 
*/
ErrorCode GetPointerToSurfaceBC(int *pid, int * ElementID,int * ReferenceSurfaceID,
  int* BoundaryConditionType, int * len);

/**
  GetPointerToVertexBC
   vertex boundary condition info
   (all arrays allocated by client, size VisibleVertexBC)
   \param (in)     pid  application id
   \param (in/out) VertexID vertex global id
   \param (in/out) BoundaryConditionType boundary condition type ( a number corresponding to Dirichlet Set ?)
   \param (in/out) len  allocated size of both arrays (on input); 
*/
ErrorCode GetPointerToVertexBC(int *pid, int * VertexID, int * BoundaryConditionType, int * len);


/**
  FIXME
   (in moab, it will create a dense, double tag ; do we care about sparse tags/ bit tags, integer tags, handle tags, etc)
   \param (in) pid  application id
   \param (in) Name  will correspond to name of the tag in moab
   \param (in) VectorType ?? 
           (in moab, tags can be associated to all entities, elements or vertices, or even sets,
            we do not restrict a tag to a specific entity type)
   \param (in) VectorDimensions : corresponding to the size of the tag (how many doubles per entity?)
*/
ErrorCode DefineVectorStorage(int * pid, char * Name, int * VectorType, int * VectorDimensions);

/**
   FIXME
   \param (in) pid  application id
   \param (in) Name  will correspond to name of the tag in moab
   \param (out) Value  double pointer for internal tag memory  
             (it assumes the entities in moab are contiguous, no gaps, and only one entity sequence)
*/
ErrorCode SolutionVectorStorage(int * pid, ,char * Name, double * Value);

/**
   FIXME 
   adjancency calls?
   \param (in) pid  application id
   \param (in) eid element global ID 
   \param (out) numadj number of adjacent elements
   \param (in/out) adjacent elements (in terms of element ID) (size of number of sides?)
*/

ErrorCode AdjacentElements (int *pid, int * eid, int * numadj, int * adjElems); 
