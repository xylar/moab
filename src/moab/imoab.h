
/**
  imoab: simple interface to moab
  callable from c, fortran77, fortran90; fortran 2003 ?
 
  pass everything by reference, so we do not have to use %VAL()
*/

/**
  this will create the moab instance, if not created already
  pass pointer to the  number of command line arguments and 
   pointer to command line arguments
*/
ErrorCode initialize_moab(int * pargc, char *** pargv);

/**
  deletes the moab instance 
*/
ErrorCode finalize_moab();

/**
  register application 
  (internally, a mesh set will be associated with this integer; all mesh 
   for this application will reside in this mesh/file set)
  whenever something is required about the mesh, this integer id will need to be passed
  
  \param app_name application name (PROTEUS, NEK5000, etc)
  \param (out) pid  application id pointer 
  \param (in) length of application name 
*/

ErrorCode register_application(char * app_name, int * pid, int len_name);



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
ErrorCode  read_header_info (int *pid, char * filename, int * GlobalVertices, int * GlobalElements, int * NumDimensions, int * NumPartitions, int len_filename);

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
ErrorCode load_mesh(int * pid, char * filename, MPI_Comm * comm, int * ghost_layers, int len_filename);

/**
  obtain local mesh size information
  \param (in) pid  application id
  \param (out) visibleVertices
  \param (out) VisibleBlocks
  \param (out) VisibleSurfaceBC
  \param (out) VisibleVertexBC
*/

ErrorCode get_mesh_info(int *pid, int * visibleVertices, int *VisibleBlocks,
int * VisibleSurfaceBC, int * VisibleVertexBC);

/**
  get vertex coordinates

  \param (in) pid 
  \param (in/out) coords  pointer to memory that will be filled with 
      interleaved coordinates
  \param (in/out) len; at input, usable memory (3*numv?); on output, actual  
*/
ErrorCode get_visible_vertices_coordinates(int *pid, double * coords, int * len);

