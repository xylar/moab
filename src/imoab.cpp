/** \file iMOAB.cpp
*/

#include "moab/Core.hpp"
using namespace moab;
#include "mpi.h"

#include "moab/imoab.h"
/*
this mhdf.h is not part of moab installation, but it is part of moab library
copy it in this folder (imoab/src/mhdf) temporarily; after imoab is part of moab, it is not neded 
*/
#include "moab/mhdf.h"
#include <stdio.h>
/*
 this is needed so far because of direct access to hdf5/mhdf
  */

#include <H5Tpublic.h>

#include <iostream>
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include <sstream>

// global variables ; should they be organized in a structure, for easier references?
// or how do we keep them global?

Interface * MBI = 0;
// we should also have the default tags stored, initialized
Tag gtags[5]; // material, neumann, dirichlet, partition tag, globalID
// should this be part of init moab?

struct appData {
  EntityHandle file_set;
  Range all_verts;
  Range primary_elems;
  Range mat_sets;
  std::map<int, int> matIndex; // map from global block id to index in mat_sets
  Range neu_sets;
  Range diri_sets;
};

// are there reasons to have multiple moab inits? Is ref count needed?
int refCountMB( 0) ;
int iArgc;
iMOAB_String * iArgv;

/*
 list of moab entity sets corresponding to each application and pid
 */
int unused_pid =0;
// std::vector<EntityHandle>  app_FileSets; // in order of creation
std::map<std::string, int> appIdMap;     // from app string (uppercase) to app id
std::vector<ParallelComm*> pcomms; // created in order of applications, one moab::ParallelComm for each
std::vector<appData> appDatas; // the same order as pcomms

/** 
  \fn ErrorCode iMOABInitialize( int argc, iMOAB_String* argv )
  \brief Initialize the iMOAB interface implementation and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective

  \param[in] argc (int)           Number of command line arguments 
  \param[in] argv (iMOAB_String*) Command line arguments
*/
ErrCode iMOABInitialize( int argc, iMOAB_String* argv )
{
   iArgc = argc;
   iArgv = argv; // shalow copy
   if (0==refCountMB)
   {
     MBI = new Core();
     // retrieve the default tags
     const char* const shared_set_tag_names[] = {MATERIAL_SET_TAG_NAME,
                                                 NEUMANN_SET_TAG_NAME,
                                                 DIRICHLET_SET_TAG_NAME,
                                                 GLOBAL_ID_TAG_NAME};
     // blocks, visible surfaceBC(neumann), vertexBC (Dirichlet), global id, parallel partition
     for (int i = 0; i < 4; i++) {

       ErrorCode rval = MBI->tag_get_handle(shared_set_tag_names[i], 1, MB_TYPE_INTEGER,
                                           gtags[i], MB_TAG_ANY);
       if (MB_SUCCESS!=rval)
         return 1;
     }
   }
   refCountMB++;
   return MB_SUCCESS;
}

#if 0
/** 
  \fn ErrorCode iMOABInitializeFortran( )
  \brief Initialize the iMOAB interface implementation from Fortran driver and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective
*/

ErrorCode iMOABInitializeFortran( );
#endif 

/**
  \fn ErrorCode iMOABFinalize()
  \brief Finalize the iMOAB interface implementation and delete the internally reference counted MOAB instance.

  <B>Operations:</B> Collective
*/
ErrCode iMOABFinalize()
{
   refCountMB--;
   if (0==refCountMB)
      delete MBI; 
   return MB_SUCCESS;
}

/**
  \fn ErrorCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid )
  \brief Register application - Create a unique application ID and bootstrap interfaces for further queries.
  
  \note
  Internally, a mesh set will be associated with the application ID and all subsequent queries on the MOAB
  instance will be directed to this mesh/file set.

  <B>Operations:</B> Collective

  \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
  \param[in]  comm (MPI_Comm*)        MPI communicator to be used for all mesh-releated queries originating from this application
  \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
*/


ErrCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid )
{
  // will create a parallel comm for this application too, so there will be a
  // mapping from *pid to file set and to parallel comm instances
  std::string name(app_name);
  if (appIdMap.find(name)!=appIdMap.end())
  {
    std::cout << " application already registered \n";
    return 1;
  }
  *pid =  unused_pid++;
  appIdMap[name] = *pid;
  // now create ParallelComm and a file set for this application
  ParallelComm * pco = new ParallelComm(MBI, *comm);

#if 1
  int index = pco->get_id(); // t could be useful to get app id from pcomm instance ...
  assert(index==*pid);
#endif
  pcomms.push_back(pco);

  // create now the file set that will be used for loading the model in
  EntityHandle file_set;
  ErrorCode rval = MBI->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != rval )
    return 1;
  appData app_data;
  app_data.file_set=file_set;
  appDatas.push_back(app_data); // it will correspond to app_FileSets[*pid] will be the file set of interest
  return 0;
}
#if 0
/**
  \fn ErrorCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length )
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


ErrorCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length );
#endif
/**
  \fn ErrorCode DeregisterApplication( iMOAB_AppID pid )
  \brief De-Register application: delete mesh (set) associated with the application ID

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID) The unique pointer to the application ID
*/
ErrCode DeregisterApplication( iMOAB_AppID pid )
{
	// the file set , parallel comm are all in vectors indexed by *pid
  // assume we did not delete anything yet
  // *pid will not be reused if we register another application
  ParallelComm * pco = pcomms[*pid];
  // we could get the pco also with
  // ParallelComm * pcomm = ParallelComm::get_pcomm(MBI, *pid);
  EntityHandle fileSet = appDatas[*pid].file_set;
  // get all entities part of the file set
  Range fileents;
  ErrorCode rval = MBI->get_entities_by_handle(fileSet, fileents, /*recursive */true);
  if (MB_SUCCESS != rval )
    return 1;

  fileents.insert(fileSet);

  rval = MBI->get_entities_by_type(fileSet, MBENTITYSET, fileents); // append all mesh sets
  if (MB_SUCCESS != rval )
    return 1;
  delete pco;
  rval = MBI->delete_entities(fileents);

  if (MB_SUCCESS != rval )
    return 1;

  return 0;
}

/**
  \fn ErrorCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
  \brief Get global information from the file

  <B>Operations:</B> Not collective

  \param[in]  filename (iMOAB_String)    The MOAB mesh file (H5M) to probe for header information
  \param[out] num_global_vertices (int*) The total number of vertices in the mesh file
  \param[out] num_global_elements (int*) The total number of elements (of highest dimension only) 
  \param[out] num_dimension (int*)       The highest dimension of elements in the mesh (Edge=1, Tri/Quad=2, Tet/Hex/Prism/Pyramid=3)
  \param[out] num_parts (int*)           The total number of partitions available in the mesh file, typically partitioned with mbpart during pre-processing
  \param[in]  filename_length (int)      Length of the file name string
*/

ErrCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
{
  std::string filen(filename);
  if (filename_length< (int)filen.length())
  {
    filen = filen.substr(0,filename_length);
  }
  *num_global_vertices = 0;
  int edges = 0;
  int faces = 0;
  int regions = 0;
  *num_global_elements =0;
  *num_dimension = 0;
  *num_parts = 0;

  mhdf_FileHandle file;
  mhdf_Status status;
  unsigned long max_id;
  struct mhdf_FileDesc* data;
  /* find PARALLEL_PARTITION tag index */
  const char * pname = "PARALLEL_PARTITION";

  long int nval, junk;
  hid_t table[3];


  file = mhdf_openFile( filen.c_str() , 0, &max_id, -1, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }

  data = mhdf_getFileSummary( file, H5T_NATIVE_ULONG, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }
  *num_dimension = data->nodes.vals_per_ent;
  *num_global_vertices = (int)data->nodes.count;

  for (int i=0; i<data->num_elem_desc; i++)
  {
    struct mhdf_ElemDesc * el_desc = &(data->elems[i]);
    struct mhdf_EntDesc * ent_d = &(el_desc->desc);
    if (0==strcmp(el_desc->type, mhdf_EDGE_TYPE_NAME)) edges += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TRI_TYPE_NAME))  faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_QUAD_TYPE_NAME)) faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYGON_TYPE_NAME)) faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TET_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PYRAMID_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PRISM_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_KNIFE_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_HEX_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYHEDRON_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_SEPTAHEDRON_TYPE_NAME)) regions += ent_d->count;
  }
  for (int i=0; i<data->num_tag_desc; i++)
  {
    struct mhdf_TagDesc * tag_desc = &(data->tags[i]);
    if (strcmp(pname,tag_desc->name)==0)
    {
      /*printf(" tag index %d is parallel partition tag\n", i);*/
      if (tag_desc->have_sparse) {
        mhdf_openSparseTagData(file, pname, &nval, &junk, table, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }
      else
      {
        /* could be dense tags on sets */
        mhdf_openDenseTagData(file, pname, mhdf_set_type_handle(), &nval, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }

      *num_parts = (int)nval;
    }
  }

  // is this required?
  if (edges >0 ){
    *num_dimension = 1; // I don't think it will ever return 1
    *num_global_elements = edges;
  }
  if (faces >0 ){
    *num_dimension = 2;
    *num_global_elements = faces;
  }
  if (regions>0){
    *num_dimension = 3;
    *num_global_elements = regions;
  }
  mhdf_closeFile( file, &status );

  free( data );
  return 0;
}



/**
  \fn ErrorCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int num_ghost_layers, int filename_length, int read_options_length )
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


ErrCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int * num_ghost_layers, int filename_length, int read_options_length )
{


  // make sure we use the file set and pcomm associated with the *pid
  std::ostringstream newopts;
  newopts  << read_options;
  newopts << ";PARALLEL_COMM="<<*pid;
  if (*num_ghost_layers>=1)
  {
    // if we want ghosts, we will want additional entities, the last .1
    // because the addl ents can be edges, faces that are part of the neumann sets
    newopts << ";PARALLEL_GHOSTS=3.0."<<*num_ghost_layers<<".3";
  }
  ErrorCode rval = MBI->load_file(filename, &appDatas[*pid].file_set, newopts.str().c_str());
  if (MB_SUCCESS!=rval)
    return 1;
  int rank = pcomms[*pid]->rank();
  int nprocs=pcomms[*pid]->size();

  // do a global id exchange, for all entities
  // when I run the driver on 4 procs, some global ids are 0, unless I do this
  // so there is a bug somewhere; after this change, global ids seem correct
  // FIXME
  Range ents;
  rval = MBI->get_entities_by_handle(appDatas[*pid].file_set, ents, true);
  if (MB_SUCCESS!=rval)
    return 1;
  rval = pcomms[*pid]->exchange_tags(gtags[3], ents);
  if (MB_SUCCESS!=rval)
    return 1;

#if 1
  // some debugging stuff
  std::ostringstream outfile;
  outfile <<"TaskMesh_n" <<nprocs<<"."<< rank<<".h5m";
  // the mesh contains ghosts too, but they are not part of mat/neumann set
  // write in serial the file, to see what tags are missing
  rval = MBI->write_file(outfile.str().c_str()); // everything on root
  if (MB_SUCCESS!=rval)
    return 1;
#endif
  return 0;
}


/**
  \fn ErrorCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
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
ErrCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
{
  // maybe do some processing of strings and lengths
  // maybe do some options processing?
  ErrorCode rval = MBI->write_file(filename,0, write_options,  &appDatas[*pid].file_set, 1);
  if (MB_SUCCESS!=rval)
    return 1;
  return 0;
}



/**
  \fn ErrorCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
  \brief Obtain local mesh size information based on the loaded file

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[out] num_visible_vertices (int*)  The number of vertices in the current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_elements (int*)  The number of elements in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_blocks (int*)    The number of material sets in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_surfaceBC (int*) The number of surfaces that have a NEUMANN_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_vertexBC (int*)  The number of vertices that have a DIRICHLET_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
*/
ErrCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
{

  // this will include ghost elements
  // we should keep a data structure with mesh, sets, etc, for each pid
  //
  EntityHandle fileSet=appDatas[*pid].file_set;
  ErrorCode rval = MBI->get_entities_by_type(fileSet, MBVERTEX, appDatas[*pid].all_verts, true); // recursive
  if (MB_SUCCESS!=rval)
    return 1;
  *num_visible_vertices = (int) appDatas[*pid].all_verts.size();
  // is dimension 3?
  rval = MBI->get_entities_by_dimension(fileSet, 3, appDatas[*pid].primary_elems, true); // recursive
  if (MB_SUCCESS!=rval)
    return 1;
  *num_visible_elements = (int) appDatas[*pid].primary_elems.size();

  // get all blocks, BCs, etc
  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[0]), 0, 1, appDatas[*pid].mat_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  *num_visible_blocks = (int)appDatas[*pid].mat_sets.size();
  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[1]), 0, 1, appDatas[*pid].neu_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  *num_visible_surfaceBC = (int)appDatas[*pid].neu_sets.size();

  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[2]), 0, 1, appDatas[*pid].diri_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  *num_visible_vertexBC= (int)appDatas[*pid].diri_sets.size();

  return 0;
}


/**
  \fn ErrorCode GetVertexID( iMOAB_AppID pid, int vertices_length, iMOAB_GlobalID* global_vertex_ID, iMOAB_LocalID* local_vertex_ID )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertices_length (int)               The allocated size of array (typical <TT>size := num_visible_vertices</TT>)
  \param[out] global_vertex_ID (iMOAB_GlobalID*)  The global IDs for all locally visible vertices (array allocated by client)
  \param[out] local_vertex_ID (iMOAB_LocalID*)    (<I><TT>Optional</TT></I>) The local IDs for all locally visible vertices (array allocated by client)
*/
ErrCode GetVertexID( iMOAB_AppID pid, int vertices_length, iMOAB_GlobalID* global_vertex_ID, iMOAB_LocalID* local_vertex_ID )
{
//
  Range & verts = appDatas[*pid].all_verts;
  // global id tag is gtags[3]
  ErrorCode rval = MBI->tag_get_data(gtags[3], verts, global_vertex_ID);
  if (MB_SUCCESS!=rval)
    return 1;
  int i=0;
  for (Range::iterator vit=verts.begin(); vit!=verts.end(); vit++, i++)
    local_vertex_ID[i]=i;
  if (i!=vertices_length)
    return 1; // problem with array length
  return 0;
}
/**
  \fn ErrorCode GetVertexOwnership( iMOAB_AppID pid, int vertices_length, int* visible_global_rank_ID )
  \brief Get vertex ownership information i.e., for each vertex based on the local ID, return the process that owns the vertex (local, shared or ghost)

  \note
  Do we need to implement this for owned ? That doesn't make sense.
  If we query only for shared, how do we relate the ordering ?

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)             The unique pointer to the application ID
  \param[in]  vertices_length (int)         The allocated size of array (typically <TT>size := num_visible_vertices</TT>)
  \param[out] visible_global_rank_ID (int*) The processor rank owning each of the local vertices 
*/
ErrCode GetVertexOwnership( iMOAB_AppID pid, int vertices_length, int* visible_global_rank_ID )
{
  Range & verts = appDatas[*pid].all_verts;
  ParallelComm * pco = pcomms[*pid];
  int i=0;
  for (Range::iterator vit=verts.begin(); vit!=verts.end(); vit++, i++)
  {
    ErrorCode rval = pco->  get_owner(*vit, visible_global_rank_ID[i]);
    if (MB_SUCCESS!=rval)
      return 1;
  }
  if (i!=vertices_length)
    return 1; // warning array allocation problem

  return 0;
}

/**
  \fn ErrorCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int coords_length, double* coordinates )
  \brief Get vertex coordinates for all local (owned and ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)     The unique pointer to the application ID
  \param[in]  coords_length (int)   The size of the allocated coordinate array (array allocated by client, <TT>size := 3*num_visible_vertices</TT>)
  \param[out] coordinates (double*) The pointer to client allocated memory that will be filled with interleaved coordinates (do need an option for blocked coordinates ?)
*/

ErrCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int coords_length, double* coordinates )
{
  Range & verts = appDatas[*pid].all_verts;
  // interleaved coordinates, so that means deep copy anyway
  if (coords_length!=3*(int)verts.size())
    return 1;
  ErrorCode rval = MBI->get_coords(verts, coordinates);
  if (MB_SUCCESS!=rval)
    return 1;
  return 0;
}

/**
  \fn ErrorCode GetBlockID( iMOAB_AppID pid, int block_length, iMOAB_GlobalID* global_block_IDs, iMOAB_LocalID* local_block_IDs )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID
  \param[in]  block_length (int)                 The allocated size of array (typical <TT>size := num_visible_blocks</TT>)
  \param[out] global_block_IDs (iMOAB_GlobalID*) The global IDs for all locally visible blocks (array allocated by client)
  \param[out] local_block_IDs (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible blocks (array allocated by client)
*/
ErrCode GetBlockID( iMOAB_AppID pid, int block_length, iMOAB_GlobalID* global_block_IDs, iMOAB_LocalID* local_block_IDs )
{
  // local id blocks? they are counted from 0 to number of visible blocks ...
  // will actually return material set tag value for global
  Range & matSets = appDatas[*pid].mat_sets;
  if (block_length!=(int)matSets.size())
    return 1;
  // return material set tag gtags[0 is material set tag
  ErrorCode rval = MBI->tag_get_data(gtags[0], matSets, global_block_IDs);
  if (MB_SUCCESS!=rval)
    return 1;
  // populate map with index
  std::map <int, int> & matIdx = appDatas[*pid].matIndex;
  //
  for (int i=0; i<(int)matSets.size(); i++)
  {
    local_block_IDs[i]= i; // TODO: do we really need this?
    matIdx[global_block_IDs[i]] = i;
  }
  return 0;
}


/**
  \fn ErrorCode  GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int* vertices_per_element, int* num_elements_in_block)
  \brief Get the global block information and elements of certain type or belonging to MATERIAL_SET

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set to be queried
  \param[out] vertices_per_element (int*)       The number of vertices per element
  \param[out] num_elements_in_block (int*)      The number of elements in block
*/
ErrCode GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID,
    int* vertices_per_element, int* num_elements_in_block)
{
  std::map<int, int> & matMap = appDatas[*pid].matIndex;
  std::map<int,int>::iterator it = matMap.find(global_block_ID);
  if (it==matMap.end())
    return 1; // error in finding block with id
  int blockIndex = matMap[global_block_ID];
  EntityHandle matMeshSet = appDatas[*pid].mat_sets[blockIndex];
  Range blo_elems;
  ErrorCode rval = MBI-> get_entities_by_handle(matMeshSet, blo_elems);
  if (MB_SUCCESS!=rval ||  blo_elems.empty() )
    return 1;

  EntityType type = MBI->type_from_handle(blo_elems[0]);
  if (!blo_elems.all_of_type(type))
    return 1; //not all of same  type

  const EntityHandle * conn=NULL;
  int num_verts=0;
  rval = MBI->get_connectivity(blo_elems[0], conn, num_verts);
  if (MB_SUCCESS!=rval)
    return 1;
  *vertices_per_element=num_verts;
  *num_elements_in_block = (int)blo_elems.size();

  return 0;
}
#if 0
/** 
  \fn ErrorCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int connectivity_length, int* element_connectivity)
  \brief Get the connectivity for elements within a certain block, ordered based on global element IDs

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  connectivity_length (int)         The allocated size of array (typical <TT>size := vertices_per_element*num_visible_elements</TT>)
  \param[out] element_connectivity (int*)       The connectivity array to store element ordering in MOAB canonical numbering scheme (array allocated by client); array contains vertex identifiers with global ID numbering
*/
ErrorCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int connectivity_length, int* element_connectivity);

/**
  \fn ErrorCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, int* element_ownership)
  \brief Get the element ownership within a certain block i.e., processor ID of the element owner

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)       The allocated size of ownership array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] element_ownership (int*)          The ownership array to store processor ID for all elements (array allocated by client) 
*/
ErrorCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, int* element_ownership);

/**
  \fn ErrorCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID)
  \brief Get the global IDs for all locally visible elements belonging to a particular block 

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)    The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)         The allocated size of global element ID array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] global_element_ID (iMOAB_GlobalID*) The global IDs for all locally visible elements (array allocated by client)
  \param[out] local_element_ID (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible elements (array allocated by client)
*/
ErrorCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID);

/**
  \fn ErrorCode GetPointerToSurfaceBC(iMOAB_AppID pid, int surface_BC_length, iMOAB_GlobalID* global_element_ID, int* reference_surface_ID, int* boundary_condition_value)
  \brief Get the surface boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  surface_BC_length (int)             The allocated size of surface boundary condition array, same as <TT>num_visible_surfaceBC</TT> returned by GetMeshInfo()
  \param[out] global_element_ID (iMOAB_GlobalID*) The global element IDs that contains the side with the surface BC
  \param[out] reference_surface_ID (int*)         The surface number with the BC in the canonical reference element (e.g., 1 to 6 for HEX, 1-4 for TET)
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the NeumannSet defined on the element)
*/
ErrorCode GetPointerToSurfaceBC(iMOAB_AppID pid, int surface_BC_length, iMOAB_GlobalID* global_element_ID, int* reference_surface_ID, int* boundary_condition_value);

/**
  \fn ErrorCode GetPointerToVertexBC(iMOAB_AppID pid, int vertex_BC_length, iMOAB_GlobalID* global_vertext_ID, int* num_vertex_BC, int* boundary_condition_value)
  \brief Get the vertex boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertex_BC_length (int)              The allocated size of vertex boundary condition array, same as <TT>num_visible_vertexBC</TT> returned by GetMeshInfo()
  \param[out] global_vertext_ID (iMOAB_GlobalID*) The global vertex ID that has Dirichlet BC defined
  \param[out] num_vertex_BC (int*)                The allocated size of vertex boundary condition array, same as num_visible_vertexBC
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the DirichletSet defined on the vertex)
*/
ErrorCode GetPointerToVertexBC(iMOAB_AppID pid, int vertex_BC_length, iMOAB_GlobalID* global_vertext_ID, int* num_vertex_BC, int* boundary_condition_value);

/**
  \fn ErrorCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int tag_storage_name_length)
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
ErrorCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int tag_storage_name_length);

/**
   \fn ErrorCode SetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to replace the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrorCode SetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrorCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrorCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrorCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrorCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrorCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)  The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)     The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[out] tag_storage_data (double*)       The array data of type <I>double</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (int)    The length of the tag_storage_name string
*/
ErrorCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length);

/**
   \fn ErrorCode GetNeighborElements(iMOAB_AppID pid, iMOAB_GlobalID global_element_ID, int* num_adjacent_elements, iMOAB_GlobalID* adjacent_element_IDs)
   \brief Compute the adjacencies for the element entities

   <B>Operations:</B> Collective

   \bug 1) How do we decide the elements/vertices to query if we are not passing EntityHandles ?
   \bug 2) We need a dimension for up/down adjacencies. Or do we care only about 1-ring neighbors ?

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  global_element_ID (iMOAB_GlobalID)     The global element ID for which adjacency information is needed
   \param[out] num_adjacent_elements (int*)           The total number of adjacent elements
   \param[out] adjacent_element_IDs (iMOAB_GlobalID*) The global element IDs of all adjacent elements to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrorCode GetNeighborElements(iMOAB_AppID pid, iMOAB_GlobalID global_element_ID, int* num_adjacent_elements, iMOAB_GlobalID* adjacent_element_IDs);

/**
   \fn ErrorCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_GlobalID global_vertex_ID, int* num_adjacent_vertices, iMOAB_GlobalID* adjacent_vertex_IDs)
   \brief Compute the adjacencies for the vertex entities

   <B>Operations:</B> Collective

   \bug 1) How do we decide the elements/vertices to query if we are not passing EntityHandles ?
   \bug 2) We need a dimension for up/down adjacencies. Or do we care only about 1-ring neighbors ?

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  global_vertex_ID (iMOAB_GlobalID)      The global vertex ID for which adjacency information is needed
   \param[out] num_adjacent_vertices (int*)           The total number of adjacent vertices
   \param[out] adjacent_vertex_IDs (iMOAB_GlobalID*)  The global element IDs of all adjacent vertices to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrorCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_GlobalID global_vertex_ID, int* num_adjacent_vertices, iMOAB_GlobalID* adjacent_vertex_IDs);

#endif
