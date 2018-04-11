/** \file iMOAB.cpp
*/

#include "moab/MOABConfig.h"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MOAB_HAVE_MPI
    #include "moab_mpi.h"
    #include "moab/ParallelComm.hpp"
    #include "moab/ParCommGraph.hpp"
#endif

#include <assert.h>
#include "moab/iMOAB.h"

/*
 this is needed so far because of direct access to hdf5/mhdf
  */
#ifdef MOAB_HAVE_HDF5
    #include "mhdf.h"
    #include <H5Tpublic.h>
#endif

#include "MBTagConventions.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/MergeMesh.hpp"

#ifdef MOAB_HAVE_TEMPESTREMAP
    #include "moab/IntxMesh/IntxUtils.hpp"

    #include "moab/Remapping/TempestRemapper.hpp"
    #include "moab/Remapping/TempestOfflineMap.hpp"
#endif

// C++ includes
#include <stdio.h>
#include <sstream>
#include <iostream>

// global variables ; should they be organized in a structure, for easier references?
// or how do we keep them global?

#ifdef __cplusplus
extern "C" {
#endif

#define CHKERRVAL(ierr) { if ( moab::MB_SUCCESS != ierr ) return 1; }
#define CHKIERRVAL(ierr) { if ( 0 != ierr ) return 1; }

struct appData
{
    EntityHandle file_set;
    int   external_id;  // external component id, unique for application
    Range all_verts;
    Range local_verts; // it could include shared, but not owned at the interface
    // these vertices would be all_verts if no ghosting was required
    Range ghost_vertices; // locally ghosted from other processors
    Range primary_elems;
    Range owned_elems;
    Range ghost_elems;
    int dimension; // 2 or 3, dimension of primary elements (redundant?)
    long num_global_elements; // reunion of all elements in primary_elements; either from hdf5 reading or from reduce
    long num_global_vertices; // reunion of all nodes, after sharing is resolved; it could be determined from hdf5 reading
    Range mat_sets;
    std::map<int, int> matIndex; // map from global block id to index in mat_sets
    Range neu_sets;
    Range diri_sets;
    std::map< std::string, Tag> tagMap;
    std::vector<Tag>  tagList;

#ifdef MOAB_HAVE_MPI
    std::vector<ParCommGraph*> pgraph; // created in order of other applications that communicate with this one
    // constructor for this ParCommGraph takes the joint comm and the MPI groups for each application
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
	moab::TempestRemapper* remapper;
    moab::TempestOfflineMap* weightMap;
	iMOAB_AppID pid_src;
	iMOAB_AppID pid_dest;
#endif
};

struct GlobalContext
{
    // are there reasons to have multiple moab inits? Is ref count needed?
    Interface* MBI;
    // we should also have the default tags stored, initialized
    Tag material_tag, neumann_tag, dirichlet_tag, globalID_tag; // material, neumann, dirichlet,  globalID
    int refCountMB;
    int iArgc;
    iMOAB_String* iArgv;
    int unused_pid;

    std::map<std::string, int> appIdMap;     // from app string (uppercase) to app id
    std::map<int, int> appIdCompMap;         // from component id to app id

#ifdef MOAB_HAVE_MPI
    std::vector<ParallelComm*> pcomms; // created in order of applications, one moab::ParallelComm for each
#endif

    std::vector<appData> appDatas; // the same order as pcomms

    GlobalContext() {MBI = 0; refCountMB = 0; unused_pid = 0; }
}  ;

static struct GlobalContext context;


ErrCode iMOAB_Initialize ( int argc, iMOAB_String* argv )
{
    context.iArgc = argc;
    context.iArgv = argv; // shallow copy

    if ( 0 == context.refCountMB )
    {
        context.MBI = new ( std::nothrow ) moab::Core;
        // retrieve the default tags
        const char* const shared_set_tag_names[] = {MATERIAL_SET_TAG_NAME,
                                                    NEUMANN_SET_TAG_NAME,
                                                    DIRICHLET_SET_TAG_NAME,
                                                    GLOBAL_ID_TAG_NAME
                                                   };
        // blocks, visible surfaceBC(neumann), vertexBC (Dirichlet), global id, parallel partition
        Tag gtags[4];

        for ( int i = 0; i < 4; i++ )
        {

            ErrorCode rval = context.MBI->tag_get_handle ( shared_set_tag_names[i], 1, MB_TYPE_INTEGER,
                             gtags[i], MB_TAG_ANY );
            CHKERRVAL(rval);
        }

        context.material_tag = gtags[0];
        context.neumann_tag = gtags[1];
        context.dirichlet_tag = gtags[2];
        context.globalID_tag = gtags[3];
    }

    context.refCountMB++;
    return 0;
}

ErrCode iMOAB_InitializeFortran()
{
    return iMOAB_Initialize ( 0, 0 );
}

ErrCode iMOAB_Finalize()
{
    context.refCountMB--;

    if ( 0 == context.refCountMB )
    { delete context.MBI; }

    return MB_SUCCESS;
}

ErrCode iMOAB_RegisterApplication ( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
                                    MPI_Comm* comm,
#endif
                                    int* compid,
                                    iMOAB_AppID pid )
{
    // will create a parallel comm for this application too, so there will be a
    // mapping from *pid to file set and to parallel comm instances
    std::string name ( app_name );

    if ( context.appIdMap.find ( name ) != context.appIdMap.end() )
    {
        std::cout << " application " << name << " already registered \n";
        return 1;
    }

    *pid =  context.unused_pid++;
    context.appIdMap[name] = *pid;

	std::cout << " application " << name << " with ID = " << *pid << " is registered now \n";

    if ( *compid <= 0 )
    {
        std::cout << " convention for external application is to have its id positive \n";
        return 1;
    }

    if ( context.appIdCompMap.find ( *compid ) != context.appIdCompMap.end() )
    {
        std::cout << " external application with comp id " << *compid << " is already registered\n";
        return 1;
    }

    context.appIdCompMap[*compid] = *pid;

    // now create ParallelComm and a file set for this application
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = new ParallelComm ( context.MBI, *comm );

#ifndef NDEBUG
    int index = pco->get_id(); // it could be useful to get app id from pcomm instance ...
    assert ( index == *pid );
    // here, we assert the the pid is the same as the id of the ParallelComm instance
    // useful for writing in parallel
#endif
    context.pcomms.push_back ( pco );
#endif

    // create now the file set that will be used for loading the model in
    EntityHandle file_set;
    ErrorCode rval = context.MBI->create_meshset ( MESHSET_SET, file_set );
    CHKERRVAL(rval);

    appData app_data;
    app_data.file_set = file_set;
    app_data.external_id = * compid; // will be used mostly for par comm graph

#ifdef MOAB_HAVE_TEMPESTREMAP
	app_data.remapper = NULL; // Only allocate as needed
#endif

    context.appDatas.push_back ( app_data ); // it will correspond to app_FileSets[*pid] will be the file set of interest
    return 0;
}

ErrCode iMOAB_RegisterFortranApplication ( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
        int* comm,
#endif
        int* compid, iMOAB_AppID pid, int app_name_length )
{
    std::string name ( app_name );

    if ( ( int ) strlen ( app_name ) > app_name_length )
    {
        std::cout << " length of string issue \n";
        return 1;
    }

    if ( context.appIdMap.find ( name ) != context.appIdMap.end() )
    {
        std::cout << " application " << name << " already registered \n";
        return 1;
    }

    *pid =  context.unused_pid++;
    context.appIdMap[name] = *pid;

    if ( context.appIdCompMap.find ( *compid ) != context.appIdCompMap.end() )
    {
        std::cout << " external application with comp id " << *compid << " is already registered\n";
        return 1;
    }

    context.appIdCompMap[*compid] = *pid;

#ifdef MOAB_HAVE_MPI
    // now create ParallelComm and a file set for this application
    // convert from fortran communicator to a c communicator
    // see transfer of handles
    // http://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node361.htm
    MPI_Comm ccomm = MPI_Comm_f2c ( ( MPI_Fint ) * comm );
    ParallelComm* pco = new ParallelComm ( context.MBI, ccomm );

#ifndef NDEBUG
    int index = pco->get_id(); // it could be useful to get app id from pcomm instance ...
    assert ( index == *pid );
#endif
    context.pcomms.push_back ( pco );
#endif

    // create now the file set that will be used for loading the model in
    EntityHandle file_set;
    ErrorCode rval = context.MBI->create_meshset ( MESHSET_SET, file_set );CHKERRVAL(rval);

    appData app_data;
    app_data.file_set = file_set;
    app_data.external_id = * compid; // will be used mostly for par comm graph
#ifdef MOAB_HAVE_TEMPESTREMAP
	app_data.remapper = NULL; // Only allocate as needed
#endif

    context.appDatas.push_back ( app_data ); // it will correspond to app_FileSets[*pid] will be the file set of interest
    return 0;
}

ErrCode iMOAB_DeregisterApplication ( iMOAB_AppID pid )
{
    // the file set , parallel comm are all in vectors indexed by *pid
    // assume we did not delete anything yet
    // *pid will not be reused if we register another application

    EntityHandle fileSet = context.appDatas[*pid].file_set;
    // get all entities part of the file set
    Range fileents;
    ErrorCode rval = context.MBI->get_entities_by_handle ( fileSet, fileents, /*recursive */true );
	CHKERRVAL(rval);

    fileents.insert ( fileSet );

    rval = context.MBI->get_entities_by_type ( fileSet, MBENTITYSET, fileents ); // append all mesh sets
	CHKERRVAL(rval);

#ifdef MOAB_HAVE_TEMPESTREMAP
  if (context.appDatas[*pid].remapper) delete context.appDatas[*pid].remapper;
#endif

#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.pcomms[*pid];
    // we could get the pco also with
    // ParallelComm * pcomm = ParallelComm::get_pcomm(context.MBI, *pid);
    delete pco;
    std::vector<ParCommGraph*>& pargs = context.appDatas[*pid].pgraph;

    // free the parallel comm graphs associated with this app
    for ( size_t k = 0; k < pargs.size(); k++ )
    { delete pargs[k]; }

#endif

    // delete first all except vertices
    Range vertices = fileents.subset_by_type ( MBVERTEX );
    Range noverts = subtract ( fileents, vertices );

    rval = context.MBI->delete_entities ( noverts );CHKERRVAL(rval);
    // now retrieve connected elements that still exist (maybe in other sets, pids?)
    Range adj_ents_left;
    rval = context.MBI->get_adjacencies(vertices, 1, false, adj_ents_left, Interface::UNION);CHKERRVAL(rval);
    rval = context.MBI->get_adjacencies(vertices, 2, false, adj_ents_left, Interface::UNION);CHKERRVAL(rval);
    rval = context.MBI->get_adjacencies(vertices, 3, false, adj_ents_left, Interface::UNION);CHKERRVAL(rval);

    if (!adj_ents_left.empty())
    {
      Range conn_verts;
      rval = context.MBI->get_connectivity(adj_ents_left, conn_verts);CHKERRVAL(rval);
      vertices = subtract(vertices, conn_verts);
    }

    rval = context.MBI->delete_entities ( vertices );CHKERRVAL(rval);

    std::map<std::string, int>::iterator mit;

    for ( mit = context.appIdMap.begin(); mit != context.appIdMap.end(); mit++ )
    {
        int pidx = mit->second;

        if ( *pid == pidx )
        {
            break;
        }
    }

    context.appIdMap.erase ( mit );
    std::map<int, int>::iterator mit1;

    for ( mit1 = context.appIdCompMap.begin(); mit1 != context.appIdCompMap.end(); mit1++ )
    {
        int pidx = mit1->second;

        if ( *pid == pidx )
        {
            break;
        }
    }

    context.appIdCompMap.erase ( mit1 );

    context.unused_pid--;
    context.appDatas.pop_back();
#ifdef MOAB_HAVE_MPI
    context.pcomms.pop_back();
#endif
    return 0;
}

ErrCode iMOAB_ReadHeaderInfo ( const iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
{
#ifdef MOAB_HAVE_HDF5
    std::string filen ( filename );

    if ( filename_length < ( int ) strlen ( filename ) )
    {
        filen = filen.substr ( 0, filename_length );
    }

    *num_global_vertices = 0;
    int edges = 0;
    int faces = 0;
    int regions = 0;
    *num_global_elements = 0;
    *num_dimension = 0;
    *num_parts = 0;

    mhdf_FileHandle file;
    mhdf_Status status;
    unsigned long max_id;
    struct mhdf_FileDesc* data;

    file = mhdf_openFile ( filen.c_str(), 0, &max_id, -1, &status );

    if ( mhdf_isError ( &status ) )
    {
        fprintf ( stderr, "%s: %s\n", filename, mhdf_message ( &status ) );
        return 1;
    }

    data = mhdf_getFileSummary ( file, H5T_NATIVE_ULONG, &status, 1 ); // will use extra set info; will get parallel partition tag info too!

    if ( mhdf_isError ( &status ) )
    {
        fprintf ( stderr, "%s: %s\n", filename, mhdf_message ( &status ) );
        return 1;
    }

    *num_dimension = data->nodes.vals_per_ent;
    *num_global_vertices = ( int ) data->nodes.count;

    for ( int i = 0; i < data->num_elem_desc; i++ )
    {
        struct mhdf_ElemDesc* el_desc = & ( data->elems[i] );
        struct mhdf_EntDesc* ent_d = & ( el_desc->desc );

        if ( 0 == strcmp ( el_desc->type, mhdf_EDGE_TYPE_NAME ) ) { edges += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_TRI_TYPE_NAME ) )  { faces += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_QUAD_TYPE_NAME ) ) { faces += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_POLYGON_TYPE_NAME ) ) { faces += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_TET_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_PYRAMID_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_PRISM_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mdhf_KNIFE_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mdhf_HEX_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_POLYHEDRON_TYPE_NAME ) ) { regions += ent_d->count; }

        if ( 0 == strcmp ( el_desc->type, mhdf_SEPTAHEDRON_TYPE_NAME ) ) { regions += ent_d->count; }
    }

    *num_parts = data->numEntSets[0];

    // is this required?
    if ( edges > 0 )
    {
        *num_dimension = 1; // I don't think it will ever return 1
        *num_global_elements = edges;
    }

    if ( faces > 0 )
    {
        *num_dimension = 2;
        *num_global_elements = faces;
    }

    if ( regions > 0 )
    {
        *num_dimension = 3;
        *num_global_elements = regions;
    }

    mhdf_closeFile ( file, &status );

    free ( data );

#else
    std::cout << " cannot retrieve header info except for h5m file \n";
#endif

    return 0;
}


ErrCode iMOAB_LoadMesh ( iMOAB_AppID pid, const iMOAB_String filename, const iMOAB_String read_options, int* num_ghost_layers, int filename_length, int read_options_length )
{
    if ( ( int ) strlen ( filename ) > filename_length )
    {
        std::cout << " filename length issue\n";
        return 1;
    }

    if ( ( int ) strlen ( read_options ) > read_options_length )
    {
        std::cout << " read options length issue\n";
        return 1;
    }

    // make sure we use the file set and pcomm associated with the *pid
    std::ostringstream newopts;
    newopts  << read_options;
#ifdef MOAB_HAVE_MPI
    std::string opts ( read_options );
    std::string pcid ( "PARALLEL_COMM=" );
    std::size_t found = opts.find ( pcid );

    if ( found != std::string::npos )
    {
        std::cerr << " cannot specify PARALLEL_COMM option, it is implicit \n";
        return 1;
    }

    // in serial, apply PARALLEL_COMM option only for h5m files; it does not work for .g files (used in test_remapping)
    std::string filen(filename);
    std::string::size_type idx = filen.rfind('.');

    if(idx != std::string::npos)
    {
      std::string extension = filen.substr(idx+1);
      if (extension == std::string("h5m"))
          newopts << ";;PARALLEL_COMM=" << *pid;
    }


    if ( *num_ghost_layers >= 1 )
    {
      // if we want ghosts, we will want additional entities, the last .1
      // because the addl ents can be edges, faces that are part of the neumann sets
      std::string pcid2 ( "PARALLEL_GHOSTS=" );
      std::size_t found2 = opts.find ( pcid2 );

      if ( found2 != std::string::npos )
      {
          std::cout << " PARALLEL_GHOSTS option is already specified, ignore passed number of layers \n";
      }
      else
      {
        // dimension of primary entities is 3 here, but it could be 2 for climate meshes; we would need to pass
        // PARALLEL_GHOSTS explicitly for 2d meshes, for example:  ";PARALLEL_GHOSTS=2.0.1"
        newopts << ";PARALLEL_GHOSTS=3.0." << *num_ghost_layers << ".3";
      }
    }

#else
    *num_ghost_layers = 0; // do not use in case of serial run
#endif
    ErrorCode rval = context.MBI->load_file ( filename, &context.appDatas[*pid].file_set, newopts.str().c_str() );
    CHKERRVAL(rval);

#ifdef VERBOSE

    // some debugging stuff
    std::ostringstream outfile;
#ifdef MOAB_HAVE_MPI
    int rank = context.pcomms[*pid]->rank();
    int nprocs = context.pcomms[*pid]->size();
    outfile << "TaskMesh_n" << nprocs << "." << rank << ".h5m";
#else
    outfile << "TaskMesh_n1.0.h5m";
#endif
    // the mesh contains ghosts too, but they are not part of mat/neumann set
    // write in serial the file, to see what tags are missing
    rval = context.MBI->write_file ( outfile.str().c_str() ); // everything on current task, written in serial
    CHKERRVAL(rval);

#endif
    int rc = iMOAB_UpdateMeshInfo ( pid );
    return rc;
}


ErrCode iMOAB_WriteMesh ( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
{
    // maybe do some processing of strings and lengths
    if ( ( int ) strlen ( filename ) > filename_length )
    {
        std::cout << " file name length issue\n";
        return 1;
    }

    if ( ( int ) strlen ( write_options ) > write_options_length )
    {
        std::cout << " write options issue\n";
        return 1;
    }

    std::ostringstream newopts;
#ifdef MOAB_HAVE_MPI
    std::string write_opts ( write_options );
    std::string pcid ( "PARALLEL_COMM=" );
    std::size_t found = write_opts.find ( pcid );

    if ( found != std::string::npos )
    {
        std::cerr << " cannot specify PARALLEL_COMM option, it is implicit \n";
        return 1;
    }

    // if write in parallel, add pc option, to be sure about which ParallelComm instance is used
    std::string pw ( "PARALLEL=WRITE_PART" );
    found = write_opts.find ( pw );

    if ( found != std::string::npos )
    {
        newopts << "PARALLEL_COMM=" << *pid << ";";
    }

#endif
    newopts  << write_options;
    ErrorCode rval = context.MBI->write_file ( filename, 0, newopts.str().c_str(),  &context.appDatas[*pid].file_set, 1 );

    if ( MB_SUCCESS != rval )
    { return 1; }

    return 0;
}


ErrCode iMOAB_UpdateMeshInfo ( iMOAB_AppID pid )
{
    // this will include ghost elements info
    //
    appData& data = context.appDatas[*pid];
    EntityHandle fileSet = data.file_set;
    // first clear all data ranges; this can be called after ghosting
    data.all_verts.clear();
    data.primary_elems.clear();
    data.local_verts.clear();
    data.ghost_vertices.clear();
    data.owned_elems.clear();
    data.ghost_elems.clear();
    data.mat_sets.clear();
    data.neu_sets.clear();
    data.diri_sets.clear();

    ErrorCode rval = context.MBI->get_entities_by_type ( fileSet, MBVERTEX, data.all_verts, true ); // recursive
	CHKERRVAL(rval);

    rval = context.MBI->get_entities_by_dimension ( fileSet, 3, data.primary_elems, true ); // recursive
	CHKERRVAL(rval);

    data.dimension = 3;

    if ( data.primary_elems.empty() )
    {
        context.appDatas[*pid].dimension = 2;
        rval = context.MBI->get_entities_by_dimension ( fileSet, 2, data.primary_elems, true ); // recursive
        CHKERRVAL(rval);

        if ( data.primary_elems.empty() )
        {
            context.appDatas[*pid].dimension = 1;
            rval = context.MBI->get_entities_by_dimension ( fileSet, 1, data.primary_elems, true ); // recursive
            CHKERRVAL(rval);

            if ( data.primary_elems.empty() )
            { return 1; } // no elements of dimension 1 or 2 or 3
        }
    }

#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.pcomms[*pid];

    // filter ghost vertices, from local
    rval = pco -> filter_pstatus ( data.all_verts, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.local_verts );
    CHKERRVAL(rval);

    data.ghost_vertices = subtract ( data.all_verts, data.local_verts );

    // get all blocks, BCs, etc

    // filter ghost elements, from local
    rval = pco -> filter_pstatus ( data.primary_elems, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.owned_elems );
	CHKERRVAL(rval);

    data.ghost_elems = subtract ( data.primary_elems, data.owned_elems );

#else

    data.local_verts = data.all_verts;
    data.owned_elems = data.primary_elems;

#endif
    rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.material_tag ), 0, 1, data.mat_sets, Interface::UNION );
	CHKERRVAL(rval);

    rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.neumann_tag ), 0, 1, data.neu_sets, Interface::UNION );
	CHKERRVAL(rval);

    rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.dirichlet_tag ), 0, 1, data.diri_sets, Interface::UNION );
	CHKERRVAL(rval);

    return 0;
}


ErrCode iMOAB_GetMeshInfo ( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int* num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
{
	ErrorCode rval;
    // this will include ghost elements
    appData& data = context.appDatas[*pid];
    EntityHandle fileSet = data.file_set;
    // first clear all data ranges; this can be called after ghosting

	if (num_visible_elements) {
		num_visible_elements[2] = ( int ) data.primary_elems.size();
		// separate ghost and local/owned primary elements
		num_visible_elements[0] = ( int ) data.owned_elems.size();
		num_visible_elements[1] = ( int ) data.ghost_elems.size();
    }
    if (num_visible_vertices) {
		num_visible_vertices[2] = ( int ) data.all_verts.size();
		num_visible_vertices[1] = ( int ) data.ghost_vertices.size();
		num_visible_vertices[0] =  num_visible_vertices[2] - num_visible_vertices[1]; // local are those that are not ghosts
	}

	if (num_visible_blocks) {	
		rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.material_tag ), 0, 1, data.mat_sets, Interface::UNION );
		CHKERRVAL(rval);

		num_visible_blocks[2] = data.mat_sets.size();
		num_visible_blocks[0] = num_visible_blocks[2];
		num_visible_blocks[1] = 0;
    }

	if (num_visible_surfaceBC) {
		rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.neumann_tag ), 0, 1, data.neu_sets, Interface::UNION );
		CHKERRVAL(rval);

		num_visible_surfaceBC[2] = 0;
		// count how many faces are in each neu set, and how many regions are
		// adjacent to them;
		int numNeuSets = ( int ) data.neu_sets.size();

		for ( int i = 0; i < numNeuSets; i++ )
		{
			Range subents;
			EntityHandle nset = data.neu_sets[i];
			rval = context.MBI->get_entities_by_dimension ( nset, data.dimension - 1, subents );
			CHKERRVAL(rval);

			for ( Range::iterator it = subents.begin(); it != subents.end(); ++it )
			{
				EntityHandle subent = *it;
				Range adjPrimaryEnts;
				rval = context.MBI->get_adjacencies ( &subent, 1, data.dimension, false, adjPrimaryEnts );
				CHKERRVAL(rval);

				num_visible_surfaceBC[2] += ( int ) adjPrimaryEnts.size();
			}
		}

		num_visible_surfaceBC[0] = num_visible_surfaceBC[2];
		num_visible_surfaceBC[1] = 0; //
    }

	if (num_visible_vertexBC) {
		rval = context.MBI->get_entities_by_type_and_tag ( fileSet, MBENTITYSET, & ( context.dirichlet_tag ), 0, 1, data.diri_sets, Interface::UNION );
		CHKERRVAL(rval);

		num_visible_vertexBC[2] = 0;
		int numDiriSets = ( int ) data.diri_sets.size();

		for ( int i = 0; i < numDiriSets; i++ )
		{
			Range verts;
			EntityHandle diset = data.diri_sets[i];
			rval = context.MBI->get_entities_by_dimension ( diset, 0, verts );
			CHKERRVAL(rval);

			num_visible_vertexBC[2] += ( int ) verts.size();
		}

		num_visible_vertexBC[0] = num_visible_vertexBC[2];
		num_visible_vertexBC[1] = 0;
    }

    return 0;
}

ErrCode iMOAB_GetVertexID ( iMOAB_AppID pid, int* vertices_length, iMOAB_GlobalID* global_vertex_ID )
{
    //
    Range& verts = context.appDatas[*pid].all_verts;

    if ( ( int ) verts.size() != *vertices_length )
    { return 1; } // problem with array length

    // global id tag is context.globalID_tag
    ErrorCode rval = context.MBI->tag_get_data ( context.globalID_tag, verts, global_vertex_ID );
	CHKERRVAL(rval);

    return 0;
}

ErrCode iMOAB_GetVertexOwnership ( iMOAB_AppID pid, int* vertices_length, int* visible_global_rank_ID )
{
    Range& verts = context.appDatas[*pid].all_verts;
    int i = 0;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.pcomms[*pid];

    for ( Range::iterator vit = verts.begin(); vit != verts.end(); vit++, i++ )
    {
        ErrorCode rval = pco->  get_owner ( *vit, visible_global_rank_ID[i] );
        CHKERRVAL(rval);
    }

    if ( i != *vertices_length )
    { return 1; } // warning array allocation problem

#else

    /* everything owned by proc 0 */
    if ( ( int ) verts.size() != *vertices_length )
    { return 1; } // warning array allocation problem

    for ( Range::iterator vit = verts.begin(); vit != verts.end(); vit++, i++ )
    { visible_global_rank_ID[i] = 0; } // all vertices are owned by processor 0, as this is serial run

#endif
    return 0;
}

ErrCode iMOAB_GetVisibleVerticesCoordinates ( iMOAB_AppID pid, int* coords_length, double* coordinates )
{
    Range& verts = context.appDatas[*pid].all_verts;

    // interleaved coordinates, so that means deep copy anyway
    if ( *coords_length != 3 * ( int ) verts.size() )
    { return 1; }

    ErrorCode rval = context.MBI->get_coords ( verts, coordinates );
    CHKERRVAL(rval);

    return 0;
}

ErrCode iMOAB_GetBlockID ( iMOAB_AppID pid, int* block_length, iMOAB_GlobalID* global_block_IDs )
{
    // local id blocks? they are counted from 0 to number of visible blocks ...
    // will actually return material set tag value for global
    Range& matSets = context.appDatas[*pid].mat_sets;

    if ( *block_length != ( int ) matSets.size() )
    { return 1; }

    // return material set tag gtags[0 is material set tag
    ErrorCode rval = context.MBI->tag_get_data ( context.material_tag, matSets, global_block_IDs );
    CHKERRVAL(rval);

    // populate map with index
    std::map <int, int>& matIdx = context.appDatas[*pid].matIndex;

    //
    for ( int i = 0; i < ( int ) matSets.size(); i++ )
    {
        matIdx[global_block_IDs[i]] = i;
    }

    return 0;
}

ErrCode iMOAB_GetBlockInfo ( iMOAB_AppID pid, iMOAB_GlobalID* global_block_ID,
                             int* vertices_per_element, int* num_elements_in_block )
{
    std::map<int, int>& matMap = context.appDatas[*pid].matIndex;
    std::map<int, int>::iterator it = matMap.find ( *global_block_ID );

    if ( it == matMap.end() )
    { return 1; } // error in finding block with id

    int blockIndex = matMap[*global_block_ID];
    EntityHandle matMeshSet = context.appDatas[*pid].mat_sets[blockIndex];
    Range blo_elems;
    ErrorCode rval = context.MBI-> get_entities_by_handle ( matMeshSet, blo_elems );

    if ( MB_SUCCESS != rval ||  blo_elems.empty() )
    { return 1; }

    EntityType type = context.MBI->type_from_handle ( blo_elems[0] );

    if ( !blo_elems.all_of_type ( type ) )
    { return 1; } //not all of same  type

    const EntityHandle* conn = NULL;
    int num_verts = 0;
    rval = context.MBI->get_connectivity ( blo_elems[0], conn, num_verts );
	CHKERRVAL(rval);

    *vertices_per_element = num_verts;
    *num_elements_in_block = ( int ) blo_elems.size();

    return 0;
}

ErrCode iMOAB_GetVisibleElementsInfo ( iMOAB_AppID pid, int* num_visible_elements,
                                       iMOAB_GlobalID* element_global_IDs, int* ranks, iMOAB_GlobalID* block_IDs )
{
    appData& data =  context.appDatas[*pid];
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.pcomms[*pid];
#endif

    ErrorCode rval = context.MBI-> tag_get_data ( context.globalID_tag, data.primary_elems, element_global_IDs );
	CHKERRVAL(rval);

    int i = 0;

    for ( Range::iterator eit = data.primary_elems.begin(); eit != data.primary_elems.end(); ++eit, ++i )
    {
#ifdef MOAB_HAVE_MPI
        rval = pco->get_owner ( *eit, ranks[i] );
        CHKERRVAL(rval);

#else
        /* everything owned by task 0 */
        ranks[i] = 0;
#endif
    }

    for ( Range::iterator mit = data.mat_sets.begin(); mit != data.mat_sets.end(); ++mit )
    {
        EntityHandle matMeshSet = *mit;
        Range elems;
        rval = context.MBI-> get_entities_by_handle ( matMeshSet, elems );
        CHKERRVAL(rval);

        int valMatTag;
        rval = context.MBI->tag_get_data ( context.material_tag, &matMeshSet, 1, &valMatTag );
        CHKERRVAL(rval);

        for ( Range::iterator eit = elems.begin(); eit != elems.end(); ++eit )
        {
            EntityHandle eh = *eit;
            int index = data.primary_elems.index ( eh );

            if ( -1 == index )
            { return 1; }

            if ( -1 >= *num_visible_elements )
            { return 1; }

            block_IDs[index] = valMatTag;
        }
    }


    return 0;
}

ErrCode iMOAB_GetBlockElementConnectivities ( iMOAB_AppID pid, iMOAB_GlobalID* global_block_ID, int* connectivity_length, int* element_connectivity )
{
    appData& data =  context.appDatas[*pid];
    std::map<int, int>& matMap = data.matIndex;
    std::map<int, int>::iterator it = matMap.find ( *global_block_ID );

    if ( it == matMap.end() )
    { return 1; } // error in finding block with id

    int blockIndex = matMap[*global_block_ID];
    EntityHandle matMeshSet = data.mat_sets[blockIndex];
    std::vector<EntityHandle> elems;

    ErrorCode rval = context.MBI-> get_entities_by_handle ( matMeshSet, elems );

    if ( MB_SUCCESS != rval ||  elems.empty() )
    { return 1; }

    std::vector<EntityHandle> vconnect;
    rval = context.MBI->get_connectivity ( &elems[0], elems.size(), vconnect );
    CHKERRVAL(rval);

    if ( *connectivity_length != ( int ) vconnect.size() )
    { return 1; } // mismatched sizes

    for ( int i = 0; i < *connectivity_length; i++ )
    {
        int inx = data.all_verts.index ( vconnect[i] );

        if ( -1 == inx )
        { return 1; } // error, vertex not in local range

        element_connectivity[i] = inx;
    }

    return 0;
}

ErrCode iMOAB_GetElementConnectivity ( iMOAB_AppID pid, iMOAB_LocalID* elem_index, int* connectivity_length, int* element_connectivity )
{
    appData& data =  context.appDatas[*pid];
    assert ( ( *elem_index >= 0 )  && ( *elem_index < ( int ) data.primary_elems.size() ) );

    int num_nodes;
    const EntityHandle* conn;

    EntityHandle eh = data.primary_elems[*elem_index];

    ErrorCode rval = context.MBI->get_connectivity ( eh, conn, num_nodes );
    CHKERRVAL(rval);

    if ( * connectivity_length < num_nodes )
    { return 1; } // wrong number of vertices

    for ( int i = 0; i < num_nodes; i++ )
    {
        int index = data.all_verts.index ( conn[i] );

        if ( -1 == index )
        { return 1; }

        element_connectivity[i] = index;
    }

    * connectivity_length = num_nodes;
    return 0;
}

ErrCode iMOAB_GetElementOwnership ( iMOAB_AppID pid, iMOAB_GlobalID* global_block_ID, int* num_elements_in_block, int* element_ownership )
{
    std::map<int, int>& matMap = context.appDatas[*pid].matIndex;

    std::map<int, int>::iterator it = matMap.find ( *global_block_ID );

    if ( it == matMap.end() )
    { return 1; } // error in finding block with id

    int blockIndex = matMap[*global_block_ID];
    EntityHandle matMeshSet = context.appDatas[*pid].mat_sets[blockIndex];
    Range elems;

    ErrorCode rval = context.MBI-> get_entities_by_handle ( matMeshSet, elems );

    if ( MB_SUCCESS != rval ||  elems.empty() )
    { return 1; }

    if ( *num_elements_in_block != ( int ) elems.size() )
    { return 1; } // bad memory allocation

    int i = 0;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.pcomms[*pid];
#endif

    for ( Range::iterator vit = elems.begin(); vit != elems.end(); vit++, i++ )
    {
#ifdef MOAB_HAVE_MPI
        rval = pco->  get_owner ( *vit, element_ownership[i] );
        CHKERRVAL(rval);
#else
        element_ownership[i] = 0; /* owned by 0 */
#endif
    }

    return 0;
}

ErrCode iMOAB_GetElementID ( iMOAB_AppID pid, iMOAB_GlobalID* global_block_ID, int* num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID )
{
    appData& data = context.appDatas[*pid];
    std::map<int, int>& matMap = data.matIndex;

    std::map<int, int>::iterator it = matMap.find ( *global_block_ID );

    if ( it == matMap.end() )
    { return 1; } // error in finding block with id

    int blockIndex = matMap[*global_block_ID];
    EntityHandle matMeshSet = data.mat_sets[blockIndex];
    Range elems;
    ErrorCode rval = context.MBI-> get_entities_by_handle ( matMeshSet, elems );

    if ( MB_SUCCESS != rval ||  elems.empty() )
    { return 1; }

    if ( *num_elements_in_block != ( int ) elems.size() )
    { return 1; } // bad memory allocation

    rval = context.MBI->tag_get_data ( context.globalID_tag, elems, global_element_ID );

    if ( MB_SUCCESS != rval )
    { return 1; }

    // check that elems are among primary_elems in data
    for ( int i = 0; i < *num_elements_in_block; i++ )
    {
        local_element_ID[i] = data.primary_elems.index ( elems[i] );

        if ( -1 == local_element_ID[i] )
        { return 1; }// error, not in local primary elements
    }

    return 0;
}

ErrCode iMOAB_GetPointerToSurfaceBC ( iMOAB_AppID pid, int* surface_BC_length, iMOAB_LocalID* local_element_ID, int* reference_surface_ID, int* boundary_condition_value )
{
    // we have to fill bc data for neumann sets;/

    // it was counted above, in GetMeshInfo
    appData& data = context.appDatas[*pid];
    int numNeuSets = ( int ) data.neu_sets.size();

    int index = 0; // index [0, surface_BC_length) for the arrays returned

    for ( int i = 0; i < numNeuSets; i++ )
    {
        Range subents;
        EntityHandle nset = data.neu_sets[i];
        ErrorCode rval = context.MBI->get_entities_by_dimension ( nset, data.dimension - 1, subents );
        CHKERRVAL(rval);

        int neuVal ;
        rval = context.MBI->tag_get_data ( context.neumann_tag, &nset, 1, &neuVal );
        CHKERRVAL(rval);

        for ( Range::iterator it = subents.begin(); it != subents.end(); ++it )
        {
            EntityHandle subent = *it;
            Range adjPrimaryEnts;
            rval = context.MBI->get_adjacencies ( &subent, 1, data.dimension, false, adjPrimaryEnts );
            CHKERRVAL(rval);

            // get global id of the primary ents, and side number of the quad/subentity
            // this is moab ordering
            for ( Range::iterator pit = adjPrimaryEnts.begin(); pit != adjPrimaryEnts.end(); pit++ )
            {
                EntityHandle primaryEnt = *pit;
                // get global id
                /*int globalID;
                rval = context.MBI->tag_get_data(gtags[3], &primaryEnt, 1, &globalID);
                if (MB_SUCCESS!=rval)
                  return 1;
                global_element_ID[index] = globalID;*/
                // get local element id
                local_element_ID[index] = data.primary_elems.index ( primaryEnt );

                if ( -1 == local_element_ID[index] )
                { return 1; } // did not find the element locally

                int side_number, sense, offset;
                rval = context.MBI->side_number ( primaryEnt, subent,  side_number, sense, offset );
                CHKERRVAL(rval);

                reference_surface_ID[index] = side_number + 1; // moab is from 0 to 5, it needs 1 to 6
                boundary_condition_value[index] = neuVal;
                index++;
            }
        }
    }

    if ( index != *surface_BC_length )
    { return 1; } // error in array allocations

    return 0;
}

ErrCode iMOAB_GetPointerToVertexBC ( iMOAB_AppID pid, int* vertex_BC_length,
                                     iMOAB_LocalID* local_vertex_ID, int* boundary_condition_value )
{
    // it was counted above, in GetMeshInfo
    appData& data = context.appDatas[*pid];
    int numDiriSets = ( int ) data.diri_sets.size();
    int index = 0; // index [0, *vertex_BC_length) for the arrays returned

    for ( int i = 0; i < numDiriSets; i++ )
    {
        Range verts;
        EntityHandle diset = data.diri_sets[i];
        ErrorCode rval = context.MBI->get_entities_by_dimension ( diset, 0, verts );
        CHKERRVAL(rval);

        int diriVal;
        rval = context.MBI->tag_get_data ( context.dirichlet_tag, &diset, 1, &diriVal );
        CHKERRVAL(rval);

        for ( Range::iterator vit = verts.begin(); vit != verts.end(); ++vit )
        {
            EntityHandle vt = *vit;
            /*int vgid;
            rval = context.MBI->tag_get_data(gtags[3], &vt, 1, &vgid);
            if (MB_SUCCESS!=rval)
              return 1;
            global_vertext_ID[index] = vgid;*/
            local_vertex_ID[index] = data.all_verts.index ( vt );

            if ( -1 == local_vertex_ID[index] )
            { return 1; } // vertex was not found

            boundary_condition_value[index] = diriVal;
            index++;
        }
    }

    if ( *vertex_BC_length != index )
    { return 1; } // array allocation issue

    return 0;
}

ErrCode iMOAB_DefineTagStorage ( iMOAB_AppID pid, const iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int* tag_index,  int tag_storage_name_length )
{
    // see if the tag is already existing, and if yes, check the type, length
    if ( *tag_type < 0 || *tag_type > 5 )
    { return 1; } // we have 6 types of tags supported so far

    DataType tagDataType;
    TagType tagType;
    void* defaultVal = NULL;
    int* defInt = new int [*components_per_entity];
    double* defDouble = new double [*components_per_entity];
    EntityHandle* defHandle = new EntityHandle[*components_per_entity];

    for ( int i = 0; i < *components_per_entity; i++ )
    {
        defInt[i] = 0;
        defDouble[i] = 0.;
        defHandle[i] = ( EntityHandle ) 0;
    }

    switch ( *tag_type )
    {
        case 0: tagDataType = MB_TYPE_INTEGER; tagType = MB_TAG_DENSE; defaultVal = defInt; break;

        case 1: tagDataType = MB_TYPE_DOUBLE;  tagType = MB_TAG_DENSE; defaultVal = defDouble; break;

        case 2: tagDataType = MB_TYPE_HANDLE;  tagType = MB_TAG_DENSE; defaultVal = defHandle; break;

        case 3: tagDataType = MB_TYPE_INTEGER; tagType = MB_TAG_SPARSE; defaultVal = defInt; break;

        case 4: tagDataType = MB_TYPE_DOUBLE;  tagType = MB_TAG_SPARSE; defaultVal = defDouble; break;

        case 5: tagDataType = MB_TYPE_HANDLE;  tagType = MB_TAG_SPARSE; defaultVal = defHandle; break;

        default :
        {
            delete [] defInt;
            delete [] defDouble;
            delete [] defHandle;
            return 1;
        } // error
    }

    std::string tag_name ( tag_storage_name );

    if ( tag_storage_name_length < ( int ) strlen ( tag_storage_name ) )
    {
        tag_name = tag_name.substr ( 0, tag_storage_name_length );
    }

    Tag tagHandle;
    ErrorCode rval = context.MBI->tag_get_handle ( tag_name.c_str(), *components_per_entity,
                     tagDataType,
                     tagHandle, tagType, defaultVal );

    if ( MB_TAG_NOT_FOUND == rval )
    {    	
		rval = context.MBI->tag_get_handle ( tag_name.c_str(), *components_per_entity,
						 tagDataType,
						 tagHandle, tagType|MB_TAG_CREAT, defaultVal );
    }

    // we don't need default values anymore, avoid leaks
    delete [] defInt;
    delete [] defDouble;
    delete [] defHandle;

    appData& data = context.appDatas[*pid];

    if ( MB_ALREADY_ALLOCATED == rval )
    {
        std::map<std::string, Tag>& mTags = data.tagMap;
        std::map<std::string, Tag>::iterator mit = mTags.find ( tag_name );

        if ( mit == mTags.end() )
        {
            // add it to the map
            mTags[tag_name] = tagHandle;
            // push it to the list of tags, too
            *tag_index = ( int ) data.tagList.size();
            data.tagList.push_back ( tagHandle ) ;
        }

        return 0; // OK, we found it, and we have it stored in the map tag
    }
    else if ( MB_SUCCESS == rval )
    {
        data.tagMap[tag_name] = tagHandle;
        *tag_index = ( int ) data.tagList.size();
        data.tagList.push_back ( tagHandle ) ;
        return 0;
    }
    else 
	    return 1; // some error, maybe the tag was not created
}

ErrCode iMOAB_SetIntTagStorage ( iMOAB_AppID pid, const iMOAB_String tag_storage_name,
                                 int* num_tag_storage_length, int* ent_type, int* tag_storage_data,
                                 int tag_storage_name_length )
{
    std::string tag_name ( tag_storage_name );

    if ( tag_storage_name_length < ( int ) strlen ( tag_storage_name ) )
    {
        tag_name = tag_name.substr ( 0, tag_storage_name_length );
    }

    appData& data = context.appDatas[*pid];

    if ( data.tagMap.find ( tag_name ) == data.tagMap.end() )
    { return 1; } // tag not defined

    Tag tag =  data.tagMap[tag_name];

    int tagLength = 0;
    ErrorCode rval = context.MBI->tag_get_length ( tag, tagLength );
    CHKERRVAL(rval);

    DataType  dtype;
    rval = context.MBI->tag_get_data_type ( tag, dtype );

    if ( MB_SUCCESS != rval || dtype != MB_TYPE_INTEGER )
    { return 1; }

    // set it on a subset of entities, based on type and length
    Range* ents_to_set;

    if ( *ent_type == 0 ) // vertices
    { ents_to_set = &data.all_verts; }
    else  // if (*ent_type == 1) // *ent_type can be 0 (vertices) or 1 (elements)
    { ents_to_set = &data.primary_elems; }

    int nents_to_be_set = *num_tag_storage_length / tagLength;

    if ( nents_to_be_set > ( int ) ents_to_set->size() || nents_to_be_set < 1 )
    { return 1; } // to many entities to be set or too few

    // restrict the range; everything is contiguous; or not?

    Range contig_range ( * ( ents_to_set->begin() ), * ( ents_to_set->begin() + nents_to_be_set - 1 ) );
    rval = context.MBI->tag_set_data ( tag, contig_range, tag_storage_data );
    CHKERRVAL(rval);

    return 0; // no error
}

ErrCode iMOAB_GetIntTagStorage ( iMOAB_AppID pid, const iMOAB_String tag_storage_name, int* num_tag_storage_length, int* ent_type, int* tag_storage_data, int tag_storage_name_length )
{
    std::string tag_name ( tag_storage_name );

    if ( tag_storage_name_length < ( int ) tag_name.length() )
    {
        tag_name = tag_name.substr ( 0, tag_storage_name_length );
    }

    appData& data = context.appDatas[*pid];

    if ( data.tagMap.find ( tag_name ) == data.tagMap.end() )
    { return 1; } // tag not defined

    Tag tag =  data.tagMap[tag_name];

    int tagLength = 0;
    ErrorCode rval = context.MBI->tag_get_length ( tag, tagLength );
    CHKERRVAL(rval);

    DataType  dtype;
    rval = context.MBI->tag_get_data_type ( tag, dtype );

    if ( MB_SUCCESS != rval || dtype != MB_TYPE_INTEGER )
    { return 1; }

    // set it on a subset of entities, based on type and length
    Range* ents_to_get;

    if ( *ent_type == 0 ) // vertices
    { ents_to_get = &data.all_verts; }
    else  // if (*ent_type == 1)
    { ents_to_get = &data.primary_elems; }

    int nents_to_get = *num_tag_storage_length / tagLength;

    if ( nents_to_get > ( int ) ents_to_get->size() || nents_to_get < 1 )
    { return 1; } // to many entities to get, or too little

    // restrict the range; everything is contiguous; or not?

    Range contig_range ( * ( ents_to_get->begin() ), * ( ents_to_get->begin() + nents_to_get - 1 ) );

    rval = context.MBI->tag_get_data ( tag, contig_range, tag_storage_data );
	CHKERRVAL(rval);

    return 0; // no error
}

ErrCode iMOAB_SetDoubleTagStorage ( iMOAB_AppID pid, const iMOAB_String tag_storage_name, int* num_tag_storage_length, int* ent_type, double* tag_storage_data, int tag_storage_name_length )
{
    // exactly the same code as for int tag :) maybe should check the type of tag too
    std::string tag_name ( tag_storage_name );

    if ( tag_storage_name_length < ( int ) tag_name.length() )
    {
        tag_name = tag_name.substr ( 0, tag_storage_name_length );
    }

    appData& data = context.appDatas[*pid];

    if ( data.tagMap.find ( tag_name ) == data.tagMap.end() )
    { return 1; } // tag not defined

    Tag tag =  data.tagMap[tag_name];

    int tagLength = 0;
    ErrorCode rval = context.MBI->tag_get_length ( tag, tagLength );
	CHKERRVAL(rval);

    DataType  dtype;
    rval = context.MBI->tag_get_data_type ( tag, dtype );

    if ( MB_SUCCESS != rval || dtype != MB_TYPE_DOUBLE )
    { return 1; }

    // set it on a subset of entities, based on type and length
    Range* ents_to_set = NULL;

    if ( * ent_type == 0 ) // vertices
    { ents_to_set = &data.all_verts; }
    else if ( * ent_type == 1 )
    { ents_to_set = &data.primary_elems; }

    int nents_to_be_set = *num_tag_storage_length / tagLength;

    if ( nents_to_be_set > ( int ) ents_to_set->size() || nents_to_be_set < 1 )
    { return 1; } // to many entities to be set

    // restrict the range; everything is contiguous; or not?

    Range contig_range ( * ( ents_to_set->begin() ), * ( ents_to_set->begin() + nents_to_be_set - 1 ) );

    rval = context.MBI->tag_set_data ( tag, contig_range, tag_storage_data );
	CHKERRVAL(rval);

    return 0; // no error
}

ErrCode iMOAB_GetDoubleTagStorage ( iMOAB_AppID pid, const iMOAB_String tag_storage_name, int* num_tag_storage_length, int* ent_type, double* tag_storage_data, int tag_storage_name_length )
{
    // exactly the same code, except tag type check
    std::string tag_name ( tag_storage_name );

    if ( tag_storage_name_length < ( int ) tag_name.length() )
    {
        tag_name = tag_name.substr ( 0, tag_storage_name_length );
    }

    appData& data = context.appDatas[*pid];

    if ( data.tagMap.find ( tag_name ) == data.tagMap.end() )
    { return 1; } // tag not defined

    Tag tag =  data.tagMap[tag_name];

    int tagLength = 0;
    ErrorCode rval = context.MBI->tag_get_length ( tag, tagLength );
	CHKERRVAL(rval);

    DataType  dtype;
    rval = context.MBI->tag_get_data_type ( tag, dtype );

    if ( MB_SUCCESS != rval || dtype != MB_TYPE_DOUBLE )
    { return 1; }

    // set it on a subset of entities, based on type and length
    Range* ents_to_get = NULL;

    if ( * ent_type == 0 ) // vertices
    { ents_to_get = &data.all_verts; }
    else if ( * ent_type == 1 )
    { ents_to_get = &data.primary_elems; }

    int nents_to_get = *num_tag_storage_length / tagLength;

    if ( nents_to_get > ( int ) ents_to_get->size() || nents_to_get < 1 )
    { return 1; } // to many entities to get

    // restrict the range; everything is contiguous; or not?

    Range contig_range ( * ( ents_to_get->begin() ), * ( ents_to_get->begin() + nents_to_get - 1 ) );
    rval = context.MBI->tag_get_data ( tag, contig_range, tag_storage_data );
	CHKERRVAL(rval);

    return 0; // no error
}

ErrCode iMOAB_SynchronizeTags ( iMOAB_AppID pid, int* num_tag, int* tag_indices, int* ent_type )
{
#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    Range ent_exchange;
    std::vector<Tag> tags;

    for ( int i = 0; i < * num_tag; i++ )
    {
        if ( tag_indices[i] < 0 || tag_indices[i] >= ( int ) data.tagList.size() )
        { return 1 ; } // error in tag index

        tags.push_back ( data.tagList[tag_indices[i]] );
    }

    if ( * ent_type == 0 )
    { ent_exchange = data.all_verts; }
    else if ( *ent_type == 1 )
    { ent_exchange = data.primary_elems; }
    else
    { return 1; } // unexpected type

    ParallelComm* pco = context.pcomms[*pid];

    ErrorCode rval = pco->exchange_tags ( tags, tags, ent_exchange );
	CHKERRVAL(rval);

#else
    /* do nothing if serial */
    // just silence the warning
    // do not call sync tags in serial!
    int k = *pid + *num_tag + *tag_indices + *ent_type; k++;
#endif

    return 0;
}

ErrCode iMOAB_ReduceTagsMax ( iMOAB_AppID pid, int* tag_index, int* ent_type )
{

#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    Range ent_exchange;
    std::vector<Tag> tags;


	if ( *tag_index < 0 || *tag_index >= ( int ) data.tagList.size() )
        { return 1 ; } // error in tag index

    tags.push_back ( data.tagList[*tag_index] );


    if ( * ent_type == 0 )
    { ent_exchange = data.all_verts; }
    else if ( *ent_type == 1 )
    { ent_exchange = data.primary_elems; }
    else
    { return 1; } // unexpected type

    ParallelComm* pco = context.pcomms[*pid];
    // we could do different MPI_Op; do not bother now, we will call from fortran, and
    ErrorCode rval = pco->reduce_tags ( tags, tags, MPI_MAX, ent_exchange );
	CHKERRVAL(rval);

#else
    /* do nothing if serial */
    // just silence the warning
    // do not call sync tags in serial!
    int k = *pid + *tag_index + *ent_type; k++; // just do junk, to avoid complaints
#endif
    return 0;
}


ErrCode iMOAB_GetNeighborElements ( iMOAB_AppID pid, iMOAB_LocalID* local_index, int* num_adjacent_elements, iMOAB_LocalID* adjacent_element_IDs )
{
    ErrorCode rval;

    // one neighbor for each subentity of dimension-1
    MeshTopoUtil mtu ( context.MBI );
    appData& data = context.appDatas[*pid];
    EntityHandle eh = data.primary_elems[*local_index];
    Range adjs;
    rval = mtu.get_bridge_adjacencies ( eh, data.dimension - 1, data.dimension, adjs );CHKERRVAL(rval);

    if ( * num_adjacent_elements < ( int ) adjs.size() )
    { return 1; } // not dimensioned correctly

    *num_adjacent_elements = ( int ) adjs.size();

    for ( int i = 0; i < * num_adjacent_elements; i++ )
    {
        adjacent_element_IDs[i] = data.primary_elems.index ( adjs[i] );
    }

    return 0;
}

#if 0

ErrCode iMOAB_GetNeighborVertices ( iMOAB_AppID pid, iMOAB_LocalID* local_vertex_ID, int* num_adjacent_vertices, iMOAB_LocalID* adjacent_vertex_IDs )
{
    return 0;
}

#endif


ErrCode iMOAB_CreateVertices ( iMOAB_AppID pid, int* coords_len, int* dim, double* coordinates )
{
    ErrorCode rval;
    appData& data = context.appDatas[*pid];

    if ( !data.local_verts.empty() ) // we should have no vertices in the app
    { return 1; }

    int nverts = *coords_len / *dim;

    rval = context.MBI->create_vertices ( coordinates, nverts, data.local_verts );CHKERRVAL(rval);

    rval = context.MBI->add_entities ( data.file_set, data.local_verts );CHKERRVAL(rval);

    // also add the vertices to the all_verts range
    data.all_verts.merge ( data.local_verts );
    return 0;
}


ErrCode iMOAB_CreateElements ( iMOAB_AppID pid, int* num_elem, int* type,  int* num_nodes_per_element,  int* connectivity,
                               int* block_ID )
{
    // Create elements
    appData& data = context.appDatas[*pid];

    ReadUtilIface* read_iface;
    ErrorCode rval = context.MBI->query_interface ( read_iface );CHKERRVAL(rval);

    EntityType mbtype = ( EntityType ) ( *type );
    EntityHandle actual_start_handle;
    EntityHandle* array = NULL;
    rval = read_iface->get_element_connect ( *num_elem,
            *num_nodes_per_element,
            mbtype,
            1,
            actual_start_handle,
            array );CHKERRVAL(rval);

    // fill up with actual connectivity from input; assume the vertices are in order, and start vertex is
    // the first in the current data vertex range
    EntityHandle firstVertex = data.local_verts[0];

    for ( int j = 0; j < *num_elem * ( *num_nodes_per_element ); j++ )
    { array[j] = connectivity[j] + firstVertex - 1; } // assumes connectivity uses 1 based array (from fortran, mostly)

    Range new_elems ( actual_start_handle, actual_start_handle + *num_elem - 1 );

    rval = context.MBI->add_entities ( data.file_set, new_elems );CHKERRVAL(rval);

    data.primary_elems.merge ( new_elems );

    // add to adjacency
    rval  = read_iface->update_adjacencies(actual_start_handle,
        *num_elem,
        *num_nodes_per_element,
        array); CHKERRVAL(rval);
    // organize all new elements in block, with the given block ID; if the block set is not existing, create  a
    // new mesh set;
    Range sets;
    int set_no = *block_ID;
    const void* setno_ptr = &set_no;
    rval = context.MBI->get_entities_by_type_and_tag ( data.file_set, MBENTITYSET,
            &context.material_tag, &setno_ptr, 1, sets );
    EntityHandle block_set;

    if ( MB_FAILURE == rval || sets.empty() )
    {
        // create a new set, with this block ID
        rval = context.MBI->create_meshset ( MESHSET_SET, block_set );CHKERRVAL(rval);

        rval = context.MBI->tag_set_data ( context.material_tag, &block_set, 1, &set_no );CHKERRVAL(rval);

        // add the material set to file set
        rval = context.MBI->add_entities ( data.file_set, &block_set, 1 );CHKERRVAL(rval);
    }
    else
    { block_set = sets[0]; } // first set is the one we want

    /// add the new ents to the clock set
    rval = context.MBI->add_entities ( block_set, new_elems );CHKERRVAL(rval);

    return 0;
}


// this makes sense only for parallel runs
ErrCode iMOAB_ResolveSharedEntities (  iMOAB_AppID pid, int* num_verts, int* marker )
{
#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    ParallelComm* pco = context.pcomms[*pid];

    // create an integer tag for resolving ; maybe it can be a long tag in the future
    // (more than 2 B vertices;)
    int dum_id = 0;
    Tag stag;
    ErrorCode rval = context.MBI->tag_get_handle ( "__sharedmarker", 1,  MB_TYPE_INTEGER, stag,
                     MB_TAG_CREAT | MB_TAG_DENSE, &dum_id );CHKERRVAL(rval);

    if ( *num_verts > ( int ) data.local_verts.size() )
    { return 1; } // we are not setting the size

    rval = context.MBI->tag_set_data ( stag, data.local_verts, ( void* ) marker ); // assumes integer tag
    EntityHandle cset = data.file_set;
    rval = pco->resolve_shared_ents ( cset, -1, -1, &stag );CHKERRVAL(rval);

    rval = context.MBI->tag_delete(stag); CHKERRVAL(rval);
    // provide partition tag equal to rank
    Tag part_tag;
    dum_id = -1;
    rval = context.MBI->tag_get_handle ( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                                         part_tag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );CHKERRVAL(rval);

    int rank = pco->rank();
    rval = context.MBI->tag_set_data ( part_tag, &cset, 1, &rank );CHKERRVAL(rval);

#endif
    return 0;
}


// this assumes that this was not called before
ErrCode iMOAB_DetermineGhostEntities (  iMOAB_AppID pid, int* ghost_dim, int* num_ghost_layers, int* bridge_dim )
{
    if ( *num_ghost_layers <= 0 )
    { return 0; } // nothing to do

#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    ParallelComm* pco = context.pcomms[*pid];

    int addl_ents = 0; //maybe we should be passing this too; most of the time we do not need additional ents
    ErrorCode rval = pco->exchange_ghost_cells ( *ghost_dim, *bridge_dim,
                     *num_ghost_layers, addl_ents, true, true, &data.file_set ); // collective call

    if ( rval != MB_SUCCESS )
    { return 1; }

    // now re-establish all mesh info; will reconstruct mesh info, based solely on what is in the file set
    int rc = iMOAB_UpdateMeshInfo ( pid );
    return rc;
#endif
    return 0;
}


ErrCode iMOAB_SetGlobalInfo ( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems )
{
    appData& data = context.appDatas[*pid];
    data.num_global_vertices = *num_global_verts;
    data.num_global_elements = *num_global_elems;
    return 0;
}


ErrCode iMOAB_GetGlobalInfo ( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems )
{
    appData& data = context.appDatas[*pid];

    if ( NULL != num_global_verts ) { *num_global_verts = data.num_global_vertices; }

    if ( NULL != num_global_elems ) { *num_global_elems = data.num_global_elements; }

    return 0;
}


#ifdef MOAB_HAVE_MPI

ErrCode iMOAB_SendMesh ( iMOAB_AppID pid, MPI_Comm* global, MPI_Group* receivingGroup, int* rcompid )
{
    int ierr;
    //appData & data = context.appDatas[*pid];
    ParallelComm* pco = context.pcomms[*pid];

    MPI_Comm sender = pco->comm(); // the sender comm is obtained from parallel comm in moab;
    // no need to pass it along
    // first see what are the processors in each group; get the sender group too, from the sender communicator
    MPI_Group senderGroup;
    ierr = MPI_Comm_group ( sender, &senderGroup );
    if ( ierr != 0 ) return 1;

    // instantiate the par comm graph
    // ParCommGraph::ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2, int coid1, int coid2)
    ParCommGraph* cgraph = new ParCommGraph ( *global, senderGroup, *receivingGroup, context.appDatas[*pid].external_id, *rcompid );
    // we should search if we have another pcomm with the same comp ids in the list already
    // sort of check existing comm graphs in the list context.appDatas[*pid].pgraph
    context.appDatas[*pid].pgraph.push_back ( cgraph );

    int sender_rank = -1;
    MPI_Comm_rank ( sender, &sender_rank );

    // decide how to distribute elements to each processor
    // now, get the entities on local processor, and pack them into a buffer for various processors
    // we will do trivial partition: first get the total number of elements from "sender"
    std::vector<int> number_elems_per_part;
    // how to distribute local elements to receiving tasks?
    // trivial partition: compute first the total number of elements need to be sent
    Range& owned = context.appDatas[*pid].owned_elems;

    int local_owned_elem = ( int ) owned.size();
    int size = pco->size();
    int rank = pco->rank();
    number_elems_per_part.resize ( size ); //
    number_elems_per_part[rank] = local_owned_elem;
#if (MPI_VERSION >= 2)
    // Use "in place" option
    ierr = MPI_Allgather ( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                           &number_elems_per_part[0], 1, MPI_INTEGER,
                           sender );
#else
    {
        std::vector<int> all_tmp ( size );
        ierr = MPI_Allgather ( &number_elems_per_part[rank], 1, MPI_INTEGER,
                               &all_tmp[0], 1, MPI_INTEGER,
                               sender );
        number_elems_per_part = all_tmp;
    }
#endif

    if ( ierr != 0 )
    { return 1; }

    // every sender computes the trivial partition, it is cheap, and we need to send it anyway to each sender
    ErrorCode rval = cgraph->compute_trivial_partition ( number_elems_per_part );

    if ( MB_SUCCESS != rval )
    { return 1; }

    rval = cgraph->send_graph ( *global );

    if ( MB_SUCCESS != rval )
    { return 1; }

    // pco is needed to pack, not for communication
    rval = cgraph->send_mesh_parts ( *global, pco, owned );

    if ( rval != MB_SUCCESS ) { return 1; }

    // mark for deletion
    MPI_Group_free(&senderGroup);
    return 0;
}


ErrCode iMOAB_ReceiveMesh ( iMOAB_AppID pid, MPI_Comm* global, MPI_Group* sendingGroup,
                            int* scompid )
{
    appData& data = context.appDatas[*pid];
    ParallelComm* pco = context.pcomms[*pid];
    MPI_Comm receive = pco->comm();
    EntityHandle local_set = data.file_set;
    ErrorCode rval;

    // first see what are the processors in each group; get the sender group too, from the sender communicator
    MPI_Group receiverGroup;
    int ierr = MPI_Comm_group ( receive, &receiverGroup );

    if ( ierr != 0 )
    { return 1; }

    // instantiate the par comm graph
    ParCommGraph* cgraph = new ParCommGraph ( *global, *sendingGroup, receiverGroup, *scompid, context.appDatas[*pid].external_id );
    // TODO we should search if we have another pcomm with the same comp ids in the list already
    // sort of check existing comm graphs in the list context.appDatas[*pid].pgraph
    context.appDatas[*pid].pgraph.push_back ( cgraph );

    int receiver_rank = -1;
    MPI_Comm_rank ( receive, &receiver_rank );

    // first, receive from sender_rank 0, the communication graph (matrix), so each receiver
    // knows what data to expect
    std::vector<int> pack_array;
    rval = cgraph->receive_comm_graph ( *global, pco, pack_array );
    if ( MB_SUCCESS != rval ) { return 1; }

    // senders across for the current receiver
    int current_receiver = cgraph->receiver ( receiver_rank );

    std::vector<int> senders_local;
    size_t n = 0;

    while ( n < pack_array.size() )
    {
        if ( current_receiver == pack_array[n] )
        {
            for ( int j = 0; j < pack_array[n + 1]; j++ )
            { senders_local.push_back ( pack_array[n + 2 + j] ); }

            break;
        }

        n = n + 2 + pack_array[n + 1];
    }

#ifdef VERBOSE
    std:: cout << " receiver " << current_receiver << " at rank " <<
               receiver_rank << " will receive from " << senders_local.size() << " tasks: ";

    for ( int k = 0; k < ( int ) senders_local.size(); k++ )
    { std::cout << " " << senders_local[k]; }

    std::cout << "\n";
#endif

    rval = cgraph->receive_mesh ( *global, pco, local_set, senders_local );

    if ( MB_SUCCESS != rval ) { return 1; }

    // after we are done, we could merge vertices that come from different senders, but
    // have the same global id
    Tag idtag;
    rval = context.MBI->tag_get_handle ( "GLOBAL_ID", idtag );

    if ( MB_SUCCESS != rval ) { return 1; }

    if ( ( int ) senders_local.size() >= 2 ) // need to remove duplicate vertices
        // that might come from different senders
    {
        Range local_ents;
        rval = context.MBI->get_entities_by_handle ( local_set, local_ents );

        if ( MB_SUCCESS != rval ) { return 1; }

        Range local_verts = local_ents.subset_by_type ( MBVERTEX );
        Range local_elems = subtract ( local_ents, local_verts );

        // remove from local set the vertices
        rval = context.MBI->remove_entities ( local_set, local_verts );

        if ( MB_SUCCESS != rval ) { return 1; }

#ifdef VERBOSE
        std::cout << "current_receiver " << current_receiver << " local verts: " << local_verts.size() << "\n";
#endif
        MergeMesh mm ( context.MBI );

        rval = mm.merge_using_integer_tag ( local_verts, idtag );

        if ( MB_SUCCESS != rval ) { return 1; }

        Range new_verts; // local elems are local entities without vertices
        rval = context.MBI->get_connectivity ( local_elems, new_verts );

        if ( MB_SUCCESS != rval ) { return 1; }

#ifdef VERBOSE
        std::cout << "after merging: new verts: " << new_verts.size() << "\n";
#endif
        rval = context.MBI->add_entities ( local_set, new_verts );

        if ( MB_SUCCESS != rval ) { return 1; }
    }

    // still need to resolve shared entities (in this case, vertices )
    rval = pco->resolve_shared_ents ( local_set, -1, -1, &idtag );

    if ( rval != MB_SUCCESS ) { return 1; }

    // set the parallel partition tag
    Tag part_tag;
	int dum_id = -1;
	rval = context.MBI->tag_get_handle ( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
										 part_tag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );CHKERRVAL(rval);

	int rank = pco->rank();
	rval = context.MBI->tag_set_data ( part_tag, &local_set, 1, &rank );CHKERRVAL(rval);

    // populate the mesh with current data info
    ierr = iMOAB_UpdateMeshInfo ( pid );

    if ( 0 != ierr ) { return 1; }

    // mark for deletion
    MPI_Group_free(&receiverGroup);

    return 0;
}

ErrCode FindParCommGraph(iMOAB_AppID pid, int *scompid, int *rcompid, ParCommGraph *& cgraph, int * sense)
{
  //appData& data = context.appDatas[*pid];
  cgraph = NULL;
  //ParallelComm* pco = context.pcomms[*pid];
  std::vector<ParCommGraph*> & vpg = context.appDatas[*pid].pgraph;
  size_t i = -1;
  *sense = 0;
  for (i=0; i<vpg.size(); i++)
  {
    ParCommGraph * pg=vpg[i];
    if ( (pg-> get_component_id1() == *scompid )&& (pg-> get_component_id2() == *rcompid ))
    {
      cgraph = pg;
      *sense = 1;
      break;
    }
    if ( (pg-> get_component_id2() == *scompid )&& (pg-> get_component_id1() == *rcompid ))
    {
      cgraph = pg;
      *sense = -1;
      break;
    }
  }
  if ( i< vpg.size() && NULL!=cgraph )
    return 0;

  return 1; // error, we did not find cgraph
}

ErrCode iMOAB_SendElementTag(iMOAB_AppID pid, int* scompid, int* rcompid, const iMOAB_String tag_storage_name,
    MPI_Comm* join, int tag_storage_name_length)
{
  // first, based on the scompid and rcompid, find the parCommGraph corresponding to this exchange
  // instantiate the par comm graph
  // ParCommGraph::ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2, int coid1, int coid2)
  ParCommGraph* cgraph = NULL;
  int sense  = 0;
  int ierr = FindParCommGraph(pid, scompid, rcompid, cgraph, &sense);
  if ( 0 != ierr || NULL == cgraph ) { return 1; }

  ParallelComm* pco = context.pcomms[*pid];
  Range& owned = context.appDatas[*pid].owned_elems;

  std::string tag_name ( tag_storage_name );

  if ( tag_storage_name_length < ( int ) strlen ( tag_storage_name ) )
  {
      tag_name = tag_name.substr ( 0, tag_storage_name_length );
  }
  Tag tagHandle;
  // basically, we assume everything is defined already on the tag,
  //   and we can get the tag just by its name
  ErrorCode rval = context.MBI->tag_get_handle ( tag_name.c_str(), tagHandle);
  if ( MB_SUCCESS != rval || NULL == tagHandle) { return 1; }
  // pco is needed to pack, and for moab instance, not for communication!
  // still use nonblocking communication, over the
  rval = cgraph->send_tag_values ( *join, pco, owned, tagHandle );

  if ( MB_SUCCESS != rval ) { return 1; }
  // now, send to each corr_tasks[i] tag data for corr_sizes[i] primary entities

  return 0;
}

ErrCode iMOAB_ReceiveElementTag(iMOAB_AppID pid, int* scompid, int* rcompid, const iMOAB_String tag_storage_name,
    MPI_Comm* join, int tag_storage_name_length)
{
  // first, based on the scompid and rcompid, find the parCommGraph corresponding to this exchange
  // instantiate the par comm graph
  // ParCommGraph::ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2, int coid1, int coid2)
  ParCommGraph* cgraph = NULL;
  int sense  = 0;
  int ierr = FindParCommGraph(pid, scompid, rcompid, cgraph, &sense);
  if ( 0 != ierr || NULL == cgraph ) { return 1; }

  ParallelComm* pco = context.pcomms[*pid];
  Range& owned = context.appDatas[*pid].owned_elems;

  std::string tag_name ( tag_storage_name );

  if ( tag_storage_name_length < ( int ) strlen ( tag_storage_name ) )
  {
      tag_name = tag_name.substr ( 0, tag_storage_name_length );
  }
  Tag tagHandle;
  // basically, we assume everything is defined already on the tag,
  //   and we can get the tag just by its name
  ErrorCode rval = context.MBI->tag_get_handle ( tag_name.c_str(), tagHandle);
  if ( MB_SUCCESS != rval || NULL == tagHandle) { return 1; }
  // pco is needed to pack, and for moab instance, not for communication!
  // still use nonblocking communication, over the
  rval = cgraph->receive_tag_values ( *join, pco, owned, tagHandle );

  if ( MB_SUCCESS != rval ) { return 1; }
  // now, send to each corr_tasks[i] tag data for corr_sizes[i] primary entities

  return 0;
}


ErrCode iMOAB_FreeSenderBuffers ( iMOAB_AppID pid, MPI_Comm* join, int* rcompid )
{
    // need first to find the pgraph that holds the information we need
    // this will be called on sender side only
    appData& data = context.appDatas[*pid];
    std::vector<ParCommGraph*>& pgrs = context.appDatas[*pid].pgraph;
    int ext_id = data.external_id;
    ParCommGraph* pg = NULL;

    for ( size_t i = 0; i < pgrs.size(); i++ )
    {
        if (  ( pgrs[i]->get_component_id2() == *rcompid  && ext_id ==  pgrs[i]->get_component_id1() ) ||
            ( pgrs[i]->get_component_id2() == ext_id  && *rcompid ==  pgrs[i]->get_component_id1())) // sense -1
        {
            pg =  pgrs[i];
            break;
        }
    }

    // if not found, problem
    if ( pg == NULL )
    { return 1; } // cannot find the graph

    pg->release_send_buffers ( *join );
    return 0;
}

#ifdef MOAB_HAVE_TEMPESTREMAP

#define USE_API


static ErrCode ComputeSphereRadius ( iMOAB_AppID pid, double* radius)
{
    ErrorCode rval;
    double coordinates[3];

    Range& verts = context.appDatas[*pid].all_verts;
    moab::EntityHandle firstVertex = (verts[0]);

    // coordinate data
    rval = context.MBI->get_coords ( &(firstVertex), 1, coordinates );CHKERRVAL(rval);

    // compute the distance from origin
    *radius = coordinates[0]*coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];

    // TODO: we could do this in a loop to verify if the pid represents a spherical mesh

    return 0;
}

ErrCode iMOAB_ComputeMeshIntersectionOnSphere ( iMOAB_AppID pid_src, iMOAB_AppID pid_tgt, iMOAB_AppID pid_intx)
{
    ErrorCode rval;

    double radius_source=1.0;
    double radius_target=1.0;
    const double epsrel=1e-8;
    const double boxeps=1.e-5;

    // Get the source and target data and pcomm objects
    appData& data_src = context.appDatas[*pid_src];
    appData& data_tgt = context.appDatas[*pid_tgt];
	appData& data_intx = context.appDatas[*pid_intx];
    ParallelComm* pco_intx = context.pcomms[*pid_intx];

	//  Sanity check: Check that the source and target meshes belong to the same pes. 
    //  assert(pco_src->get_id() == pco_tgt->get_id());
    //  assert(pco_src->get_id() == pco_intx->get_id());

    // Mesh intersection has already been computed; Return early.
    if(data_intx.remapper != NULL) return 0;

    rval = pco_intx->check_all_shared_handles();CHKERRVAL(rval);
    
	// print verbosely about the problem setting
	{
		moab::Range rintxverts, rintxelems;
		rval = context.MBI->get_entities_by_dimension ( data_src.file_set, 0, rintxverts );CHKERRVAL(rval);
		rval = context.MBI->get_entities_by_dimension ( data_src.file_set, 2, rintxelems );CHKERRVAL(rval);
		rval = fix_degenerate_quads ( context.MBI, data_src.file_set );CHKERRVAL(rval);
		rval = positive_orientation ( context.MBI, data_src.file_set, radius_source );CHKERRVAL(rval);
		ErrCode ierr = iMOAB_UpdateMeshInfo(pid_src); CHKIERRVAL(ierr);
#ifdef VERBOSE
 		std::cout << "The red set contains " << rintxverts.size() << " vertices and " << rintxelems.size() << " elements \n";
#endif

		moab::Range bintxverts, bintxelems;
		rval = context.MBI->get_entities_by_dimension ( data_tgt.file_set, 0, bintxverts );CHKERRVAL(rval);
		rval = context.MBI->get_entities_by_dimension ( data_tgt.file_set, 2, bintxelems );CHKERRVAL(rval);
		rval = fix_degenerate_quads ( context.MBI, data_tgt.file_set );CHKERRVAL(rval);
		rval = positive_orientation ( context.MBI, data_tgt.file_set, radius_target );CHKERRVAL(rval);
		ierr = iMOAB_UpdateMeshInfo(pid_tgt); CHKIERRVAL(ierr);
#ifdef VERBOSE
 		std::cout << "The blue set contains " << bintxverts.size() << " vertices and " << bintxelems.size() << " elements \n";
#endif
	}

	// set the context for the source and destination applications
	data_intx.pid_src = pid_src;
	data_intx.pid_dest = pid_tgt;

	// Now allocate and initialize the remapper object
    data_intx.remapper = new moab::TempestRemapper ( context.MBI, pco_intx );
    data_intx.remapper->meshValidate = true;
    data_intx.remapper->constructEdgeMap = true;
    
    // Do not create new filesets; Use the sets from our respective applications
    data_intx.remapper->initialize(false);
    data_intx.remapper->GetMeshSet ( moab::Remapper::SourceMesh ) = data_src.file_set;
    data_intx.remapper->GetMeshSet ( moab::Remapper::TargetMesh ) = data_tgt.file_set;
    data_intx.remapper->GetMeshSet ( moab::Remapper::IntersectedMesh ) = data_intx.file_set;

    // Rescale the radius of both to compute the intersection
    ComputeSphereRadius(pid_src, &radius_source);
    ComputeSphereRadius(pid_tgt, &radius_target);
    std::cout << "Radius of spheres: source = " << radius_source << " and target = " << radius_target << "\n";

    /* Let make sure that the radius match for source and target meshes. If not, rescale now and unscale later. */
    bool radii_scaled = false;
    if (fabs(radius_source - radius_target) > 1e-10) { /* the radii are different */
        radii_scaled = true;
        rval = ScaleToRadius(context.MBI, data_src.file_set, 1.0);CHKERRVAL(rval);
        rval = ScaleToRadius(context.MBI, data_tgt.file_set, 1.0);CHKERRVAL(rval);
    }

    // Rescale the radius of both to compute the intersection
    ComputeSphereRadius(pid_src, &radius_source);
    ComputeSphereRadius(pid_tgt, &radius_target);
    std::cout << "Radius of spheres: source = " << radius_source << " and target = " << radius_target << "\n";
#ifdef USE_API

    rval = data_intx.remapper->ConvertMeshToTempest ( moab::Remapper::SourceMesh );CHKERRVAL(rval);
    rval = data_intx.remapper->ConvertMeshToTempest ( moab::Remapper::TargetMesh );CHKERRVAL(rval);

	// Compute intersections with MOAB
	rval = data_intx.remapper->ComputeOverlapMesh ( epsrel, 1.0, 1.0, boxeps, false );CHKERRVAL(rval);
    // rval = data_intx.remapper->ConvertMeshToTempest ( moab::Remapper::IntersectedMesh );CHKERRVAL(rval);

#else    
    // Create the intersection object on the sphere between two meshes
	// setup the intersector 
	moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere( context.MBI );
	mbintx->set_error_tolerance ( epsrel );
	mbintx->set_box_error ( boxeps );
	mbintx->set_radius_source_mesh ( 1.0 );
    mbintx->set_radius_destination_mesh ( 1.0 );
	mbintx->set_parallel_comm ( pco_intx );

	rval = mbintx->FindMaxEdges ( data_src.file_set, data_tgt.file_set );CHKERRVAL(rval);

	// Migrate the meshes locally so that we have full coverage of the source meshset
	moab::Range local_verts;
	rval = mbintx->build_processor_euler_boxes ( data_tgt.file_set, local_verts );CHKERRVAL(rval);

	// Compute the covering set
	moab::EntityHandle covering_set;
    rval = context.MBI->create_meshset ( moab::MESHSET_SET, covering_set );CHKERRVAL(rval);

	// This step involves lots of communication if mesh is distributed very differently
	rval = mbintx->construct_covering_set ( data_src.file_set, covering_set );CHKERRVAL(rval);

	// Now let's invoke the MOAB intersection algorithm in parallel with a
	// source and target mesh set representing two different decompositions
	rval = mbintx->intersect_meshes ( covering_set, data_tgt.file_set, data_intx.file_set );CHKERRVAL(rval);

	// free the memory
	delete mbintx;
#endif

    {
        // Now let us re-convert the MOAB mesh back to Tempest representation
        rval = data_intx.remapper->AssociateSrcTargetInOverlap();CHKERRVAL(rval);
        rval = data_intx.remapper->ConvertMOABMesh_WithSortedEntitiesBySource();CHKERRVAL(rval);
    }

    // Set the context for the OfflineMap computation
    data_intx.weightMap = new moab::TempestOfflineMap ( data_intx.remapper );

    if (radii_scaled) { /* the radii are different, so lets rescale back */
        rval = ScaleToRadius(context.MBI, data_src.file_set, radius_source);CHKERRVAL(rval);
        rval = ScaleToRadius(context.MBI, data_tgt.file_set, radius_target);CHKERRVAL(rval);
    }

    return 0;
}

// this call must be collective on the joint communicator
//  intersection tasks on coupler will need to send to the components tasks the list of
// id elements that are relevant: they intersected some of the target elements (which are not needed here)
//  in the intersection
ErrCode iMOAB_CoverageGraph(MPI_Comm* join, iMOAB_AppID pid_src, int* scompid, iMOAB_AppID pid_migr,
    int* migrcomp, iMOAB_AppID pid_intx)
{
  // first, based on the scompid and migrcomp, find the parCommGraph corresponding to this exchange

  ErrorCode rval;
  std::vector<int> srcSenders;
  std::vector<int> receivers;
  ParCommGraph* sendGraph = NULL;
  int sense = 0;
  int ierr = FindParCommGraph(pid_src, scompid, migrcomp, sendGraph, &sense);
  if ( 0 != ierr || NULL == sendGraph || sense != 1) {
    std::cout<<" probably not on component source PEs \n";
  }
  else
  {
    // report the sender and receiver tasks in the joint comm
    srcSenders = sendGraph->senders();
    receivers = sendGraph->receivers();
    std::cout << "senders: " << srcSenders.size() << " first sender: "<< srcSenders[0] << std::endl;
  }
  ParCommGraph * recvGraph = NULL;
  int senseRec = 0;
  ierr = FindParCommGraph(pid_migr, scompid, migrcomp, recvGraph, &senseRec);
  if ( 0 != ierr || NULL == recvGraph ) {
    std::cout << " not on receive PEs for migrated mesh \n";
  }
  else
  {
    // report the sender and receiver tasks in the joint comm, from migrated mesh pt of view
    srcSenders = recvGraph->senders();
    receivers = recvGraph->receivers();
    std::cout << "receivers: " << receivers.size() << " first receiver: "<< receivers[0] << std::endl;
  }

  // loop over pid_intx elements, to see what original processors in joint comm have sent the coverage mesh
  // if we are on intx tasks, send coverage info towards original component tasks, about needed cells
  //
  TupleList TLcovIDs;
  TLcovIDs.initialize(2, 0, 0, 0, 0);// to proc, GLOBAL ID; estimate about 100 IDs to be sent
  // will push_back a new tuple, if needed
  TLcovIDs.enableWriteAccess();
  // the crystal router will send ID cell to the original source, on the component task
  // if we are on intx tasks, loop over all intx elements and

  int currentRankInJointComm = -1;
  ierr = MPI_Comm_rank(*join, &currentRankInJointComm);
  if (MPI_SUCCESS != ierr)
    return 1; // fatal error, abort

  // if currentRankInJointComm is in receivers list, it means that we are on intx tasks too, we need to
  // send information towards component tasks
  if ( find(receivers.begin(), receivers.end(), currentRankInJointComm) != receivers.end() ) // we are on receivers tasks, we can request intx info
  {
    // find the pcomm for the intx pid
    if (*pid_intx >= context.appDatas.size() )
      return 1;
    appData& dataIntx = context.appDatas[*pid_intx];
    Tag parentTag ;
    rval = context.MBI->tag_get_handle("BlueParent", parentTag); // id of the blue, source element
    if (MB_SUCCESS != rval || !parentTag)
      return 1; // fatal error, abort
    Tag orgSendProcTag ;
    rval = context.MBI->tag_get_handle("orig_sending_processor", orgSendProcTag);
    if ( MB_SUCCESS != rval || !orgSendProcTag)
      return 1; // fatal error, abort
    // find the file set, red parents for intx cells, and put them in tuples
    EntityHandle intxSet=dataIntx.file_set;
    // get all entities from the set, and look at their RedParent
    Range cells;
    rval = context.MBI->get_entities_by_dimension(intxSet, 2, cells);
    if (MB_SUCCESS != rval)
      return 1; // fatal error, abort
    std::map<int, std::set<int> > idsFromProcs; // send that info back to enhance parCommGraph cache
    for (Range::iterator it=cells.begin(); it!=cells.end(); it++)
    {
      EntityHandle intx_cell = *it;
      int gidCell, origProc; // look at receivers
      rval = context.MBI->tag_get_data(parentTag, &intx_cell, 1, &gidCell);
      if (MB_SUCCESS != rval)
        return 1;
      rval = context.MBI->tag_get_data(orgSendProcTag, &intx_cell, 1, &origProc); // in the
      if (MB_SUCCESS != rval)
        return 1;
      std::set<int> &setInts = idsFromProcs[origProc];
      setInts.insert(gidCell);
      //std::cout << origProc << " id:" << gidCell << " size: " << setInts.size() << std::endl;
    }

    std::cout<<" map size:" << idsFromProcs.size() << std::endl;
    // arrange in tuples , use map iterators to send the ids
    for (std::map<int, std::set<int> >::iterator mit = idsFromProcs.begin(); mit!=idsFromProcs.end(); mit++)
    {
      int procToSendTo = mit->first;
      std::set<int> & idSet = mit->second;
      for (std::set<int>::iterator sit=idSet.begin(); sit!=idSet.end(); sit++)
      {
        int n=TLcovIDs.get_n();
        TLcovIDs.reserve();
        TLcovIDs.vi_wr[2*n] = procToSendTo; // send to processor
        TLcovIDs.vi_wr[2*n+1] = *sit; // global id needs index in the local_verts range
      }
    }
  }
  ProcConfig pc(*join); // proc config does the crystal router
  pc.crystal_router()->gs_transfer(1, TLcovIDs, 0); // communication towards component tasks, with what ids are needed
  // for each task from receiver

  // a test to know if we are on the sender tasks (original component, in this case, atmosphere)
  if (NULL != sendGraph)
  {
    // collect TLcovIDs tuple, will set in a local map/set, the ids that are sent to each receiver task
    rval = sendGraph->settle_send_graph(TLcovIDs);
    if (MB_SUCCESS != rval)
      return 1;
  }
  return 0;// success

}

ErrCode iMOAB_ComputeScalarProjectionWeights ( iMOAB_AppID pid_intx, 
                                               const iMOAB_String disc_method_source, int* disc_order_source,
                                               const iMOAB_String disc_method_target, int* disc_order_target,
                                               int* fVolumetric, int* fNoConservation,
                                               int* fValidate,
                                               const iMOAB_String source_soln_tag_dof_name,
                                               const iMOAB_String target_soln_tag_dof_name,
                                               int disc_method_source_length,
                                               int disc_method_target_length,
                                               int source_soln_tag_dof_name_length,
                                               int target_soln_tag_dof_name_length)
{
	moab::ErrorCode rval;
	
	assert(disc_order_source && disc_order_target && *disc_order_source > 0 && *disc_order_target > 0);
	assert(disc_method_source_length > 0 && disc_method_target_length > 0);
    assert(source_soln_tag_dof_name_length > 0 && target_soln_tag_dof_name_length > 0);
    
    // Get the source and target data and pcomm objects
	appData& data_intx = context.appDatas[*pid_intx];
    ParallelComm* pco_intx = context.pcomms[*pid_intx];

	// Now allocate and initialize the remapper object
    moab::TempestRemapper* remapper = data_intx.remapper;

	// Setup computation of weights
	// Call to generate an offline map with the tempest meshes
    moab::TempestOfflineMap* weightMap = data_intx.weightMap;
    assert(weightMap != NULL);

    // Now let us compute the local-global mapping and store it in the context
    // We need this mapping when computing matvec products and to do reductions in parallel
	// Additionally, the call below will also compute weights with TempestRemap
	rval = weightMap->GenerateOfflineMap ( std::string(disc_method_source), std::string(disc_method_target),        // std::string strInputType, std::string strOutputType,
										   (*disc_order_source), (*disc_order_target),    // const int nPin, const int nPout,
                                           false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
										   (fVolumetric ? *fVolumetric > 0 : false),  // bool fVolumetric=false, 
                                           (fNoConservation ? *fNoConservation > 0 : false), // bool fNoConservation=false, 
                                           false, // bool fNoCheck=false,
                                           source_soln_tag_dof_name, target_soln_tag_dof_name,
										   "", //"",   // std::string strVariables="", std::string strOutputMap="",
										   "", "",   // std::string strInputData="", std::string strOutputData="",
										   "", false,  // std::string strNColName="", bool fOutputDouble=false,
										   "", false, 0.0,   // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0,
										   false, false   // bool fInputConcave = false, bool fOutputConcave = false
										 );CHKERRVAL(rval);

    // Mapping computation done
	if (fValidate && *fValidate)
	{
		const double radius = 1.0 /*2.0*acos(-1.0)*/;
		double local_areas[3], global_areas[3]; // Array for Initial area, and through Method 1 and Method 2
		local_areas[0] = area_on_sphere_lHuiller ( context.MBI, context.appDatas[*(context.appDatas[*pid_intx].pid_src)].file_set, radius );
		local_areas[1] = area_on_sphere_lHuiller ( context.MBI, context.appDatas[*pid_intx].file_set, radius );
		local_areas[2] = area_on_sphere ( context.MBI, context.appDatas[*pid_intx].file_set, radius );

		MPI_Allreduce ( &local_areas, &global_areas, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

		if ( !pco_intx->rank() )
		{
			printf ( "initial area: %12.10f\n", global_areas[0] );
			printf ( "area with l'Huiller: %12.10f with Girard: %12.10f\n", global_areas[1], global_areas[2] );
			printf ( "  relative difference areas = %12.10e\n", fabs ( global_areas[1] - global_areas[2] ) / global_areas[1] );
			printf ( "  relative error = %12.10e\n", fabs ( global_areas[1] - global_areas[0] ) / global_areas[1] );
		}
	}
	
	return 0;
}


ErrCode iMOAB_ApplyScalarProjectionWeights (   iMOAB_AppID pid_intersection, 
                                               const iMOAB_String solution_tag_name,
                                               int solution_tag_name_length )
{
    moab::ErrorCode rval;

    assert(solution_tag_name_length > 0);

    // Get the source and target data and pcomm objects
    appData& data_intx = context.appDatas[*pid_intersection];
    ParallelComm* pco_intx = context.pcomms[*pid_intersection];

    // Now allocate and initialize the remapper object
    moab::TempestRemapper* remapper = data_intx.remapper;
    moab::TempestOfflineMap* weightMap = data_intx.weightMap;

    /* Global ID - exchange data for covering data */
    Tag solnTag;
    rval = context.MBI->tag_get_handle ( solution_tag_name, solnTag );CHKERRVAL(rval);

    // moab::HypreParVector sVals(weightMap->GetRowVector()), tVals(weightMap->GetColVector());
    std::vector<double> solSTagVals(weightMap->GetSourceLocalNDofs());
    std::vector<double> solTTagVals(weightMap->GetDestinationLocalNDofs());

    moab::Range& covSrcEnts = remapper->GetMeshEntities(moab::Remapper::CoveringMesh);
    moab::Range& tgtEnts = remapper->GetMeshEntities(moab::Remapper::TargetMesh);
    // assert(covSrcEnts.size() == )

    rval = context.MBI->tag_get_data ( solnTag, covSrcEnts, &solSTagVals[0] );CHKERRVAL(rval);

    rval = weightMap->ApplyWeights(solSTagVals, solTTagVals);CHKERRVAL(rval);

    rval = context.MBI->tag_set_data ( solnTag, tgtEnts, &solTTagVals[0] );CHKERRVAL(rval);

    return 0;
}


#endif

#endif

#ifdef __cplusplus
}
#endif
