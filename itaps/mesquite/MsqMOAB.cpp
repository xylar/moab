/** \file MsqMOAB.cpp
 *  \brief
 *  \author Vijay Mahadevan
 */
#include "MsqMOAB.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include <assert.h>
#include "MsqIBase.hpp"
#include <algorithm>

#ifdef IMESH_MAJOR_VERSION
# define IMESH_VERSION_ATLEAST(MAJOR,MINOR) \
           1000*IMESH_MAJOR_VERSION+IMESH_MINOR_VERSION <= \
           1000*MAJOR+MINOR
#else
# define IMESH_VERSION_ATLEAST(MAJOR,MINOR) 0
#endif

namespace MBMesquite
{


    /*************************************************************************
     *                          Mesh Definition
     ************************************************************************/

    MsqMOAB::MsqMOAB ( moab::Core* mesh,
                       moab::EntityHandle meshset,
                       moab::EntityType type,
                       MsqError& err,
                       const moab::Tag* fixed_tag,
                       const moab::Tag* slaved_tag )
        : meshInstance ( mesh ),
          inputSetType ( moab::MBMAXTYPE ),
          inputSet ( 0 ),
          byteTag ( 0 ),
          createdByteTag ( false ),
          geometricDimension ( 0 )
    {
        init_active_mesh ( mesh, err, fixed_tag, slaved_tag );
        MSQ_ERRRTN ( err );
        set_active_set ( meshset, type, err );
        MSQ_ERRRTN ( err );
    }

    MsqMOAB::MsqMOAB ( moab::Core* mesh,
                       moab::EntityType type,
                       MsqError& err,
                       const moab::Tag* fixed_tag,
                       const moab::Tag* slaved_tag )
        : meshInstance ( mesh ),
          inputSetType ( moab::MBMAXTYPE ),
          inputSet ( 0 ),
          byteTag ( 0 ),
          createdByteTag ( false ),
          geometricDimension ( 0 )
    {
        init_active_mesh ( mesh, err, fixed_tag, slaved_tag );
        MSQ_ERRRTN ( err );

        moab::EntityHandle root_set = 0 /*meshInstance->get_root_set()*/;
        set_active_set ( root_set, type, err ); MSQ_ERRRTN ( err );
    }

    MsqMOAB::~MsqMOAB()
    {
        moab::ErrorCode ierr;

        if ( createdByteTag )
        {
            ierr = meshInstance->tag_delete ( byteTag ); MB_CHK_ERR_RET ( ierr );
        }
    }

    moab::Core* MsqMOAB::get_interface() const
    {
        return meshInstance;
    }

    moab::EntityHandle MsqMOAB::get_entity_set() const
    {
        return inputSet;
    }

    moab::DataType MsqMOAB::check_valid_flag_tag ( moab::Tag tag,
            const char* /*which_flag*/,
            MsqError& err )
    {
        moab::ErrorCode ierr;
        int size;
        std::string name;
        moab::DataType type = moab::MB_MAX_DATA_TYPE;

        ierr = meshInstance->tag_get_data_type ( tag, type ); MB_CHK_ERR_RET_VAL ( ierr, type );
        ierr = meshInstance->tag_get_name ( tag, name ); MB_CHK_ERR_RET_VAL ( ierr, type );
        ierr = meshInstance->tag_get_length ( tag, size ); MB_CHK_ERR_RET_VAL ( ierr, type );

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
        return type;
    }

    void MsqMOAB::init_active_mesh ( moab::Core* /*mesh*/,
                                     MsqError& err,
                                     const moab::Tag* fixed_tag,
                                     const moab::Tag* slaved_tag )
    {
        moab::ErrorCode ierr;

        // Initialize topology map
        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );

        const size_t mapsize = sizeof ( topologyMap ) / sizeof ( MBMesquite::EntityTopology );

        if ( mapsize < moab::MBMAXTYPE )
        {
            MSQ_SETERR ( err ) ( "MsqMOAB needs to be updated for new iMesh element topologies.",
                                 MsqError::INTERNAL_ERROR );
        }

        for ( size_t i = 0; i <= moab::MBMAXTYPE; ++i )
        { topologyMap[i] = MBMesquite::MIXED; }

        topologyMap[moab::MBTRI     ] = MBMesquite::TRIANGLE;
        topologyMap[moab::MBQUAD    ] = MBMesquite::QUADRILATERAL;
        topologyMap[moab::MBTET     ] = MBMesquite::TETRAHEDRON;
        topologyMap[moab::MBHEX     ] = MBMesquite::HEXAHEDRON;
        topologyMap[moab::MBPRISM   ] = MBMesquite::PRISM;
        topologyMap[moab::MBPYRAMID ] = MBMesquite::PYRAMID;

        // Check that fixed tag is valid
        haveFixedTag = false;

        if ( fixed_tag )
        {
            fixedTagType = check_valid_flag_tag ( *fixed_tag, "fixed", err );
            MSQ_ERRRTN ( err );
            haveFixedTag = true;
            fixedTag = *fixed_tag;
        }

        // Check that slaved tag is valid
        haveSlavedTag = false;

        if ( slaved_tag )
        {
            slavedTagType = check_valid_flag_tag ( *slaved_tag, "slaved", err );
            MSQ_ERRRTN ( err );
            haveSlavedTag = true;
            slavedTag = *slaved_tag;
        }

        // Get/create tag for vertex byte
        ierr = meshInstance->tag_get_handle ( VERTEX_BYTE_TAG_NAME, byteTag );

        if ( moab::MB_SUCCESS != ierr )
        {
            int defval = 0;
            ierr = meshInstance->tag_get_handle ( VERTEX_BYTE_TAG_NAME, 1, moab::MB_TYPE_INTEGER, byteTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &defval, &createdByteTag ); //MB_CHK_ERR_RET(ierr);
        }
        else
        {
            int size;
            moab::DataType type;
            ierr = meshInstance->tag_get_length ( byteTag, size ); MB_CHK_ERR_RET ( ierr );
            ierr = meshInstance->tag_get_data_type ( byteTag, type ); MB_CHK_ERR_RET ( ierr );
        }

        ierr = meshInstance->get_dimension ( geometricDimension ); MB_CHK_ERR_RET ( ierr );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }


    void MsqMOAB::set_fixed_tag ( moab::Tag tag, MsqError& err )
    {
        moab::DataType t = check_valid_flag_tag ( tag, "fixed", err );
        MSQ_ERRRTN ( err );
        fixedTag = tag;
        fixedTagType = t;
        haveFixedTag = true;
    }

    void MsqMOAB::clear_fixed_tag()
    {
        haveFixedTag = false;
    }

    const moab::Tag* MsqMOAB::get_fixed_tag() const
    {
        return haveFixedTag ? &fixedTag : 0;
    }

    void MsqMOAB::set_slaved_tag ( moab::Tag tag, MsqError& err )
    {
        moab::DataType t = check_valid_flag_tag ( tag, "slaved", err );
        MSQ_ERRRTN ( err );
        slavedTag = tag;
        slavedTagType = t;
        haveSlavedTag = true;
    }

    void MsqMOAB::clear_slaved_tag()
    {
        haveSlavedTag = false;
    }

    const moab::Tag* MsqMOAB::get_slaved_tag() const
    {
        return haveSlavedTag ? &slavedTag : 0;
    }



    void MsqMOAB::set_active_set ( moab::EntityHandle elem_set,
                                   moab::EntityType type_in,
                                   MsqError& err )
    {
        inputSetType = type_in;
        inputSet = elem_set;

        // clear vertex byte
        std::vector<VertexHandle> verts;
        get_all_vertices ( verts, err ); MSQ_ERRRTN ( err );

        if ( !verts.empty() )
        {
            std::vector<unsigned char> zeros ( verts.size(), 0 );
            vertices_set_byte ( arrptr ( verts ), arrptr ( zeros ), verts.size(), err );
            MSQ_CHKERR ( err );
        }
    }



    // Returns whether this mesh lies in a 2D or 3D coordinate system.
    int MsqMOAB::get_geometric_dimension ( MBMesquite::MsqError& /*err*/ )
    {
        return geometricDimension;
    }


    //************ Vertex Properties ********************

    void MsqMOAB::get_flag_data ( moab::Tag tag,
                                  bool have_tag,
                                  moab::DataType type,
                                  const VertexHandle vert_array[],
                                  std::vector<bool>& flag_array,
                                  size_t num_vtx,
                                  MsqError& err )
    {
        if ( !num_vtx )
        { return; }

        if ( !have_tag )
        {
            flag_array.clear();
            flag_array.resize ( num_vtx, false );
            return;
        }

        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        flag_array.resize ( num_vtx );

        assert ( sizeof ( VertexHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( vert_array );

        moab::ErrorCode ierr;
        int alloc = num_vtx;
        assert ( ( size_t ) alloc == num_vtx ); // size_t can hold larger values than int if 64-bit

        if ( type == moab::MB_TYPE_INTEGER )
        {
            std::vector<int> values ( num_vtx );
            int* ptr = arrptr ( values );
            ierr = meshInstance->tag_get_data ( tag, arr, num_vtx, ptr ); MB_CHK_ERR_RET ( ierr );

            for ( size_t i = 0; i < num_vtx; ++i )
            { flag_array[i] = !!values[i]; }
        }
        else if ( type == moab::MB_TYPE_OPAQUE )
        {
            std::vector<char> values ( num_vtx );
            void* ptr = arrptr ( values );
            ierr = meshInstance->tag_get_data ( tag, arr, num_vtx, ptr ); MB_CHK_ERR_RET ( ierr );

            for ( size_t i = 0; i < num_vtx; ++i )
            { flag_array[i] = !!values[i]; }
        }
        else
        {
            MSQ_SETERR ( err ) ( "Invalid tag type for vertex flag data", MsqError::INVALID_STATE );
            return ;
        }

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    // Returns true or false, indicating whether the vertex
    // is allowed to be repositioned.  True indicates that the vertex
    // is fixed and cannot be moved.  Note that this is a read-only
    // property; this flag can't be modified by users of the
    // MBMesquite::Mesh interface.
    void MsqMOAB::vertices_get_fixed_flag (
        const VertexHandle vert_array[],
        std::vector<bool>& bool_array,
        size_t num_vtx, MsqError& err )
    {
        get_flag_data ( fixedTag, haveFixedTag, fixedTagType, vert_array, bool_array, num_vtx, err );
    }


    void MsqMOAB::vertices_get_slaved_flag (
        const VertexHandle vert_array[],
        std::vector<bool>& bool_array,
        size_t num_vtx, MsqError& err )
    {
        get_flag_data ( slavedTag, haveSlavedTag, slavedTagType, vert_array, bool_array, num_vtx, err );
    }

    // Get vertex coordinates
    void MsqMOAB::vertices_get_coordinates (
        const MBMesquite::Mesh::VertexHandle vert_array[],
        MsqVertex* coordinates,
        size_t num_vtx,
        MsqError& err )
    {
        if ( !num_vtx )
        { return; }

        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        std::vector<double> dbl_store ( 3 * num_vtx );
        double* dbl_array = arrptr ( dbl_store );

        moab::ErrorCode ierr;
        // int junk = 3*num_vtx;
        assert ( sizeof ( VertexHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( vert_array );
        ierr = meshInstance->get_coords ( arr, num_vtx, dbl_array ); MB_CHK_ERR_RET ( ierr );

        if ( geometricDimension == 2 )
        {
            double* iter = dbl_array;

            for ( size_t i = 0; i < num_vtx; ++i )
            {
                coordinates[i].x ( *iter ); ++iter;
                coordinates[i].y ( *iter ); ++iter;
                coordinates[i].z ( 0 );
            }
        }
        else
        {
            double* iter = dbl_array;

            for ( size_t i = 0; i < num_vtx; ++i )
            {
                coordinates[i].set ( iter );
                iter += 3;
            }
        }

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    void MsqMOAB::vertex_set_coordinates (
        MBMesquite::Mesh::VertexHandle vertex,
        const Vector3D& coords, MsqError& err )
    {
        moab::ErrorCode ierr;
        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        moab::EntityHandle* bh = reinterpret_cast<moab::EntityHandle*> ( &vertex );
        double xyz[3] = {coords[0], coords[1], coords[2]};
        ierr = meshInstance->set_coords ( bh, 1, xyz ); MB_CHK_ERR_RET ( ierr );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    // Each vertex has a byte-sized flag that can be used to store
    // flags.  This byte's value is neither set nor used by the mesh
    // implementation.  It is intended to be used by Mesquite algorithms.
    // Until a vertex's byte has been explicitly set, its value is 0.
    void MsqMOAB::vertex_set_byte (
        MBMesquite::Mesh::VertexHandle vertex,
        unsigned char byte, MsqError& err )
    {
        moab::ErrorCode ierr;
        int value = byte;
        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        moab::EntityHandle* bh = reinterpret_cast<moab::EntityHandle*> ( &vertex );
        ierr = meshInstance->tag_set_data ( byteTag, bh, 1, &value ); MB_CHK_ERR_RET ( ierr );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    void MsqMOAB::vertices_set_byte (
        const VertexHandle* vert_array,
        const unsigned char* byte_array,
        size_t array_size, MsqError& err )
    {
        if ( !array_size )
        { return; }

        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        std::vector<int> data ( array_size );
        std::copy ( byte_array, byte_array + array_size, data.begin() );
        moab::ErrorCode ierr;
        assert ( sizeof ( VertexHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( vert_array );
        ierr = meshInstance->tag_set_data ( byteTag, arr, array_size, arrptr ( data ) ); MB_CHK_ERR_RET ( ierr );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    // Retrieve the byte value for the specified vertex or vertices.
    // The byte value is 0 if it has not yet been set via one of the
    // *_set_byte() functions.
    void MsqMOAB::vertex_get_byte (
        MBMesquite::Mesh::VertexHandle vertex,
        unsigned char* byte, MsqError& err )
    {
        moab::ErrorCode ierr;
        int value;
        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        moab::EntityHandle* bh = reinterpret_cast<moab::EntityHandle*> ( &vertex );
        ierr = meshInstance->tag_get_data ( byteTag, bh, 1, &value ); MB_CHK_ERR_RET ( ierr );
        *byte = value;
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    void MsqMOAB::vertices_get_byte (
        const VertexHandle* vert_array,
        unsigned char* byte_array,
        size_t array_size, MsqError& err )
    {
        if ( !array_size )
        { return; }

        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );
        std::vector<int> data ( array_size );
        moab::ErrorCode ierr;
        int* ptr = arrptr ( data );
        assert ( sizeof ( VertexHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( vert_array );
        ierr = meshInstance->tag_get_data ( byteTag, arr, array_size, ptr ); MB_CHK_ERR_RET ( ierr );
        std::copy ( data.begin(), data.end(), byte_array );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }


    //**************** Topology *****************

    void MsqMOAB::get_adjacent_entities ( const moab::EntityHandle* source,
                                          size_t num_source,
                                          moab::EntityType target_type,
                                          std::vector<EntityHandle>& target,
                                          std::vector<size_t>& offsets,
                                          MsqError& err )
    {
        if ( num_source == 0 )
        {
            target.clear();
            offsets.clear();
            offsets.reserve ( 1 );
            offsets.push_back ( 0 );
            return;
        }

        err.set_error ( MBMesquite::MsqError::UNKNOWN_ERROR );

        moab::ErrorCode ierr;
        int num_adj = 0, num_offset = 0;
        unsigned iadjoff = 0;

        assert ( sizeof ( size_t ) >= sizeof ( int ) );
        offsets.resize ( num_source + 1 );
        int* ptr2;
        bool expand = false;

        if ( sizeof ( size_t ) > sizeof ( int ) )
        {
            ptr2 = ( int* ) malloc ( sizeof ( int ) * ( num_source + 1 ) );
            expand = true;
            num_offset = num_source + 1;
        }
        else
        {
            // sizeof(int) == sizeof(size_t)
            ptr2 = ( int* ) arrptr ( offsets );
            num_offset = offsets.size();
        }

        // std::cout << "Target capacity: " << target.capacity() << " and num sources = " <<  num_source << std::endl;

        assert ( sizeof ( moab::EntityHandle ) == sizeof ( EntityHandle ) );
        bool have_adj = false;

        // If passed vector has allocated storage, try to use existing space
        if ( target.capacity() >= num_source )
        {
            target.resize ( target.capacity() );
            // target.clear();
            moab::EntityHandle* ptr = reinterpret_cast<moab::EntityHandle*> ( arrptr ( target ) );

            ptr2[0] = 0;

            for ( unsigned is = 0; is < num_source; ++is )
            {
                moab::Range adjents;

                if ( target_type == moab::MBVERTEX )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 0, false, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type == moab::MBEDGE )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 1, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type <= moab::MBPOLYGON )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 2, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type < moab::MBENTITYSET )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 3, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else   /* Either EntitySet or MBMaxType -- Failures */
                {
                    MSQ_SETERR ( err ) ( process_itaps_error ( 1 ), MsqError::NOT_IMPLEMENTED ); // "Invalid Target entity type specified"
                    return;
                }

                ptr2[is + 1] = iadjoff + adjents.size();

                for ( unsigned iadj = 0; iadj < adjents.size(); ++iadj, ++iadjoff )
                { ptr[iadjoff] = adjents[iadj]; }

                // target.push_back( static_cast<EntityHandle>(&adjents[iadj]) );
                // std::cout << "1. Source entity [ " << is << "]: n(adjacencies) = " << offsets[is+1] << std::endl;
            }

            if ( moab::MB_SUCCESS == ierr )
            {
                have_adj = true;
                num_adj = iadjoff;
            }
        }

        // If implementation passed back a size, try that
        if ( !have_adj && num_adj && ( unsigned ) num_adj > target.capacity() )
        {
            target.resize ( num_adj );
            // int junk1 = target.capacity(), junk3 = offsets.size();
            moab::EntityHandle* ptr = reinterpret_cast<moab::EntityHandle*> ( arrptr ( target ) );

            ptr2[0] = 0;

            for ( unsigned is = 0; is < num_source; ++is )
            {
                moab::Range adjents;

                if ( target_type == moab::MBVERTEX )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 0, false, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type == moab::MBEDGE )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 1, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type <= moab::MBPOLYGON )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 2, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type < moab::MBENTITYSET )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 3, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else   /* Either EntitySet or MBMaxType -- Failures */
                {
                    MSQ_SETERR ( err ) ( process_itaps_error ( 1 ), MsqError::NOT_IMPLEMENTED ); // "Invalid Target entity type specified"
                    return;
                }

                ptr2[is + 1] = iadjoff + adjents.size();

                for ( unsigned iadj = 0; iadj < adjents.size(); ++iadj, ++iadjoff )
                { ptr[iadjoff] = adjents[iadj]; }
            }

            if ( moab::MB_SUCCESS == ierr )
            {
                have_adj = true;
            }
        }

        // Try with empty sidl array, and copy into elements vector
        if ( !have_adj )
        {
            // If implementation passed back a size, try that
            std::vector<moab::EntityHandle> mArray;

            ptr2[0] = iadjoff;

            for ( unsigned is = 0; is < num_source; ++is )
            {
                moab::Range adjents;

                if ( target_type == moab::MBVERTEX )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 0, false, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type == moab::MBEDGE )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 1, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type <= moab::MBPOLYGON )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 2, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else if ( target_type < moab::MBENTITYSET )
                {
                    ierr = meshInstance->get_adjacencies ( &source[is], 1, 3, true, adjents, moab::Interface::INTERSECT ); //MB_CHK_ERR_RET(ierr);
                }
                else   /* Either EntitySet or MBMaxType -- Failures */
                {
                    MSQ_SETERR ( err ) ( process_itaps_error ( 1 ), MsqError::NOT_IMPLEMENTED ); // "Invalid Target entity type specified"
                    return;
                }

                ptr2[is + 1] = iadjoff + adjents.size();

                for ( unsigned iadj = 0; iadj < adjents.size(); ++iadj, ++iadjoff )
                { mArray.push_back ( adjents[iadj] ); }
            }

            num_adj = mArray.size();
            target.resize ( num_adj );
            std::copy ( mArray.begin(), mArray.end(), reinterpret_cast<moab::EntityHandle*> ( arrptr ( target ) ) );
            mArray.clear();
        }

        if ( expand )
        {
            for ( size_t i = num_offset; i > 0; --i )
            { offsets[i - 1] = static_cast<size_t> ( ptr2[i - 1] ); }

            free ( ptr2 );
        }

        // assert( (unsigned)num_offset == offsets.size() );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }


    void MsqMOAB::vertices_get_attached_elements (
        const VertexHandle* vertices,
        size_t num_vertex,
        std::vector<ElementHandle>& elements,
        std::vector<size_t>& offsets,
        MsqError& err )
    {
        // moab::ErrorCode ierr;
        bool cont;
        assert ( sizeof ( EntityHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* verts = reinterpret_cast<const moab::EntityHandle*> ( vertices );
        get_adjacent_entities ( verts, num_vertex, inputSetType, elements, offsets, err );
        MSQ_ERRRTN ( err );

        moab::EntityHandle root_set = 0 /*meshInstance->get_root_set()*/;

        // Remove all elements not in inputSet
        if ( root_set != inputSet )
        {
            std::vector<size_t>::iterator offset_iter = offsets.begin();
            size_t read_idx, write_idx;
            moab::EntityHandle* bh = reinterpret_cast<moab::EntityHandle*> ( arrptr ( elements ) );

            for ( read_idx = write_idx = 0; read_idx < elements.size(); ++read_idx )
            {
                if ( *offset_iter == read_idx )
                {
                    *offset_iter = write_idx;
                    ++offset_iter;
                }

                cont = meshInstance->contains_entities ( inputSet, &bh[read_idx], 1 );

                if ( cont )
                { elements[write_idx++] = elements[read_idx]; }
            }

            *offset_iter = write_idx;
            elements.resize ( write_idx );
        }
    }


    //**************** Element Topology *****************


    /** Get connectivity
     *\param elements - Array of length num_elems containing elements
     *                  handles of elements for which connectivity is to
     *                  be queried.
     *\param vertices - Array of vertex handles in connectivity list.
     *\param offsets  - Indices into \ref vertex_handles, one per element
     */
    void MsqMOAB::elements_get_attached_vertices (
        const ElementHandle* elements,
        size_t num_elems,
        std::vector<VertexHandle>& vertices,
        std::vector<size_t>& offsets,
        MBMesquite::MsqError& err )
    {
        moab::ErrorCode ierr;
        assert ( sizeof ( moab::EntityHandle ) == sizeof ( EntityHandle ) );
        const moab::EntityHandle* elems = reinterpret_cast<const moab::EntityHandle*> ( elements );
        // get_adjacent_entities( elems, num_elems, moab::MBVERTEX, vertices, offsets, err ); MSQ_CHKERR(err);
        offsets.resize ( num_elems + 1 );
        offsets[0] = 0;
        vertices.clear();
        std::vector<moab::EntityHandle> mbverts;

        for ( unsigned ie = 0; ie < num_elems; ++ie )
        {
            const moab::EntityHandle* conn;
            int nconn;
            ierr = meshInstance->get_connectivity ( elems[ie], conn, nconn, false ); MB_CHK_ERR_RET ( ierr );
            offsets[ie + 1] = offsets[ie] + nconn;

            for ( int iconn = 0; iconn < nconn; ++iconn )
            {
                mbverts.push_back ( conn[iconn] );
            }
        }

        vertices.resize ( mbverts.size() );
        moab::EntityHandle* verts = reinterpret_cast<moab::EntityHandle*> ( arrptr ( vertices ) );

        for ( size_t iverts = 0; iverts < mbverts.size(); ++iverts )
        { verts[iverts] = mbverts[iverts]; }

        mbverts.clear();
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }


    void MsqMOAB::get_all_elements ( std::vector<ElementHandle>& elements,
                                     MsqError& err )
    {
        moab::ErrorCode ierr;

        if ( inputSetType == moab::MBMAXTYPE )
        {
            int num_vol, num_face;

            ierr = meshInstance->get_number_entities_by_dimension ( inputSet, 2, num_face, false ); MB_CHK_ERR_RET ( ierr );
            ierr = meshInstance->get_number_entities_by_dimension ( inputSet, 3, num_vol, false ); MB_CHK_ERR_RET ( ierr );
            elements.resize ( num_face + num_vol );

            if ( elements.empty() )
            { return; }

            moab::EntityHandle* ptr = reinterpret_cast<moab::EntityHandle*> ( arrptr ( elements ) );

            if ( num_face )
            {
                std::vector<moab::EntityHandle> faces;
                ierr = meshInstance->get_entities_by_dimension ( inputSet, 2, faces, false ); MB_CHK_ERR_RET ( ierr );
                assert ( faces.size() - num_face == 0 );
                std::copy ( faces.begin(), faces.end(), ptr );
            }

            if ( num_vol )
            {
                ptr += num_face;
                std::vector<moab::EntityHandle> regions;
                ierr = meshInstance->get_entities_by_dimension ( inputSet, 3, regions, false ); MB_CHK_ERR_RET ( ierr );
                assert ( regions.size() - num_vol == 0 );
                std::copy ( regions.begin(), regions.end(), ptr );
            }
        }
        else
        {
            int count;
            ierr = meshInstance->get_number_entities_by_type ( inputSet, inputSetType, count, false ); MB_CHK_ERR_RET ( ierr );

            if ( !count )
            { return; }

            elements.resize ( count );

            moab::EntityHandle* ptr = reinterpret_cast<moab::EntityHandle*> ( arrptr ( elements ) );
            std::vector<moab::EntityHandle> entities;
            ierr = meshInstance->get_entities_by_type ( inputSet, inputSetType, entities, false ); MB_CHK_ERR_RET ( ierr );
            assert ( entities.size() - count == 0 );
            std::copy ( entities.begin(), entities.end(), ptr );
        }

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    void MsqMOAB::get_all_vertices ( std::vector<VertexHandle>& vertices,
                                     MsqError& err )
    {
        std::vector<ElementHandle> elems;
        get_all_elements ( elems, err ); MSQ_CHKERR ( err );

        if ( elems.empty() )
        { return; }

        std::vector<size_t> offsets;
        elements_get_attached_vertices ( arrptr ( elems ), elems.size(), vertices, offsets, err ); MSQ_CHKERR ( err );

        std::sort ( vertices.begin(), vertices.end() );
        vertices.erase ( std::unique ( vertices.begin(), vertices.end() ), vertices.end() );
    }


    // Returns the topologies of the given entities.  The "entity_topologies"
    // array must be at least "num_elements" in size.
    void MsqMOAB::elements_get_topologies (
        const ElementHandle* element_handle_array,
        EntityTopology* element_topologies,
        size_t num_elements, MsqError& err )
    {
        if ( !num_elements )
        { return; }

        assert ( sizeof ( ElementHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( element_handle_array );

        for ( size_t i = 0; i < num_elements; ++i )
        {
            element_topologies[i] = topologyMap[ meshInstance->type_from_handle ( arr[i] ) ];
        }

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    //**************** Memory Management ****************
    // Tells the mesh that the client is finished with a given
    // entity handle.
    void MsqMOAB::release_entity_handles (
        const MBMesquite::Mesh::EntityHandle* /*handle_array*/,
        size_t /*num_handles*/, MsqError& /*err*/ )
    {
        // Do nothing
    }

    // Instead of deleting a Mesh when you think you are done,
    // call release().  In simple cases, the implementation could
    // just call the destructor.  More sophisticated implementations
    // may want to keep the Mesh object to live longer than Mesquite
    // is using it.
    void MsqMOAB::release()
    {
    }

    //**************** Tags ****************
    TagHandle MsqMOAB::tag_create ( const std::string& name,
                                    TagType type, unsigned length,
                                    const void*,
                                    MsqError& err )
    {
        moab::DataType moab_type;

        switch ( type )
        {
            case MBMesquite::Mesh::BYTE:   moab_type = moab::MB_TYPE_OPAQUE;   break;

            case MBMesquite::Mesh::INT:    moab_type = moab::MB_TYPE_INTEGER;  break;

            case MBMesquite::Mesh::DOUBLE: moab_type = moab::MB_TYPE_DOUBLE;   break;

            case MBMesquite::Mesh::HANDLE: moab_type = moab::MB_TYPE_HANDLE;   break;

            default:
                MSQ_SETERR ( err ) ( "Invalid tag type", MsqError::INVALID_ARG );
                return 0;
        }

        moab::ErrorCode ierr;
        moab::Tag result = 0;
        ierr = meshInstance->tag_get_handle ( name.c_str(), length, moab_type, result, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE ); MB_CHK_ERR_RET_VAL ( ierr, result );

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
        return static_cast<TagHandle> ( result );
    }

    void MsqMOAB::tag_delete ( TagHandle handle, MsqError& err )
    {
        moab::ErrorCode ierr = meshInstance->tag_delete ( static_cast<moab::Tag> ( handle ) ); MB_CHK_ERR_RET ( ierr );
        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    TagHandle MsqMOAB::tag_get ( const std::string& name, MsqError& err )
    {
        moab::ErrorCode ierr;
        moab::Tag handle = 0;
        ierr = meshInstance->tag_get_handle ( name.c_str(), handle ); MB_CHK_ERR_RET_VAL ( ierr, handle );

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
        return static_cast<TagHandle> ( handle );
    }


    void MsqMOAB::tag_properties ( TagHandle handle,
                                   std::string& name_out,
                                   TagType& type_out,
                                   unsigned& length_out,
                                   MsqError& err )
    {
        std::string buffer;
        moab::ErrorCode ierr;
        moab::DataType itype;
        int size;

        moab::Tag th = static_cast<moab::Tag> ( handle );
        ierr = meshInstance->tag_get_name ( th, buffer ); MB_CHK_ERR_RET ( ierr );
        ierr = meshInstance->tag_get_length ( th, size ); MB_CHK_ERR_RET ( ierr );
        ierr = meshInstance->tag_get_data_type ( th, itype ); MB_CHK_ERR_RET ( ierr );

        name_out = buffer;
        length_out = size;

        switch ( itype )
        {
            case moab::MB_TYPE_OPAQUE   : type_out = MBMesquite::Mesh::BYTE  ; break;

            case moab::MB_TYPE_INTEGER  : type_out = MBMesquite::Mesh::INT   ; break;

            case moab::MB_TYPE_DOUBLE   : type_out = MBMesquite::Mesh::DOUBLE; break;

            case moab::MB_TYPE_HANDLE   : type_out = MBMesquite::Mesh::HANDLE; break;

            default:
                MSQ_SETERR ( err ) ( "Unsupported iMesh tag type", MsqError::NOT_IMPLEMENTED );
                return;
        }

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

    void MsqMOAB::tag_set_element_data ( TagHandle tag,
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         const void* data,
                                         MsqError& err )
    {
        tag_set_data ( tag, num_elems, array, data, err );
    }

    void MsqMOAB::tag_set_vertex_data ( TagHandle tag,
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        const void* data,
                                        MsqError& err )
    {
        tag_set_data ( tag, num_elems, array, data, err );
    }

    void MsqMOAB::tag_set_data ( TagHandle tag,
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 const void* data,
                                 MsqError& err )
    {
        moab::ErrorCode ierr;
        int size;
        moab::Tag th = static_cast<moab::Tag> ( tag );
        ierr = meshInstance->tag_get_length ( th, size ); MB_CHK_ERR_RET ( ierr );

        assert ( sizeof ( EntityHandle ) == sizeof ( moab::EntityHandle ) );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( array );
        ierr = meshInstance->tag_set_data ( th, arr, num_elems, data ); MB_CHK_ERR_RET ( ierr );

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }


    void MsqMOAB::tag_get_element_data ( TagHandle tag,
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         void* data,
                                         MsqError& err )
    {
        tag_get_data ( tag, num_elems, array, data, err );
    }

    void MsqMOAB::tag_get_vertex_data ( TagHandle tag,
                                        size_t num_verts,
                                        const VertexHandle* array,
                                        void* data,
                                        MsqError& err )
    {
        tag_get_data ( tag, num_verts, array, data, err );
    }

    void MsqMOAB::tag_get_data ( TagHandle tag,
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 void* data,
                                 MsqError& err )
    {
        moab::ErrorCode ierr;

#if IMESH_VERSION_ATLEAST(1,1)
        void* ptr = data;
#else
        char* ptr = static_cast<char*> ( data );
#endif

        assert ( sizeof ( EntityHandle ) == sizeof ( moab::EntityHandle ) );
        moab::Tag th = static_cast<moab::Tag> ( tag );
        const moab::EntityHandle* arr = reinterpret_cast<const moab::EntityHandle*> ( array );
        ierr = meshInstance->tag_get_data ( th, arr, num_elems, ptr ); MB_CHK_ERR_RET ( ierr );

        err.set_error ( MBMesquite::MsqError::NO_ERROR );
    }

} // namespace Mesquite

