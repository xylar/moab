/*
 * =====================================================================================
 *
 *       Filename:  TempestRemapper.hpp
 *
 *    Description:  Interface to the TempestRemap library to enable intersection and
 *                  high-order conservative remapping of climate solution from
 *                  arbitrary resolution of source and target grids on the sphere.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#include <string>
#include <iostream>
#include <cassert>
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/ReadUtilIface.hpp"

// Intersection includes
#include "moab/IntxMesh/Intx2MeshOnSphere.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"

namespace moab
{

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::initialize(bool initialize_fsets)
{
    ErrorCode rval;
    if (initialize_fsets) {
		rval = m_interface->create_meshset ( moab::MESHSET_SET, m_source_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
		rval = m_interface->create_meshset ( moab::MESHSET_SET, m_target_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
		rval = m_interface->create_meshset ( moab::MESHSET_SET, m_overlap_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
	}
	else {
		m_source_set = 0;
		m_target_set = 0;
		m_overlap_set = 0;
	}

    m_source = NULL;
    m_target = NULL;
    m_overlap = NULL;
    m_covering_source = NULL;

    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

TempestRemapper::~TempestRemapper()
{
    // destroy all meshes
    if ( m_source ) delete m_source;
    if ( m_target ) delete m_target;
    if ( m_overlap ) delete m_overlap;
    if ( m_covering_source && m_pcomm->size() != 1) delete m_covering_source;

    m_source_entities.clear();
    m_target_entities.clear();
    m_overlap_entities.clear();
    m_intersecting_target_entities.clear();
    gid_to_lid_src.clear(); gid_to_lid_tgt.clear(); gid_to_lid_covsrc.clear();
    lid_to_gid_src.clear(); lid_to_gid_tgt.clear(); lid_to_gid_covsrc.clear();
}

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::LoadMesh ( Remapper::IntersectionContext ctx, std::string inputFilename, TempestMeshType type )
{
    if ( ctx == Remapper::SourceMesh )
    {
        m_source_type = type;
        return LoadTempestMesh_Private ( inputFilename, &m_source );
    }
    else if ( ctx == Remapper::TargetMesh )
    {
        m_target_type = type;
        return LoadTempestMesh_Private ( inputFilename, &m_target );
    }
    else if ( ctx != Remapper::DEFAULT )
    {
        m_overlap_type = type;
        return LoadTempestMesh_Private ( inputFilename, &m_overlap );
    }
    else
    {
        MB_CHK_SET_ERR ( MB_FAILURE, "Invalid IntersectionContext context provided" );
    }
}

ErrorCode TempestRemapper::LoadTempestMesh_Private ( std::string inputFilename, Mesh** tempest_mesh )
{
    const bool outputEnabled = ( TempestRemapper::verbose && ((!m_pcomm) || !m_pcomm->rank()) );
    if ( outputEnabled ) std::cout << "\nLoading TempestRemap Mesh object from file = " << inputFilename << " ...\n";

    {
        NcError error ( NcError::silent_nonfatal );

        try
        {
            // Load input mesh
            if ( outputEnabled ) std::cout << "Loading mesh ...\n";
            Mesh* mesh = new Mesh ( inputFilename );
            mesh->RemoveZeroEdges();
            if ( outputEnabled ) std::cout << "----------------\n";

            // Validate mesh
            if ( meshValidate )
            {
                if ( outputEnabled ) std::cout << "Validating mesh ...\n";
                mesh->Validate();
                if ( outputEnabled ) std::cout << "-------------------\n";
            }

            // Construct the edge map on the mesh
            if ( constructEdgeMap )
            {
                if ( outputEnabled ) std::cout << "Constructing edge map on mesh ...\n";
                mesh->ConstructEdgeMap();
                if ( outputEnabled ) std::cout << "---------------------------------\n";
            }

            if ( tempest_mesh ) *tempest_mesh = mesh;

        }
        catch ( Exception & e )
        {
            std::cout << "TempestRemap ERROR: " << e.ToString() << "\n";
            return MB_FAILURE;

        }
        catch ( ... )
        {
            return MB_FAILURE;
        }
    }
    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::ConvertTempestMesh ( Remapper::IntersectionContext ctx )
{
    const bool outputEnabled = ( TempestRemapper::verbose && ((!m_pcomm) || !m_pcomm->rank()) );
    if ( ctx == Remapper::SourceMesh )
    {
        if ( outputEnabled ) std::cout << "\nConverting (source) TempestRemap Mesh object to MOAB representation ...\n";
        return ConvertTempestMeshToMOAB_Private ( m_source_type, m_source, m_source_set );
    }
    else if ( ctx == Remapper::TargetMesh )
    {
        if ( outputEnabled ) std::cout << "\nConverting (target) TempestRemap Mesh object to MOAB representation ...\n";
        return ConvertTempestMeshToMOAB_Private ( m_target_type, m_target, m_target_set );
    }
    else if ( ctx != Remapper::DEFAULT )
    {
        if ( outputEnabled ) std::cout << "\nConverting (overlap) TempestRemap Mesh object to MOAB representation ...\n";
        return ConvertTempestMeshToMOAB_Private ( m_overlap_type, m_overlap, m_overlap_set );
    }
    else
    {
        MB_CHK_SET_ERR ( MB_FAILURE, "Invalid IntersectionContext context provided" );
    }
}


ErrorCode TempestRemapper::ConvertTempestMeshToMOAB_Private ( TempestMeshType meshType, Mesh* mesh, EntityHandle& mesh_set )
{
    ErrorCode rval;

    const bool outputEnabled = ( TempestRemapper::verbose && ((!m_pcomm) || !m_pcomm->rank()) );
    const NodeVector& nodes = mesh->nodes;
    const FaceVector& faces = mesh->faces;

    ReadUtilIface* iface;
    rval = m_interface->query_interface ( iface ); MB_CHK_SET_ERR ( rval, "Can't get reader interface" );

    Tag tmp_mb_loc_tag;
    int locid_def = -1;
    rval = m_interface->tag_get_handle ( "TEMPEST_MOAB_LOCALID", 1, MB_TYPE_INTEGER, tmp_mb_loc_tag, MB_TAG_DENSE | MB_TAG_CREAT, &locid_def ); MB_CHK_ERR ( rval );

    // Set the data for the vertices
    std::vector<double*> arrays;
    std::vector<int> loc_id ( nodes.size() );
    EntityHandle startv;
    rval = iface->get_node_coords ( 3, nodes.size(), 0, startv, arrays ); MB_CHK_SET_ERR ( rval, "Can't get node coords" );
    for ( unsigned iverts = 0; iverts < nodes.size(); ++iverts )
    {
        const Node& node = nodes[iverts];
        arrays[0][iverts] = node.x;
        arrays[1][iverts] = node.y;
        arrays[2][iverts] = node.z;
        loc_id[iverts] = iverts;
    }
    Range mbverts ( startv, startv + nodes.size() - 1 );
    m_interface->add_entities ( mesh_set, mbverts );
    rval = m_interface->tag_set_data ( tmp_mb_loc_tag, mbverts, &loc_id[0] ); MB_CHK_ERR ( rval );
    loc_id.clear();

    // We will assume all elements are of the same type - for now;
    // need a better way to categorize without doing a full pass first
    const unsigned lnum_v_per_elem = faces[0].edges.size(); // Linear elements: nedges = nverts ?
    if ( ( meshType < OVERLAP_FILES ) && lnum_v_per_elem <= 4 )
    {
        const unsigned num_v_per_elem = lnum_v_per_elem;
        EntityHandle starte; // Connectivity
        EntityHandle* conn;
        rval = iface->get_element_connect ( faces.size(), num_v_per_elem, MBPOLYGON, 0, starte, conn ); MB_CHK_SET_ERR ( rval, "Can't get element connectivity" );
        for ( unsigned ifaces = 0, offset = 0; ifaces < faces.size(); ++ifaces )
        {
            const Face& face = faces[ifaces];
            conn[offset++] = startv + face.edges[0].node[1];
            for ( unsigned iedges = 1; iedges < face.edges.size(); ++iedges )
            {
                conn[offset++] = startv + face.edges[iedges].node[1];
            }
            m_interface->add_entities ( mesh_set, &starte, 1 );
            starte++;
        }
    }
    else
    {
        if ( outputEnabled ) std::cout << "..Mesh size: Nodes [" << nodes.size() << "] Elements [" << faces.size() << "].\n";
        const int NMAXPOLYEDGES = 15;
        std::vector<unsigned> nPolys ( NMAXPOLYEDGES, 0 );
        std::vector<std::vector<int> > typeNSeqs ( NMAXPOLYEDGES );
        for ( unsigned ifaces = 0; ifaces < faces.size(); ++ifaces )
        {
            const int iType = faces[ifaces].edges.size();
            nPolys[iType]++;
            typeNSeqs[iType].push_back ( ifaces );
        }
        int iBlock = 0;
        for ( unsigned iType = 0; iType < NMAXPOLYEDGES; ++iType )
        {
            if ( !nPolys[iType] ) continue; // Nothing to do

            if ( outputEnabled ) std::cout << "....Block " << iBlock++ << " Polygons [" << iType << "] Elements [" << nPolys[iType] << "].\n";
            const unsigned num_v_per_elem = iType;
            EntityHandle starte; // Connectivity
            EntityHandle* conn;

            // Allocate the connectivity array, depending on the element type
            switch ( num_v_per_elem )
            {
            case 3:
                rval = iface->get_element_connect ( nPolys[iType], num_v_per_elem, MBTRI, 0, starte, conn ); MB_CHK_SET_ERR ( rval, "Can't get element connectivity" );
                break;
            case 4:
                rval = iface->get_element_connect ( nPolys[iType], num_v_per_elem, MBQUAD, 0, starte, conn ); MB_CHK_SET_ERR ( rval, "Can't get element connectivity" );
                break;
            default:
                rval = iface->get_element_connect ( nPolys[iType], num_v_per_elem, MBPOLYGON, 0, starte, conn ); MB_CHK_SET_ERR ( rval, "Can't get element connectivity" );
                break;
            }
            Range mbcells ( starte, starte + nPolys[iType] - 1 );
            m_interface->add_entities ( mesh_set, mbcells );

            for ( unsigned ifaces = 0, offset = 0; ifaces < typeNSeqs[iType].size(); ++ifaces )
            {
                const Face& face = faces[typeNSeqs[iType][ifaces]];
                conn[offset++] = startv + face.edges[0].node[1];
                for ( unsigned iedges = 1; iedges < face.edges.size(); ++iedges )
                {
                    conn[offset++] = startv + face.edges[iedges].node[1];
                }
            }

            // Now let us update the adjacency data, because some elements are new
            rval = iface->update_adjacencies ( starte, nPolys[iType], num_v_per_elem, conn ); MB_CHK_SET_ERR ( rval, "Can't update adjacencies" );
            // Generate all adj entities dimension 1 and 2 (edges and faces/ tri or qua)
            Range edges;
            rval = m_interface->get_adjacencies ( mbcells, 1, true, edges,
                                                  Interface::UNION ); MB_CHK_SET_ERR ( rval, "Can't get edges" );
        }
    }

    return MB_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::ConvertMeshToTempest ( Remapper::IntersectionContext ctx )
{
    ErrorCode rval;
    const bool outputEnabled = ( TempestRemapper::verbose && ((!m_pcomm) || !m_pcomm->rank()) );

    if ( ctx == Remapper::SourceMesh )
    {
        if ( !m_source ) m_source = new Mesh();
        if ( outputEnabled ) std::cout << "\nConverting (source) MOAB to TempestRemap Mesh representation ...\n";
        rval = ConvertMOABMeshToTempest_Private ( m_source, m_source_set, m_source_entities );
    }
    else if ( ctx == Remapper::TargetMesh )
    {
        if ( !m_target ) m_target = new Mesh();
        if ( outputEnabled ) std::cout << "\nConverting (target) MOAB to TempestRemap Mesh representation ...\n";
        rval = ConvertMOABMeshToTempest_Private ( m_target, m_target_set, m_target_entities );
    }
    else if ( ctx != Remapper::DEFAULT )     // Overlap mesh
    {
        if ( !m_overlap ) m_overlap = new Mesh();
        if ( outputEnabled ) std::cout << "\nConverting (overlap) MOAB to TempestRemap Mesh representation ...\n";
        rval = ConvertMOABMeshToTempest_Private ( m_overlap, m_overlap_set, m_overlap_entities );
    }
    else
    {
        MB_CHK_SET_ERR ( MB_FAILURE, "Invalid IntersectionContext context provided" );
    }

    return rval;
}


ErrorCode TempestRemapper::ConvertMOABMeshToTempest_Private ( Mesh* mesh, EntityHandle mesh_set, moab::Range& elems )
{
    ErrorCode rval;
    NodeVector& nodes = mesh->nodes;
    FaceVector& faces = mesh->faces;

    moab::Range verts;
    rval = m_interface->get_entities_by_dimension ( mesh_set, 2, elems ); MB_CHK_ERR ( rval );

    // resize the number of elements in Tempest mesh
    faces.resize ( elems.size() );

    // let us now get the vertices from all the elements
    rval = m_interface->get_connectivity ( elems, verts ); MB_CHK_ERR ( rval );

    for ( unsigned iface = 0; iface < elems.size(); ++iface )
    {
        Face& face = faces[iface];
        EntityHandle ehandle = elems[iface];

        // compute the number of edges per faces
        std::vector< EntityHandle > face_edges;
        rval = m_interface->get_adjacencies ( &ehandle, 1, 1, true, face_edges ); MB_CHK_ERR ( rval );
        face.edges.resize ( face_edges.size() );

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity ( ehandle, connectface, nnodesf ); MB_CHK_ERR ( rval );

        // Can be untrue for polygonal elements with mixed pentagons and hexagons
        // int nvtx;
        // if (face_edges.size() - nnodesf != 0) {
        //  nvtx = face_edges.size();
        // }

        for ( size_t iverts = 0; iverts < face_edges.size(); ++iverts )
        {
            int indx = verts.index ( connectface[iverts] );
            assert ( indx >= 0 );
            face.SetNode ( iverts, indx );
        }
    }

    unsigned nnodes = verts.size();
    nodes.resize ( nnodes );

    Tag tmp_mb_loc_tag;
    int locid_def = -1;
    rval = m_interface->tag_get_handle ( "TEMPEST_MOAB_LOCALID", 1, MB_TYPE_INTEGER, tmp_mb_loc_tag, MB_TAG_DENSE | MB_TAG_CREAT, &locid_def ); MB_CHK_ERR ( rval );

    std::vector<int> loc_id ( nnodes );
    rval = m_interface->tag_get_data ( tmp_mb_loc_tag, verts, &loc_id[0] ); MB_CHK_ERR ( rval );
    loc_id.clear();

    // Set the data for the vertices
    std::vector<double> coordx ( nnodes ), coordy ( nnodes ), coordz ( nnodes );
    rval = m_interface->get_coords ( verts, &coordx[0], &coordy[0], &coordz[0] ); MB_CHK_ERR ( rval );
    for ( unsigned inode = 0; inode < nnodes; ++inode )
    {
        Node& node = nodes[inode];
        node.x = coordx[inode];
        node.y = coordy[inode];
        node.z = coordz[inode];
    }
    coordx.clear();
    coordy.clear();
    coordz.clear();
    verts.clear();

    mesh->RemoveZeroEdges();
    mesh->RemoveCoincidentNodes();

    // Generate reverse node array and edge map
    if ( constructEdgeMap ) mesh->ConstructEdgeMap();
    mesh->ConstructReverseNodeArray();

    // mesh->Validate();
    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

bool IntPairComparator ( const std::pair<int, int> &a, const std::pair<int, int> &b )
{
    if ( a.first == b.first )
        return a.second < b.second;
    else
        return a.first < b.first;
}

ErrorCode TempestRemapper::ConvertMOABMesh_WithSortedEntitiesBySource()
{
    ErrorCode rval;

    rval = m_interface->get_entities_by_dimension ( m_overlap_set, 2, m_overlap_entities ); MB_CHK_ERR ( rval );

    // Allocate for the overlap mesh
    if ( !m_overlap ) m_overlap = new Mesh();

    std::vector<std::pair<int, int> > sorted_overlap_order ( m_overlap_entities.size() );
    {
        Tag bluePtag, redPtag;
        rval = m_interface->tag_get_handle ( "BlueParent", bluePtag ); MB_CHK_ERR ( rval );
        rval = m_interface->tag_get_handle ( "RedParent", redPtag ); MB_CHK_ERR ( rval );

        // Overlap mesh: resize the source and target connection arrays
        m_overlap->vecSourceFaceIx.resize ( m_overlap_entities.size() );
        m_overlap->vecTargetFaceIx.resize ( m_overlap_entities.size() );

        m_covering_overlap_flag.resize(m_covering_source_entities.size(), false);

        // Overlap mesh: resize the source and target connection arrays
        std::vector<int> rbids_src ( m_overlap_entities.size() ), rbids_tgt ( m_overlap_entities.size() );
        rval = m_interface->tag_get_data ( bluePtag,  m_overlap_entities, &rbids_src[0] ); MB_CHK_ERR ( rval );
        rval = m_interface->tag_get_data ( redPtag,  m_overlap_entities, &rbids_tgt[0] ); MB_CHK_ERR ( rval );

        // Let us re-sort the entities based on the vecSourceFaceIx values
        for ( size_t ix = 0; ix < m_overlap_entities.size(); ++ix )
        {
            sorted_overlap_order[ix].first = gid_to_lid_covsrc[rbids_src[ix]];
            sorted_overlap_order[ix].second = ix;
        }
        std::sort ( sorted_overlap_order.begin(), sorted_overlap_order.end(), IntPairComparator );

        unsigned participating_covering_elems = 0;
        for ( unsigned ie = 0; ie < m_overlap_entities.size(); ++ie )
        {
            m_overlap->vecSourceFaceIx[ie] = gid_to_lid_covsrc[rbids_src[sorted_overlap_order[ie].second]];
            participating_covering_elems++;
            // if ( !m_covering_overlap_flag[rbids_src[ie]] ) {
            //     m_covering_overlap_flag[rbids_src[ie]] = true;
            //     participating_covering_elems++;
            // }

            m_overlap->vecTargetFaceIx[ie] = gid_to_lid_tgt[rbids_tgt[sorted_overlap_order[ie].second]];
            // if ( !m_pcomm->rank() ) printf ( "Element %i :: Src: [%i], Tgt: [%i]\n", ie, m_overlap->vecSourceFaceIx[ie], m_overlap->vecTargetFaceIx[ie] );
        }
        // printf ( "[%i] Total participating covering elements: %i\n", m_pcomm->rank(), participating_covering_elems );
    }

    FaceVector& faces = m_overlap->faces;
    faces.resize ( m_overlap_entities.size() );

    NodeVector& nodes = m_overlap->nodes;
    Range verts;
    for ( unsigned ifac = 0; ifac < m_overlap_entities.size(); ++ifac )
    {
        const unsigned iface = sorted_overlap_order[ifac].second;

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity ( m_overlap_entities[iface], connectface, nnodesf ); MB_CHK_ERR ( rval );

        for ( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            if ( verts.index ( connectface[iverts] ) < 0 )
                verts.insert ( connectface[iverts] );
        }
    }

      // compute the number of edges per faces
    moab::Range face_edges_exist, face_edges_all, face_edges_noexist;
    rval = m_interface->get_adjacencies ( m_overlap_entities, 1, false, face_edges_exist, moab::Interface::UNION); MB_CHK_ERR ( rval );
    rval = m_interface->get_adjacencies ( m_overlap_entities, 1, true, face_edges_all, moab::Interface::UNION); MB_CHK_ERR ( rval );
    face_edges_noexist = subtract(face_edges_all, face_edges_exist);

    for ( unsigned ifac = 0; ifac < m_overlap_entities.size(); ++ifac )
    {
        const unsigned iface = sorted_overlap_order[ifac].second;
        Face& face = faces[ifac];
        EntityHandle ehandle = m_overlap_entities[iface];

        // compute the number of edges per faces
        std::vector< EntityHandle > face_edges;
        rval = m_interface->get_adjacencies ( &ehandle, 1, 1, false, face_edges ); MB_CHK_ERR ( rval );
        face.edges.resize ( face_edges.size() );

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity ( ehandle, connectface, nnodesf ); MB_CHK_ERR ( rval );

        // Can be untrue for polygonal elements with mixed pentagons and hexagons
        // assert(face_edges.size() - nnodesf == 0);

        for ( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            int indx = verts.index ( connectface[iverts] );
            assert ( indx >= 0 );
            face.SetNode ( iverts, indx );
        }
    }

    rval = m_interface->delete_entities(face_edges_noexist);MB_CHK_ERR(rval);
    // rval = m_interface->add_entities(m_overlap_set,face_edges_noexist);MB_CHK_ERR(rval);

    unsigned nnodes = verts.size();
    nodes.resize ( nnodes );

    Tag tmp_mb_loc_tag;
    int locid_def = -1;
    rval = m_interface->tag_get_handle ( "TEMPEST_MOAB_LOCALID", 1, MB_TYPE_INTEGER, tmp_mb_loc_tag, MB_TAG_DENSE | MB_TAG_CREAT, &locid_def ); MB_CHK_ERR ( rval );

    std::vector<int> loc_id ( nnodes );
    rval = m_interface->tag_get_data ( tmp_mb_loc_tag, verts, &loc_id[0] ); MB_CHK_ERR ( rval );
    loc_id.clear();

    // Set the data for the vertices
    std::vector<double> coordx ( nnodes ), coordy ( nnodes ), coordz ( nnodes );
    rval = m_interface->get_coords ( verts, &coordx[0], &coordy[0], &coordz[0] ); MB_CHK_ERR ( rval );
    for ( unsigned inode = 0; inode < nnodes; ++inode )
    {
        Node& node = nodes[inode];
        node.x = coordx[inode];
        node.y = coordy[inode];
        node.z = coordz[inode];
    }
    coordx.clear();
    coordy.clear();
    coordz.clear();
    verts.clear();

    m_overlap->RemoveZeroEdges();
    m_overlap->RemoveCoincidentNodes();

    // Generate reverse node array and edge map
    if ( constructEdgeMap ) m_overlap->ConstructEdgeMap();
    m_overlap->ConstructReverseNodeArray();

    // m_overlap->Validate();
    return MB_SUCCESS;
}

// Should be ordered as Source, Target, Overlap
ErrorCode TempestRemapper::AssociateSrcTargetInOverlap()
{
    ErrorCode rval;

    gid_to_lid_src.clear(); lid_to_gid_src.clear();
    gid_to_lid_covsrc.clear(); lid_to_gid_covsrc.clear();
    gid_to_lid_tgt.clear(); lid_to_gid_tgt.clear();
    {
        Tag gidtag;
        rval = m_interface->tag_get_handle ( "GLOBAL_ID", gidtag ); MB_CHK_ERR ( rval );

        std::vector<int> gids ( m_covering_source_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_covering_source_entities, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_covsrc[gids[ie]] = ie;
            lid_to_gid_covsrc[ie] = gids[ie];
        }

        gids.resize ( m_source_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_source_entities, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_src[gids[ie]] = ie;
            lid_to_gid_src[ie] = gids[ie];
        }

        gids.resize ( m_target_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_target_entities/*m_covering_source_set*/, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_tgt[gids[ie]] = ie;
            lid_to_gid_tgt[ie] = gids[ie];
        }
    }

    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::ComputeOverlapMesh ( double tolerance, double radius_src, double radius_tgt, double boxeps, bool use_tempest )
{
    ErrorCode rval;
    // First, split based on whether to use Tempest or MOAB
    // If Tempest
    //   1) Check for valid Mesh and pointers to objects for source/target
    //   2) Invoke GenerateOverlapWithMeshes routine from Tempest library
    // If MOAB
    //   1) Check for valid source and target meshsets (and entities)
    //   2) Build processor bounding boxes and construct a covering set
    //   3) Perform intersection between the source (covering) and target entities 
    if ( use_tempest )
    {
        // Now let us construct the overlap mesh, by calling TempestRemap interface directly
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert ( m_source != NULL );
        assert ( m_target != NULL );
        if ( m_overlap != NULL ) delete m_overlap;
        m_overlap = new Mesh();
        bool concaveMeshA=false, concaveMeshB=false;
        int err = GenerateOverlapWithMeshes ( *m_source, *m_target, *m_overlap, "" /*outFilename*/, "exact", concaveMeshA, concaveMeshB, false );
        if (err) {
            rval = MB_FAILURE;
            return rval;
        }
    }
    else
    {
        // const double radius = 1.0 /*2.0*acos(-1.0)*/;
        // const double boxeps = 0.1;
        // Create the intersection on the sphere object and set up necessary parameters
        moab::Range local_verts;
        moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere ( m_interface );

        mbintx->set_error_tolerance ( tolerance );
        mbintx->set_radius_source_mesh ( radius_src );
        mbintx->set_radius_destination_mesh ( radius_tgt );
        mbintx->set_box_error ( boxeps );
        mbintx->set_parallel_comm ( m_pcomm );

        rval = mbintx->FindMaxEdges ( m_source_set, m_target_set ); MB_CHK_ERR ( rval );

        // Note: lots of communication possible, if mesh is distributed very differently
        if ( m_pcomm->size() != 1 )
        {
            rval = mbintx->build_processor_euler_boxes ( m_target_set, local_verts ); MB_CHK_ERR ( rval );

            rval = m_interface->create_meshset ( moab::MESHSET_SET, m_covering_source_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            rval = mbintx->construct_covering_set ( m_source_set, m_covering_source_set ); MB_CHK_ERR ( rval );
        }
        else
        {
            m_covering_source_set = m_source_set;
            m_covering_source = m_source;
            m_covering_source_entities = m_source_entities; // this is a tempest mesh object; careful about incrementing the reference?
            m_intersecting_target_entities = m_source_entities; // no migration needed; source is completely covering target
        }

        // Now perform the actual parallel intersection between the source and the target meshes
        rval = mbintx->intersect_meshes ( m_covering_source_set, m_target_set, m_overlap_set ); MB_CHK_SET_ERR ( rval, "Can't compute the intersection of meshes on the sphere" );

        if (m_pcomm->size() > 1) {

            // because we do not want to work with elements in coverage set that do not participate in intersection,
            // remove them from the coverage set
            // we will not delete them yet, just remove from the set !
            Range intxCov;
            Range intxCells;
            Tag blueParentHandleTag;
            rval = m_interface->tag_get_handle("BlueParent", blueParentHandleTag);  MB_CHK_ERR ( rval );
            rval = m_interface->get_entities_by_dimension(m_overlap_set, 2, intxCells);  MB_CHK_ERR ( rval );
            for (Range::iterator it=intxCells.begin(); it!=intxCells.end(); it++)
            {
              EntityHandle intxCell= *it;
              // EntityHandle blueParent;
              int blueParent;
              rval = m_interface->tag_get_data(blueParentHandleTag, &intxCell, 1, &blueParent ); MB_CHK_ERR ( rval );
              intxCov.insert(m_covering_source_entities[gid_to_lid_covsrc[blueParent]]);
            }

            Range notNeededCovCells = moab::subtract(m_covering_source_entities, intxCov);
            // remove now from coverage set the cells that are not needed
            // rval = m_interface->remove_entities(m_covering_source_set, notNeededCovCells); MB_CHK_ERR ( rval );
            m_covering_source_entities = moab::subtract(m_covering_source_entities, notNeededCovCells);
            m_intersecting_target_entities = moab::intersect ( m_source_entities, m_covering_source_entities );
// #ifdef VERBOSE
            std::cout << " remove from coverage set elements that are not intersected: " << notNeededCovCells.size() << "\n";
// #endif

            m_covering_source = new Mesh();
            rval = ConvertMOABMeshToTempest_Private ( m_covering_source, m_covering_source_set, m_covering_source_entities ); MB_CHK_SET_ERR ( rval, "Can't convert source Tempest mesh" );
        }

        // Fix any inconsistencies in the overlap mesh
        rval = fix_degenerate_quads(m_interface, m_overlap_set);MB_CHK_ERR(rval);
        rval = positive_orientation(m_interface, m_overlap_set, radius_tgt);MB_CHK_ERR(rval);

        // Now let us re-convert the MOAB mesh back to Tempest representation
        rval = this->AssociateSrcTargetInOverlap();MB_CHK_ERR(rval);
        rval = this->ConvertMOABMesh_WithSortedEntitiesBySource();MB_CHK_ERR(rval);

        // free the memory
        delete mbintx;
    }

    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

}

