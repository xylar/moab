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
#include "moab/TempestRemapper.hpp"
#include "moab/ReadUtilIface.hpp"

// Intersection includes
#include "moab/Intx2MeshOnSphere.hpp"
#include "moab/IntxUtils.hpp"

namespace moab
{

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::initialize()
{
    ErrorCode rval;
    rval = m_interface->create_meshset ( moab::MESHSET_SET, m_source_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
    rval = m_interface->create_meshset ( moab::MESHSET_SET, m_target_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
    rval = m_interface->create_meshset ( moab::MESHSET_SET, m_overlap_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );

    m_source = NULL;
    m_target = NULL;
    m_overlap = NULL;

    return MB_SUCCESS;
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
    const bool outputEnabled = ( TempestRemapper::verbose && !m_pcomm->rank() );
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
    if ( ctx == Remapper::SourceMesh )
    {
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (source) TempestRemap Mesh object to MOAB representation ...\n";
        return ConvertTempestMeshToMOAB_Private ( m_source_type, m_source, m_source_set );
    }
    else if ( ctx == Remapper::TargetMesh )
    {
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (target) TempestRemap Mesh object to MOAB representation ...\n";
        return ConvertTempestMeshToMOAB_Private ( m_target_type, m_target, m_target_set );
    }
    else if ( ctx != Remapper::DEFAULT )
    {
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (overlap) TempestRemap Mesh object to MOAB representation ...\n";
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
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "..Mesh size: Nodes [" << nodes.size() << "] Elements [" << faces.size() << "].\n";
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

            if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "....Block " << iBlock++ << " Polygons [" << iType << "] Elements [" << nPolys[iType] << "].\n";
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

    if ( ctx == Remapper::SourceMesh )
    {
        if ( !m_source ) m_source = new Mesh();
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (source) MOAB to TempestRemap Mesh representation ...\n";
        rval = ConvertMOABMeshToTempest_Private ( m_source, m_source_set, m_source_entities );
    }
    else if ( ctx == Remapper::TargetMesh )
    {
        if ( !m_target ) m_target = new Mesh();
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (target) MOAB to TempestRemap Mesh representation ...\n";
        rval = ConvertMOABMeshToTempest_Private ( m_target, m_target_set, m_target_entities );
    }
    else if ( ctx != Remapper::DEFAULT ) // Overlap mesh
    {
        if ( !m_overlap ) m_overlap = new Mesh();
        if ( TempestRemapper::verbose && !m_pcomm->rank() ) std::cout << "\nConverting (overlap) MOAB to TempestRemap Mesh representation ...\n";
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

    Range verts;
    rval = m_interface->get_entities_by_dimension ( mesh_set, 2, elems ); MB_CHK_ERR ( rval );

    faces.resize ( elems.size() );
    // rval = m_interface->get_adjacencies(elems, 0, true, verts, Interface::UNION); MB_CHK_ERR(rval);
    // rval = m_interface->get_entities_by_dimension(mesh_set, 0, verts); MB_CHK_ERR(rval);
    // verts.intersect(vertsall);

    for ( unsigned iface = 0; iface < elems.size(); ++iface )
    {

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity ( elems[iface], connectface, nnodesf ); MB_CHK_ERR ( rval );

        for ( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            if ( verts.index ( connectface[iverts] ) < 0 )
                verts.insert ( connectface[iverts] );
        }
    }

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
        // assert(face_edges.size() - nnodesf == 0);

        for ( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            int indx = verts.index ( connectface[iverts] );
            assert ( indx >= 0 );
            face.SetNode ( iverts, indx );
        }
    }
    // std::cout << "+-- Elements = " << elems.size() << " and verts = " << verts.size() << "\n"; std::cout.flush();

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

        rval = m_pcomm->exchange_tags ( gidtag, m_covering_source_entities ); MB_CHK_ERR ( rval );

        std::vector<int> gids ( m_covering_source_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_covering_source_entities, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_covsrc[gids[ie]] = ie;
            lid_to_gid_covsrc[ie] = gids[ie];
            // if(!m_pcomm->rank()) std::cout << "Src[" << ie << "]: GID = " << gids[ie] << " and value = " << gid_to_lid_src[gids[ie]] << "\n";
        }

        gids.resize ( m_source_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_source_entities, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_src[gids[ie]] = ie;
            lid_to_gid_src[ie] = gids[ie];
            // if(!m_pcomm->rank()) std::cout << "Src[" << ie << "]: GID = " << gids[ie] << " and value = " << gid_to_lid_src[gids[ie]] << "\n";
        }

        gids.resize ( m_target_entities.size(), -1 );
        rval = m_interface->tag_get_data ( gidtag,  m_target_entities/*m_covering_source_set*/, &gids[0] ); MB_CHK_ERR ( rval );
        for ( unsigned ie = 0; ie < gids.size(); ++ie )
        {
            gid_to_lid_tgt[gids[ie]] = ie;
            lid_to_gid_tgt[ie] = gids[ie];
            // if(!m_pcomm->rank()) std::cout << "Target[" << ie << "]: GID = " << gids[ie] << " and value = " << gid_to_lid_tgt[gids[ie]] << "\n";
        }
    }

    Tag redPtag, bluePtag;
    rval = m_interface->tag_get_handle ( "RedParent", redPtag ); MB_CHK_ERR ( rval );
    rval = m_interface->tag_get_handle ( "BlueParent", bluePtag ); MB_CHK_ERR ( rval );

    // Overlap mesh: resize the source and target connection arrays
    m_overlap->vecSourceFaceIx.resize ( m_overlap_entities.size() );
    m_overlap->vecTargetFaceIx.resize ( m_overlap_entities.size() );

    std::vector<int> rbids_src ( m_overlap_entities.size() );
    rval = m_interface->tag_get_data ( bluePtag,  m_overlap_entities, &rbids_src[0] ); MB_CHK_ERR ( rval );
    for ( unsigned ie = 0; ie < m_overlap_entities.size(); ++ie )
    {
        m_overlap->vecSourceFaceIx[ie] = gid_to_lid_covsrc[rbids_src[ie]];
        // if(m_pcomm->rank()) std::cout << "Overlap vecSourceFaceIx[" << ie << "]: GID = " << rbids_src[ie] << " and value = " << m_overlap->vecSourceFaceIx[ie] << "\n";
    }
    std::vector<int> rbids_tgt ( m_overlap_entities.size() );
    rval = m_interface->tag_get_data ( redPtag,  m_overlap_entities, &rbids_tgt[0] ); MB_CHK_ERR ( rval );
    for ( unsigned ie = 0; ie < m_overlap_entities.size(); ++ie )
    {
        m_overlap->vecTargetFaceIx[ie] = gid_to_lid_tgt[rbids_tgt[ie]];
        // if(!m_pcomm->rank()) std::cout << "Overlap vecTargetFaceIx[" << ie << "]: GID = " << rbids_tgt[ie] << " and value = " << m_overlap->vecTargetFaceIx[ie] << "\n";
    }

    if ( true )
    {
        // Let us re-sort the entities based on the vecSourceFaceIx values
        m_sorted_overlap_order.resize ( m_overlap_entities.size() );
        std::vector<std::pair<int, int> > unsorted_order ( m_overlap_entities.size() );
        for ( size_t ix = 0; ix < m_overlap_entities.size(); ++ix )
        {
            m_sorted_overlap_order[ix].first = m_overlap->vecSourceFaceIx[ix];
            m_sorted_overlap_order[ix].second = ix;

            unsorted_order[ix].first = m_overlap->vecSourceFaceIx[ix];
            unsorted_order[ix].second = ix;
        }

        std::sort ( m_sorted_overlap_order.begin(), m_sorted_overlap_order.end(), IntPairComparator );

        // if (m_pcomm->rank()) {
        // 	std::cout << "New Sorted Order:" << std::endl;
        //  for (size_t ix=0; ix < m_overlap_entities.size(); ++ix) {
        //  	std::cout << "[" << ix << "]:" << m_sorted_overlap_order[ix].first << " = " << m_sorted_overlap_order[ix].second ;
        //  	std::cout << " -- Unsorted :" << unsorted_order[ix].first << " = " << unsorted_order[ix].second << std::endl;
        //  }
        // }

        m_sorted_overlap = new Mesh();
        m_sorted_overlap->type = m_overlap->type; // propagate the type
        m_sorted_overlap->nodes.resize ( m_overlap->nodes.size() ); // propagate the nodes sizes
        std::copy ( m_overlap->nodes.begin(), m_overlap->nodes.end(), m_sorted_overlap->nodes.begin() ); // copy the nodes as is (unsorted)
        m_sorted_overlap->faces.resize ( m_overlap->faces.size() ); // propagate the face sizes
        m_sorted_overlap->vecSourceFaceIx.resize ( m_overlap->vecSourceFaceIx.size() );
        m_sorted_overlap->vecTargetFaceIx.resize ( m_overlap->vecTargetFaceIx.size() );
//      m_sorted_overlap->vecFaceArea.Initialize(m_overlap->vecFaceArea.GetRows());
        for ( size_t ifac = 0; ifac < m_overlap_entities.size(); ++ifac )
        {
            int tmp;
            tmp = gid_to_lid_covsrc[ifac];
            gid_to_lid_covsrc[rbids_src[ifac]] = gid_to_lid_covsrc[rbids_src[m_sorted_overlap_order[ifac].second]];
            gid_to_lid_covsrc[rbids_src[m_sorted_overlap_order[ifac].second]] = tmp;
            tmp = gid_to_lid_tgt[rbids_tgt[ifac]];
            gid_to_lid_tgt[rbids_tgt[ifac]] = gid_to_lid_tgt[rbids_tgt[m_sorted_overlap_order[ifac].second]];
            gid_to_lid_tgt[rbids_tgt[m_sorted_overlap_order[ifac].second]] = tmp;
            // std::swap(gid_to_lid_covsrc[[ifac]], gid_to_lid_covsrc[m_sorted_overlap_order[rbids_tgt[ifac]].second]);
            // std::swap(gid_to_lid_tgt[rbids_tgt[ifac]], gid_to_lid_tgt[m_sorted_overlap_order[rbids_tgt[ifac]].second]);
            m_sorted_overlap->vecSourceFaceIx[ifac] = m_overlap->vecSourceFaceIx[m_sorted_overlap_order[ifac].second];
            m_sorted_overlap->vecTargetFaceIx[ifac] = m_overlap->vecTargetFaceIx[m_sorted_overlap_order[ifac].second];
            // m_sorted_overlap->vecFaceArea[m_sorted_overlap_order[ifac].second] = m_overlap->vecFaceArea[ifac];
//      	std::cout << "Area [" << ifac << "]:" << m_overlap->vecFaceArea[ifac] << std::endl;
            Face& sface = m_sorted_overlap->faces[ifac];
            sface.edges.resize ( m_overlap->faces[ifac].edges.size() );
            for ( size_t iedg = 0; iedg < sface.edges.size(); ++iedg )
            {
                const Face& uface = m_overlap->faces[m_sorted_overlap_order[ifac].second];
                sface.edges[iedg].node[0] = uface.edges[iedg].node[0];
                sface.edges[iedg].node[1] = uface.edges[iedg].node[1];
                sface.edges[iedg].type = uface.edges[iedg].type;
            }
            // std::copy(m_overlap->faces[ifac].edges.begin(), m_overlap->faces[ifac].edges.end(), sface.edges.begin()); // copy the nodes as is (unsorted)
        }
        // std::cout << "New SrcFaceIx Order:" << std::endl;
        // for (size_t ix=0; ix < m_overlap_entities.size(); ++ix) {
        // 	std::cout << "[" << ix << "]:" << m_overlap->vecSourceFaceIx[ix] << " = " << m_sorted_overlap->vecSourceFaceIx[ix] << std::endl;
        // }

        delete m_overlap;
        m_overlap = m_sorted_overlap;
        m_overlap->ConstructEdgeMap();
        m_overlap->CalculateFaceAreas();
    }
    rbids_src.clear();
    rbids_tgt.clear();

    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

ErrorCode TempestRemapper::ComputeOverlapMesh ( double tolerance, double radius, bool use_tempest )
{
    ErrorCode rval;
    // First, split based on whether to use Tempest or MOAB
    // If Tempest
    //   1) Check for valid Mesh and pointers to objects for source/target
    //   2) Invoke GenerateOverlapWithMeshes routine from Tempest library
    // If MOAB
    //   1) Check for valid source and target meshsets (and entities)
    //   2)
    if ( use_tempest )
    {
        // Now let us construct the overlap mesh, by calling TempestRemap interface directly
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert ( m_source != NULL );
        assert ( m_target != NULL );
        if ( m_overlap != NULL ) delete m_overlap;
        m_overlap = GenerateOverlapWithMeshes ( *m_source, *m_target, "" /*outFilename*/, "exact", false );
    }
    else
    {
        // const double radius = 1.0 /*2.0*acos(-1.0)*/;
        const double boxeps = 0.1;
        // Create the intersection on the sphere object and set up necessary parameters
        moab::Range local_verts;
        moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere ( m_interface );

        mbintx->SetErrorTolerance ( tolerance );
        mbintx->set_box_error ( boxeps );
        mbintx->SetRadius ( radius );
        mbintx->set_parallel_comm ( m_pcomm );

        rval = mbintx->FindMaxEdges ( m_source_set, m_target_set ); MB_CHK_ERR ( rval );

        // Note: lots of communication possible, if mesh is distributed very differently
        if ( m_pcomm->size() != 1 )
        {
            rval = mbintx->build_processor_euler_boxes ( m_target_set, local_verts ); MB_CHK_ERR ( rval );

            rval = m_interface->create_meshset ( moab::MESHSET_SET, m_covering_source_set ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            rval = mbintx->construct_covering_set ( m_source_set, m_covering_source_set ); MB_CHK_ERR ( rval );

            m_covering_source = new Mesh();
            rval = ConvertMOABMeshToTempest_Private ( m_covering_source, m_covering_source_set, m_covering_source_entities ); MB_CHK_SET_ERR ( rval, "Can't convert source Tempest mesh" );

            /* Global ID - exchange data for covering data */
            Tag gidTag;
            rval = m_interface->tag_get_handle ( "GLOBAL_ID", gidTag ); MB_CHK_ERR ( rval );
            rval = m_pcomm->exchange_tags ( gidTag, m_covering_source_entities ); MB_CHK_ERR ( rval );

            {
                std::vector<int> gids ( m_covering_source_entities.size(), -1 );
                rval = m_interface->tag_get_data ( gidTag,  m_covering_source_entities, &gids[0] ); MB_CHK_ERR ( rval );
#if 0
                MPI_Barrier ( m_pcomm->comm() );
                if ( !m_pcomm->rank() )
                {
                    for ( unsigned ie = 0; ie < gids.size(); ++ie )
                        std::cout << "[Proc0] CoveringSrc[" << ie << "]: GID = " << gids[ie] << "\n";
                }
                MPI_Barrier ( m_pcomm->comm() );
                if ( m_pcomm->rank() )
                {
                    for ( unsigned ie = 0; ie < gids.size(); ++ie )
                        std::cout << "[Proc1] CoveringSrc[" << ie << "]: GID = " << gids[ie] << "\n";
                }
                MPI_Barrier ( m_pcomm->comm() );
#endif
            }

            m_intersecting_target_entities = moab::intersect ( m_source_entities, m_covering_source_entities );
            std::cout << "Number of common entities = " << m_intersecting_target_entities.size() << "/ (" << m_source_entities.size() << ", " << m_covering_source_entities.size() << ")\n";

        }
        else
        {
            m_covering_source_set = m_source_set;
            m_covering_source = m_source;
            m_covering_source_entities = m_source_entities; // this is a tempest mesh object; careful about incrementing the reference?
            m_intersecting_target_entities = m_source_entities; // no migration needed; source is completely covering target
        }

        // Now perform the actual parallel intersection between the source and the target meshes
        // rval = mbintx->intersect_meshes(m_source_set, m_covering_target_set, m_overlap_set);MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
        rval = mbintx->intersect_meshes ( m_covering_source_set, m_target_set, m_overlap_set ); MB_CHK_SET_ERR ( rval, "Can't compute the intersection of meshes on the sphere" );

        // rval = m_interface->add_entities(m_overlap_set, &m_source_set, 1);MB_CHK_ERR(rval);
        // rval = m_interface->add_entities(m_overlap_set, &m_target_set, 1);MB_CHK_ERR(rval);

        // Not needed
        /*
        rval = fix_degenerate_quads(m_interface, m_overlap_set);MB_CHK_ERR(rval);
        rval = positive_orientation(m_interface, m_overlap_set, radius);MB_CHK_ERR(rval);
        */

        // free the memory
        delete mbintx;
    }

    return MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

}
