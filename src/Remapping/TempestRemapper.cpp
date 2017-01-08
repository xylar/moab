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
ErrorCode TempestRemapper::initialize()
{
	ErrorCode rval;
	rval = m_interface->create_meshset(moab::MESHSET_SET, m_source_set);MB_CHK_SET_ERR(rval, "Can't create new set");
	rval = m_interface->create_meshset(moab::MESHSET_SET, m_target_set);MB_CHK_SET_ERR(rval, "Can't create new set");
	rval = m_interface->create_meshset(moab::MESHSET_SET, m_overlap_set);MB_CHK_SET_ERR(rval, "Can't create new set");

	m_source = NULL;
	m_target = NULL;
	m_overlap = NULL;

	return MB_SUCCESS;
}

ErrorCode TempestRemapper::LoadTempestMesh_Private(std::string inputFilename, Mesh** tempest_mesh)
{
	if (TempestRemapper::verbose) std::cout << "\nLoading TempestRemap Mesh object from file = " << inputFilename << " ...\n";

	{
		NcError error(NcError::silent_nonfatal);

		try {
			// Load input mesh
			if (TempestRemapper::verbose) std::cout << "Loading mesh ...\n";
			Mesh* mesh = new Mesh(inputFilename);
			mesh->RemoveZeroEdges();
			if (TempestRemapper::verbose) std::cout << "----------------\n";

			// Validate mesh
			if (meshValidate) {
				if (TempestRemapper::verbose) std::cout << "Validating mesh ...\n";
				mesh->Validate();
				if (TempestRemapper::verbose) std::cout << "-------------------\n";
			}

			// Construct the edge map on the mesh
			if (constructEdgeMap) {
				if (TempestRemapper::verbose) std::cout << "Constructing edge map on mesh ...\n";
				mesh->ConstructEdgeMap();
				if (TempestRemapper::verbose) std::cout << "---------------------------------\n";
			}

			if (tempest_mesh) *tempest_mesh = mesh;

		} catch (Exception & e) {
			std::cout << "TempestRemap ERROR: " << e.ToString() << "\n";
			return MB_FAILURE;

		} catch (...) {
			return MB_FAILURE;
		}
	}
	return MB_SUCCESS;
}


ErrorCode TempestRemapper::LoadMesh(Remapper::IntersectionContext ctx, std::string inputFilename, TempestMeshType type)
{
	if (ctx == Remapper::SourceMesh) {
		m_source_type = type;
		return LoadTempestMesh_Private(inputFilename, &m_source);
	}
	else if (ctx == Remapper::TargetMesh) {
		m_target_type = type;
		return LoadTempestMesh_Private(inputFilename, &m_target);
	}
	else if (ctx != Remapper::DEFAULT) {
		m_overlap_type = type;
		return LoadTempestMesh_Private(inputFilename, &m_overlap);
	}
	else {
		MB_CHK_SET_ERR(MB_FAILURE, "Invalid IntersectionContext context provided");
	}
}


ErrorCode TempestRemapper::ConvertTempestMeshToMOAB_Private(TempestMeshType meshType, Mesh* mesh, EntityHandle& mesh_set)
{
	ErrorCode rval;

	if (TempestRemapper::verbose) std::cout << "\nConverting TempestRemap Mesh object to MOAB representation ...\n";
	const NodeVector& nodes = mesh->nodes;
	const FaceVector& faces = mesh->faces;

	// rval = m_interface->create_meshset(MESHSET_SET, mesh_set); MB_CHK_SET_ERR(rval, "Can't create new set");

	ReadUtilIface* iface;
	rval = m_interface->query_interface(iface); MB_CHK_SET_ERR(rval, "Can't get reader interface");

	Tag tmp_mb_loc_tag;
	int locid_def=-1;
	rval = m_interface->tag_get_handle("TEMPEST_MOAB_LOCALID", 1, MB_TYPE_INTEGER, tmp_mb_loc_tag, MB_TAG_DENSE | MB_TAG_CREAT, &locid_def);MB_CHK_ERR(rval);

	// Set the data for the vertices
	std::vector<double*> arrays;
	std::vector<int> loc_id(nodes.size());
	EntityHandle startv;
	rval = iface->get_node_coords(3, nodes.size(), 0, startv, arrays); MB_CHK_SET_ERR(rval, "Can't get node coords");
	for (unsigned iverts = 0; iverts < nodes.size(); ++iverts) {
		const Node& node = nodes[iverts];
		arrays[0][iverts] = node.x;
		arrays[1][iverts] = node.y;
		arrays[2][iverts] = node.z;
		loc_id[iverts] = iverts;
	}
	Range mbverts(startv, startv + nodes.size() - 1);
	m_interface->add_entities(mesh_set, mbverts);
	rval = m_interface->tag_set_data(tmp_mb_loc_tag, mbverts, &loc_id[0]);MB_CHK_ERR(rval);
	loc_id.clear();

	// We will assume all elements are of the same type - for now;
	// need a better way to categorize without doing a full pass first
	const unsigned lnum_v_per_elem = faces[0].edges.size(); // Linear elements: nedges = nverts ?
	if ((meshType != OVERLAP && meshType != OVERLAP_V2) && lnum_v_per_elem <= 4) {
		const unsigned num_v_per_elem = lnum_v_per_elem;
		EntityHandle starte; // Connectivity
		EntityHandle* conn;
		rval = iface->get_element_connect(faces.size(), num_v_per_elem, MBPOLYGON, 0, starte, conn); MB_CHK_SET_ERR(rval, "Can't get element connectivity");
		for (unsigned ifaces = 0, offset = 0; ifaces < faces.size(); ++ifaces) {
			const Face& face = faces[ifaces];
			conn[offset++] = startv + face.edges[0].node[1];
			for (unsigned iedges = 1; iedges < face.edges.size(); ++iedges) {
				conn[offset++] = startv + face.edges[iedges].node[1];
			}
			m_interface->add_entities(mesh_set, &starte, 1);
			starte++;
		}
	}
	else {
		if (TempestRemapper::verbose) std::cout << "..Mesh size: Nodes [" << nodes.size() << "] Elements [" << faces.size() << "].\n";
		const int NMAXPOLYEDGES = 15;
		std::vector<unsigned> nPolys(NMAXPOLYEDGES, 0);
		std::vector<std::vector<int> > typeNSeqs(NMAXPOLYEDGES);
		for (unsigned ifaces = 0; ifaces < faces.size(); ++ifaces) {
			const int iType = faces[ifaces].edges.size();
			nPolys[iType]++;
			typeNSeqs[iType].push_back(ifaces);
		}
		int iBlock = 0;
		for (unsigned iType = 0; iType < NMAXPOLYEDGES; ++iType) {
			if (!nPolys[iType]) continue; // Nothing to do

			if (TempestRemapper::verbose) std::cout << "....Block " << iBlock++ << " Polygons [" << iType << "] Elements [" << nPolys[iType] << "].\n";
			const unsigned num_v_per_elem = iType;
			EntityHandle starte; // Connectivity
			EntityHandle* conn;

			// Allocate the connectivity array, depending on the element type
			switch (num_v_per_elem) {
			case 3:
				rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBTRI, 0, starte, conn); MB_CHK_SET_ERR(rval, "Can't get element connectivity");
				break;
			case 4:
				rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBQUAD, 0, starte, conn); MB_CHK_SET_ERR(rval, "Can't get element connectivity");
				break;
			default:
				rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBPOLYGON, 0, starte, conn); MB_CHK_SET_ERR(rval, "Can't get element connectivity");
				break;
			}
			Range mbcells(starte, starte + nPolys[iType] - 1);
			m_interface->add_entities(mesh_set, mbcells);

			for (unsigned ifaces = 0, offset = 0; ifaces < typeNSeqs[iType].size(); ++ifaces) {
				const Face& face = faces[typeNSeqs[iType][ifaces]];
				conn[offset++] = startv + face.edges[0].node[1];
				for (unsigned iedges = 1; iedges < face.edges.size(); ++iedges) {
					conn[offset++] = startv + face.edges[iedges].node[1];
				}
			}

			// Now let us update the adjacency data, because some elements are new
			rval = iface->update_adjacencies(starte, nPolys[iType], num_v_per_elem, conn); MB_CHK_SET_ERR(rval, "Can't update adjacencies");
			// Generate all adj entities dimension 1 and 2 (edges and faces/ tri or qua)
			Range edges;
			rval = m_interface->get_adjacencies(mbcells, 1, true, edges,
			                           Interface::UNION); MB_CHK_SET_ERR(rval, "Can't get edges");
		}
	}

	return MB_SUCCESS;
}


ErrorCode TempestRemapper::ConvertTempestMesh(Remapper::IntersectionContext ctx)
{
	if (ctx == Remapper::SourceMesh) {
		return ConvertTempestMeshToMOAB_Private(m_source_type, m_source, m_source_set);
	}
	else if (ctx == Remapper::TargetMesh) {
		return ConvertTempestMeshToMOAB_Private(m_target_type, m_target, m_target_set);
	}
	else if (ctx != Remapper::DEFAULT) {
		return ConvertTempestMeshToMOAB_Private(m_overlap_type, m_overlap, m_overlap_set);
	}
	else {
		MB_CHK_SET_ERR(MB_FAILURE, "Invalid IntersectionContext context provided");
	}
}


ErrorCode TempestRemapper::ConvertMOABMeshToTempest_Private(Mesh* mesh, EntityHandle mesh_set)
{
	if (TempestRemapper::verbose) std::cout << "\nConverting MOAB Mesh object to TempestRemap Mesh representation ...\n";
	NodeVector& nodes = mesh->nodes;
	FaceVector& faces = mesh->faces;

	ErrorCode rval;
	Range verts, elems;
	rval = m_interface->get_entities_by_dimension(mesh_set, 2, elems); MB_CHK_ERR(rval);
	faces.resize(elems.size());

	rval = m_interface->get_adjacencies(elems, 0, true, verts, Interface::UNION); MB_CHK_ERR(rval);

	int iface = 0;
	for (Range::iterator ielems = elems.begin(); ielems != elems.end(); ++ielems, ++iface) {
		Face& face = faces[iface];

		// compute the number of edges per faces
		std::vector< EntityHandle > face_edges;
		rval = m_interface->get_adjacencies(&(*ielems), 1, 1, true, face_edges); MB_CHK_ERR(rval);
		face.edges.resize(face_edges.size());

		// get the connectivity for each edge
		const EntityHandle* connectface;
		int nnodesf;
		rval = m_interface->get_connectivity(*ielems, connectface, nnodesf); MB_CHK_ERR(rval);

		for (unsigned ifaces = 0; ifaces < face_edges.size(); ++ifaces) {
			face.SetNode( ifaces, verts.index(connectface[ifaces]) );
		}
	}

	int nnodes=verts.size();
	nodes.resize(nnodes);

	std::cout << "+-- Total vertices = " << nnodes << " and elements = " << elems.size() << std::endl;

	Tag tmp_mb_loc_tag;
	int locid_def=-1;
	rval = m_interface->tag_get_handle("TEMPEST_MOAB_LOCALID", 1, MB_TYPE_INTEGER, tmp_mb_loc_tag, MB_TAG_DENSE | MB_TAG_CREAT, &locid_def);MB_CHK_ERR(rval);

	std::vector<int> loc_id(nnodes);
	rval = m_interface->tag_get_data(tmp_mb_loc_tag, verts, &loc_id[0]);MB_CHK_ERR(rval);
	loc_id.clear();

	// Set the data for the vertices
	std::vector<double> coordx(nnodes), coordy(nnodes), coordz(nnodes);
	rval = m_interface->get_coords(verts, &coordx[0], &coordy[0], &coordz[0]); MB_CHK_ERR(rval);
	for (int inode=0; inode < nnodes; ++inode) {
		Node& node = nodes[inode];
		node.x = coordx[inode];
		node.y = coordy[inode];
		node.z = coordz[inode];
	}
	coordx.clear();
	coordy.clear();
	coordz.clear();

	// mesh->RemoveZeroEdges();
	// mesh->RemoveCoincidentNodes();
	mesh->Validate();
	return MB_SUCCESS;
}

// Should be ordered as Source, Target, Overlap
ErrorCode TempestRemapper::AssociateSrcTargetInOverlap(Mesh* mesh, EntityHandle* meshsets)
{
	ErrorCode rval;
	Tag redPtag,bluePtag;
	rval = m_interface->tag_get_handle("RedParent", redPtag);MB_CHK_ERR(rval);
	rval = m_interface->tag_get_handle("BlueParent", bluePtag);MB_CHK_ERR(rval);

	Range redEls, blueEls;
	rval = m_interface->get_entities_by_dimension(meshsets[0], 2, redEls); MB_CHK_ERR(rval);
	rval = m_interface->get_entities_by_dimension(meshsets[1], 2, blueEls); MB_CHK_ERR(rval);

	// Overlap mesh: mesh[2]
	mesh[2].vecSourceFaceIx.resize(redEls.size());
	mesh[2].vecTargetFaceIx.resize(blueEls.size());

	return MB_SUCCESS;
}


ErrorCode TempestRemapper::ExchangeGhostWeights(OfflineMap* /*weightMap*/)
{
	return MB_SUCCESS;
}


ErrorCode TempestRemapper::ConvertMeshToTempest(Remapper::IntersectionContext ctx)
{
	ErrorCode rval;

	if (ctx == Remapper::SourceMesh) {
		if (!m_source) m_source = new Mesh();
		rval = ConvertMOABMeshToTempest_Private(m_source, m_source_set);
	}
	else if (ctx == Remapper::TargetMesh) {
		if (!m_target) m_target = new Mesh();
		rval = ConvertMOABMeshToTempest_Private(m_target, m_target_set);
	}
	else if (ctx != Remapper::DEFAULT) { // Overlap mesh
		if (!m_overlap) m_overlap = new Mesh();
		rval = ConvertMOABMeshToTempest_Private(m_overlap, m_overlap_set);
	}
	else {
		MB_CHK_SET_ERR(MB_FAILURE, "Invalid IntersectionContext context provided");
	}

	return rval;
}

ErrorCode TempestRemapper::ComputeOverlapMesh(double tolerance, bool use_tempest)
{
	ErrorCode rval;
	// First, split based on whether to use Tempest or MOAB
	// If Tempest
	//   1) Check for valid Mesh and pointers to objects for source/target
	//   2) Invoke GenerateOverlapWithMeshes routine from Tempest library
	// If MOAB
	//   1) Check for valid source and target meshsets (and entities)
	//   2) 
	if (use_tempest) {
      // Now let us construct the overlap mesh, by calling TempestRemap interface directly
      // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
	  assert(m_source != NULL);
	  assert(m_target != NULL);
	  if (m_overlap != NULL) delete m_overlap;
      m_overlap = GenerateOverlapWithMeshes(*m_source, *m_target, "" /*outFilename*/, "exact", false);
	}
	else {
      const double radius = 1.0 /*2.0*acos(-1.0)*/;
      const double boxeps = 0.1;
	  // Create the intersection on the sphere object and set up necessary parameters
      moab::Range local_verts;
      moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(m_interface);
      mbintx->SetErrorTolerance(tolerance);
      mbintx->set_box_error(boxeps);
      mbintx->SetRadius(radius);

      rval = mbintx->FindMaxEdges(m_source_set, m_target_set);MB_CHK_ERR(rval);

      rval = mbintx->build_processor_euler_boxes(m_target_set, local_verts); MB_CHK_ERR(rval);
      
      moab::EntityHandle covering_set;
      // Note: lots of communication possible, if mesh is distributed very differently
      rval = mbintx->construct_covering_set(m_source_set, covering_set); MB_CHK_ERR(rval);

      // Now perform the actual parallel intersection between the source and the target meshes
      rval = mbintx->intersect_meshes(covering_set, m_target_set, m_overlap_set);MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");

      rval = m_interface->add_entities(m_overlap_set, &m_source_set, 1);MB_CHK_ERR(rval);
      rval = m_interface->add_entities(m_overlap_set, &m_target_set, 1);MB_CHK_ERR(rval);

      rval = fix_degenerate_quads(m_interface, m_overlap_set);MB_CHK_ERR(rval);
      rval = positive_orientation(m_interface, m_overlap_set, radius);MB_CHK_ERR(rval);

      // free the memory
      delete mbintx;
	}

	return MB_SUCCESS;
}


}
