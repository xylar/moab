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
    void TempestRemapper::initialize()
    {
        Range mbverts, mbelems;
    }

    ErrorCode TempestRemapper::LoadTempestMesh(std::string inputFilename, Mesh** tempest_mesh, bool meshValidate, bool constructEdgeMap)
    {
    	std::cout << "\nLoading TempestRemap Mesh object from file = " << inputFilename << " ...\n";
	    {
	        NcError error(NcError::silent_nonfatal);

	        try {
	            // Load input mesh
	            std::cout << "Loading mesh ...\n";
	            Mesh* mesh = new Mesh(inputFilename);
	            mesh->RemoveZeroEdges();
	            std::cout << "----------------\n";

	            // Validate mesh
	            if (meshValidate) {
	                std::cout << "Validating mesh ...\n";
	                mesh->Validate();
	                std::cout << "-------------------\n";
	            }

	            // Construct the edge map on the mesh
	            if (constructEdgeMap) {
	                std::cout << "Constructing edge map on mesh ...\n";
	                mesh->ConstructEdgeMap();
	                std::cout << "---------------------------------\n";
	            }

	            if (tempest_mesh) *tempest_mesh = mesh;

	        } catch(Exception & e) {
	            std::cout << "TempestRemap ERROR: " << e.ToString() << "\n";
	            return MB_FAILURE;

	        } catch(...) {
	            return MB_FAILURE;
	        }
	    }
	    return MB_SUCCESS;
    }


    ErrorCode TempestRemapper::LoadMesh(IntersectionContext ctx, std::string inputFilename, TempestMeshType type)
    {
    	Mesh* mesh;
    	if (ctx == SourceMesh) {
    		mesh = m_source;
    		m_source_type = type;
    	}
    	else if (ctx == TargetMesh) {
    		mesh = m_target;
    		m_target_type = type;
    	}
    	else {
    		mesh = m_overlap;
    	}

    	return TempestRemapper::LoadTempestMesh(inputFilename, &mesh, meshValidate, constructEdgeMap);    	
    }


	ErrorCode TempestRemapper::ConvertTempestMeshToMOAB(TempestMeshType meshType, Interface* mb, Mesh* mesh, EntityHandle& mesh_set)
	{
	  ErrorCode rval;

	  std::cout << "\nConverting TempestRemap Mesh object to MOAB representation ...\n";
	  const NodeVector& nodes = mesh->nodes;
	  const FaceVector& faces = mesh->faces;

	  rval = mb->create_meshset(MESHSET_SET, mesh_set);MB_CHK_SET_ERR(rval, "Can't create new set");

	  ReadUtilIface* iface;
	  rval = mb->query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

	  // Set the data for the vertices
	  std::vector<double*> arrays;
	  EntityHandle startv;
	  rval = iface->get_node_coords(3, nodes.size(), 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");
	  for (unsigned iverts=0; iverts < nodes.size(); ++iverts) {
	      const Node& node = nodes[iverts];
	      arrays[0][iverts] = node.x;
	      arrays[1][iverts] = node.y;
	      arrays[2][iverts] = node.z;
	  }
	  Range mbverts(startv, startv + nodes.size() - 1);
	  mb->add_entities(mesh_set, mbverts);

	  // We will assume all elements are of the same type - for now;
	  // need a better way to categorize without doing a full pass first
	  const unsigned lnum_v_per_elem = faces[0].edges.size(); // Linear elements: nedges = nverts ?
	  if ((meshType != OVERLAP && meshType != OVERLAP_V2) && lnum_v_per_elem <= 4) {
	      const unsigned num_v_per_elem = lnum_v_per_elem;
	      EntityHandle starte; // Connectivity
	      EntityHandle* conn;
	      rval = iface->get_element_connect(faces.size(), num_v_per_elem, MBPOLYGON, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
	      for (unsigned ifaces=0,offset=0; ifaces < faces.size(); ++ifaces) {
	          const Face& face = faces[ifaces];
	          conn[offset++] = startv+face.edges[0].node[1];
	          for (unsigned iedges=1; iedges < face.edges.size(); ++iedges) {
	              conn[offset++] = startv+face.edges[iedges].node[1];
	          }
	      }
	  }
	  else {
	      std::cout << "..Mesh size: Nodes [" << nodes.size() << "] Elements [" << faces.size() << "].\n";
	      const int NMAXPOLYEDGES = 15;
	      std::vector<unsigned> nPolys(NMAXPOLYEDGES,0);
	      std::vector<std::vector<int> > typeNSeqs(NMAXPOLYEDGES);
	      for (unsigned ifaces=0; ifaces < faces.size(); ++ifaces) {
	          const int iType = faces[ifaces].edges.size();
	          nPolys[iType]++;
	          typeNSeqs[iType].push_back(ifaces);
	      }
	      int iBlock=0;
	      for (unsigned iType=0; iType < NMAXPOLYEDGES; ++iType) {
	          if (!nPolys[iType]) continue; // Nothing to do

	          std::cout << "....Block " << iBlock++ << " Polygons [" << iType << "] Elements [" << nPolys[iType] << "].\n";
	          const unsigned num_v_per_elem = iType;
	          EntityHandle starte; // Connectivity
	          EntityHandle* conn;

	          // Allocate the connectivity array, depending on the element type
	          switch(num_v_per_elem) {
	            case 3:
	              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBTRI, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
	              break;
	            case 4:
	              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
	              break;
	            default:
	              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, MBPOLYGON, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
	              break;
	          }
	          Range mbcells(starte, starte + nPolys[iType] - 1);
	          mb->add_entities(mesh_set, mbcells);

	          for (unsigned ifaces=0,offset=0; ifaces < typeNSeqs[iType].size(); ++ifaces) {
	              const Face& face = faces[typeNSeqs[iType][ifaces]];
	              conn[offset++] = startv+face.edges[0].node[1];
	              for (unsigned iedges=1; iedges < face.edges.size(); ++iedges) {
	                  conn[offset++] = startv+face.edges[iedges].node[1];
	              }
	          }

	          // Now let us update the adjacency data, because some elements are new
	          rval = iface->update_adjacencies(starte, nPolys[iType], num_v_per_elem, conn);MB_CHK_SET_ERR(rval, "Can't update adjacencies");
	          // Generate all adj entities dimension 1 and 2 (edges and faces/ tri or qua)
	          Range edges;
	          rval = mb->get_adjacencies(mbcells, 1, true, edges,
	                                     Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get edges");
	      }
	  }

	  return MB_SUCCESS;
	}


    ErrorCode TempestRemapper::ConvertTempestMesh(IntersectionContext ctx, EntityHandle& meshset)
	{
	  Mesh* mesh;
	  TempestMeshType type;

		if (ctx == SourceMesh) {
			mesh = m_source;
			type = m_source_type;
		}
		else if (ctx == TargetMesh) {
			mesh = m_target;
			type = m_target_type;
		}
		else {
			mesh = m_overlap;
			type = OVERLAP;
		}
	  
	  return TempestRemapper::ConvertTempestMeshToMOAB(type, m_interface, mesh, meshset);
	}


	ErrorCode TempestRemapper::ConvertMOABMeshToTempest(Interface* mb, Mesh* mesh, EntityHandle mesh_set)
	{
	  ErrorCode rval;

	  std::cout << "\nConverting MOAB Mesh object to TempestRemap Mesh representation ...\n";
	  NodeVector& nodes = mesh->nodes;
	  FaceVector& faces = mesh->faces;

	  Range verts;
	  rval = mb->get_entities_by_dimension(mesh_set, 0, verts);MB_CHK_ERR(rval);
	  nodes.resize(verts.size());

	  // Set the data for the vertices
	  int inode=0;
	  std::vector<double> coordx(verts.size()), coordy(verts.size()), coordz(verts.size());
	  rval = mb->get_coords(verts, &coordx[0], &coordy[0], &coordz[0]);MB_CHK_ERR(rval);
	  for (Range::iterator iverts=verts.begin(); iverts!=verts.end(); ++iverts) {
	      Node& node = nodes[inode];
	      node.x = coordx[inode];
	      node.y = coordy[inode];
	      node.z = coordz[inode];
	      inode++;
	  }
	  coordx.clear();
	  coordy.clear();
	  coordz.clear();

	  Range elems;
	  rval = mb->get_entities_by_dimension(mesh_set, 2, elems);MB_CHK_ERR(rval);
	  faces.resize(elems.size());

	  int iface=0;
	  for (Range::iterator ielems=elems.begin(); ielems!=elems.end(); ++ielems) {
	    Face& face = faces[iface++];

	    // compute the number of edges per faces
	    std::vector< EntityHandle > face_edges;
	    rval = mb->get_adjacencies(&(*ielems), 1, 1, true, face_edges);MB_CHK_ERR(rval);
	    face.edges.resize(face_edges.size());

	    for (unsigned iedges=0; iedges < face_edges.size(); ++iedges) {
	        Edge& edge = face.edges[iedges];

	        // get the connectivity for each edge
	        const EntityHandle* connect;
	        int nnodes;
	        rval = mb->get_connectivity(face_edges[iedges], connect, nnodes);MB_CHK_ERR(rval);

	        // we expect only linear edges (2 nodes/edge)
	        assert(nnodes == 2);

	        // assign the edge nodes
	        edge.node[0] = connect[0];
	        edge.node[1] = connect[1];
	    }
	  }

	  mesh->RemoveZeroEdges();
	  return MB_SUCCESS;
	}


    ErrorCode TempestRemapper::ConvertMeshToTempest(EntityHandle meshset, Mesh* mesh)
	{ 
	  return TempestRemapper::ConvertMOABMeshToTempest(m_interface, mesh, meshset);
	}


}
