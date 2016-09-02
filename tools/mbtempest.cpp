/*
 * Usage: MOAB-Tempest tool
 *
 * Generate a Cubed-Sphere mesh: ./mbtempest -t 0 -res 25 -f cubed_sphere_mesh.exo
 * Generate a RLL mesh: ./mbtempest -t 1 -res 25 -f rll_mesh.exo
 * Generate a Icosahedral-Sphere mesh: ./mbtempest -t 2 -res 25 <-dual> -f icosa_mesh.exo
 *
 * Now you can compute the intersections between the meshes too!
 *
 * Generate the overlap mesh: ./mbtempest -t 3 -l cubed_sphere_mesh.exo -l rll_mesh.exo -f overlap_mesh.exo
 *
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReadUtilIface.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

// Tempest includes
#include "TempestRemapAPI.h"

enum MeshType { CS=0, RLL=1, ICO=2, OVERLAP=3 };

struct ToolContext {
    int blockSize;
    std::vector<std::string> inFilenames;
    std::string outFilename;
    int meshType;
    bool computeDual;

    ToolContext() : blockSize(5), outFilename("output.exo"), meshType(CS), computeDual(false)
    {
        inFilenames.resize(2);
    }

    void ParseCLOptions(int argc, char* argv[]) {
        ProgOptions opts;
        std::string expectedFName="output.exo";

        opts.addOpt<int>("res,r", "Resolution of the mesh (default=5)", &blockSize);
        opts.addOpt<int>("type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP=3])", &meshType);
        opts.addOpt<std::string>("file,f", "Output mesh filename (default=output.exo)", &outFilename);
        opts.addOpt<void>("dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual);
        opts.addOpt<std::string>("load,l", "Input mesh filenames (a source and target mesh)", &expectedFName);

        opts.parseCommandLine(argc, argv);

        if (meshType == OVERLAP) {
          opts.getOptAllArgs("load,l", inFilenames);
          assert(inFilenames.size() == 2);
        }
    }
};

// Forward declare some methods
moab::ErrorCode LoadTempestMesh(ToolContext& ctx, Mesh** tempest_mesh);
moab::ErrorCode ConvertTempestMeshToMOAB(ToolContext& ctx, Mesh* mesh, moab::Interface* mb);
moab::ErrorCode ConvertMOABMeshToTempest(ToolContext& ctx, moab::Interface* mb, Mesh* mesh);

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;

#ifdef MOAB_HAVE_MPI
  int proc_id = 0, size = 1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  ToolContext ctx;
  ctx.ParseCLOptions(argc, argv);

  Mesh* tempest_mesh=NULL;
  rval = LoadTempestMesh(ctx, &tempest_mesh);MB_CHK_ERR(rval);

  moab::Core* mbCore = new (std::nothrow) moab::Core;
  if (NULL == mbCore) return 1;

  // First convert the loaded Tempest mesh to MOAB
  rval = ConvertTempestMeshToMOAB(ctx, tempest_mesh, mbCore);MB_CHK_ERR(rval);
#if 0
  mbCore->print_database();
#endif
  mbCore->write_file("test.vtk");
  delete tempest_mesh;

  // Now let us re-convert the MOAB mesh back to Tempest representation
  Mesh* tempest_mesh_copy = new Mesh();
  rval = ConvertMOABMeshToTempest(ctx, mbCore, tempest_mesh_copy);MB_CHK_ERR(rval);
  tempest_mesh_copy->Write("test.exo");
  delete tempest_mesh_copy;

  delete mbCore;
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(0);
}

moab::ErrorCode LoadTempestMesh(ToolContext& ctx, Mesh** tempest_mesh)
{
    std::cout << "\nLoading TempestRemap Mesh object ...\n";
    switch(ctx.meshType) {
      case OVERLAP:
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        *tempest_mesh = GenerateOverlapMesh(ctx.inFilenames[0], ctx.inFilenames[1], ctx.outFilename, "exact", false);
        break;
      case ICO:
        *tempest_mesh = GenerateICOMesh(ctx.blockSize, ctx.computeDual, ctx.outFilename);
        break;
      case RLL:
        *tempest_mesh = GenerateRLLMesh(ctx.blockSize*2, ctx.blockSize, 0.0, 360.0, -90.0, 90.0, false, ctx.outFilename);
        break;
      case CS:
      default:
        *tempest_mesh = GenerateCSMesh(ctx.blockSize, false, ctx.outFilename);
        break;
    }

    if (!*tempest_mesh) {
      std::cout << "Tempest Mesh is not a complete object; Quitting...";
      exit(-1);
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ConvertTempestMeshToMOAB(ToolContext& ctx, Mesh* mesh, moab::Interface* mb)
{
  std::cout << "\nConverting TempestRemap Mesh object to MOAB representation ...\n";
  const NodeVector& nodes = mesh->nodes;
  const FaceVector& faces = mesh->faces;

  moab::ReadUtilIface* iface;
  moab::ErrorCode rval = mb->query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  // Set the data for the vertices
  std::vector<double*> arrays;
  moab::EntityHandle startv;
  rval = iface->get_node_coords(3, nodes.size(), 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");
  for (unsigned iverts=0; iverts < nodes.size(); ++iverts) {
      const Node& node = nodes[iverts];
      arrays[0][iverts] = node.x;
      arrays[1][iverts] = node.y;
      arrays[2][iverts] = node.z;
  }

  // We will assume all elements are of the same type - for now;
  // need a better way to categorize without doing a full pass first
  const unsigned lnum_v_per_elem = faces[0].edges.size(); // Linear elements: nedges = nverts ?
  if (ctx.meshType != OVERLAP && lnum_v_per_elem <= 4) {
      const unsigned num_v_per_elem = lnum_v_per_elem;
      moab::EntityHandle starte; // Connectivity
      moab::EntityHandle* conn;
      rval = iface->get_element_connect(faces.size(), num_v_per_elem, moab::MBPOLYGON, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
      for (unsigned ifaces=0,offset=0; ifaces < faces.size(); ++ifaces) {
          const Face& face = faces[ifaces];
          conn[offset++] = startv+face.edges[0].node[1];
          for (unsigned iedges=1; iedges < face.edges.size(); ++iedges) {
              conn[offset++] = startv+face.edges[iedges].node[1];
          }
      }
  }
  else {
      const int NMAXPOLYEDGES = 15;
      std::vector<unsigned> nPolys(NMAXPOLYEDGES,0);
      std::vector<std::vector<int> > typeNSeqs(NMAXPOLYEDGES);
      for (unsigned ifaces=0; ifaces < faces.size(); ++ifaces) {
          const int iType = faces[ifaces].edges.size();
          nPolys[iType]++;
          typeNSeqs[iType].push_back(ifaces);
      }
      for (unsigned iType=0; iType < NMAXPOLYEDGES; ++iType) {
          if (!nPolys[iType]) continue; // Nothing to do

          std::cout << "Found " << nPolys[iType] << " polygonal elements with " << iType << " edges.\n";
          const unsigned num_v_per_elem = iType;
          moab::EntityHandle starte; // Connectivity
          moab::EntityHandle* conn;

          // Allocate the connectivity array, depending on the element type
          switch(num_v_per_elem) {
            case 3:
              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, moab::MBTRI, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
              break;
            case 4:
              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, moab::MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
              break;
            default:
              rval = iface->get_element_connect(nPolys[iType], num_v_per_elem, moab::MBPOLYGON, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
              break;
          }

          for (unsigned ifaces=0,offset=0; ifaces < typeNSeqs[iType].size(); ++ifaces) {
              const Face& face = faces[typeNSeqs[iType][ifaces]];
              conn[offset++] = startv+face.edges[0].node[1];
              for (unsigned iedges=1; iedges < face.edges.size(); ++iedges) {
                  conn[offset++] = startv+face.edges[iedges].node[1];
              }
          }
      }
  }

  return moab::MB_SUCCESS;
}


moab::ErrorCode ConvertMOABMeshToTempest(ToolContext& , moab::Interface* mb, Mesh* mesh)
{
  moab::ErrorCode rval;

  std::cout << "\nConverting MOAB Mesh object to TempestRemap Mesh representation ...\n";
  NodeVector& nodes = mesh->nodes;
  FaceVector& faces = mesh->faces;

  moab::Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_ERR(rval);
  nodes.resize(verts.size());

  // Set the data for the vertices
  int inode=0;
  std::vector<double> coordx(verts.size()), coordy(verts.size()), coordz(verts.size());
  rval = mb->get_coords(verts, &coordx[0], &coordy[0], &coordz[0]);MB_CHK_ERR(rval);
  for (moab::Range::iterator iverts=verts.begin(); iverts!=verts.end(); ++iverts) {
      Node& node = nodes[inode];
      node.x = coordx[inode];
      node.y = coordy[inode];
      node.z = coordz[inode];
      inode++;
  }
  coordx.clear();
  coordy.clear();
  coordz.clear();

  moab::Range elems;
  rval = mb->get_entities_by_dimension(0, 2, elems);MB_CHK_ERR(rval);
  faces.resize(elems.size());

  int iface=0;
  for (moab::Range::iterator ielems=elems.begin(); ielems!=elems.end(); ++ielems) {
    Face& face = faces[iface++];

    // compute the number of edges per faces
    std::vector< moab::EntityHandle > face_edges;
    rval = mb->get_adjacencies(&(*ielems), 1, 1, true, face_edges);MB_CHK_ERR(rval);
    face.edges.resize(face_edges.size());

    for (unsigned iedges=0; iedges < face_edges.size(); ++iedges) {
        Edge& edge = face.edges[iedges];

        // get the connectivity for each edge
        const moab::EntityHandle* connect;
        int nnodes;
        rval = mb->get_connectivity(face_edges[iedges], connect, nnodes);MB_CHK_ERR(rval);

        // we expect only linear edges (2 nodes/edge)
        assert(nnodes == 2);

        // assign the edge nodes
        edge.node[0] = connect[0];
        edge.node[1] = connect[1];
    }
  }

  return moab::MB_SUCCESS;
}
