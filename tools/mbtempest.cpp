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

/* Exit values */
#define SUCCESS 0
#define USAGE_ERROR 1
#define TEMPEST_ERROR 2
#define NOT_IMPLEMENTED 3

using namespace moab;

moab::ErrorCode translate_tempest_mesh(Mesh* mesh, moab::Interface* mb, int meshType);

enum MeshType { CS=0, RLL=1, ICO=2, OVERLAP=3 };

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;
  int blockSize = 5;
  std::string expectedFName="output.exo";
  std::vector<std::string> inFilenames;
  std::string outFilename="output.exo";
  int meshType=0;
  bool computeDual=false;

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  ProgOptions opts;

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

  std::vector<char*> progargs;
  Mesh* tempest_mesh;
  switch(meshType) {
    case OVERLAP:
      // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
      tempest_mesh = GenerateOverlapMesh(inFilenames[0], inFilenames[1], outFilename, "exact", false);
      break;
    case ICO:
      tempest_mesh = GenerateICOMesh(blockSize, computeDual, outFilename);
      break;
    case RLL:
      tempest_mesh = GenerateRLLMesh(blockSize*2, blockSize, 0.0, 360.0, -90.0, 90.0, false, outFilename);
      break;
    case CS:
    default:
      tempest_mesh = GenerateCSMesh(blockSize, false, outFilename);
      break;
  }

  if (!tempest_mesh) {
    std::cout << "Tempest Mesh is not a complete object; Quitting...";
    exit(TEMPEST_ERROR);
  }

  Core* mbCore = new (std::nothrow) Core;
  if (NULL == mbCore) {
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  ErrorCode rval = translate_tempest_mesh(tempest_mesh, mbCore, meshType);MB_CHK_ERR(rval);
#if 1
  mbCore->print_database();
  mbCore->write_file("test.vtk");
#endif

  delete mbCore;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(SUCCESS);
}

moab::ErrorCode translate_tempest_mesh(Mesh* mesh, moab::Interface* mb, int meshType)
{
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
  if (meshType != OVERLAP && lnum_v_per_elem <= 4) {
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
