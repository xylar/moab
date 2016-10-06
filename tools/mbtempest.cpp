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
#include "moab/CpuTimer.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

// Intersection includes
#include "moab/Intx2MeshOnSphere.hpp"
#include "moab/IntxUtils.hpp"

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

enum TempestMeshType { CS=0, RLL=1, ICO=2, OVERLAP=3, OVERLAP_V2=4 };

struct ToolContext {
    int blockSize;
    std::vector<std::string> inFilenames;
    std::vector<Mesh*> meshes;
    std::string outFilename;
    TempestMeshType meshType;
    bool computeDual;
    bool computeWeights;

    ToolContext() : blockSize(5), outFilename("output.exo"), meshType(CS), computeDual(false), computeWeights(false)
    {
      inFilenames.resize(2);
      timer = new moab::CpuTimer();
    }

    ~ToolContext()
    {
      inFilenames.clear();
      outFilename.clear();
      delete timer;
    }

    void clear()
    {
      for (unsigned i=0; i < meshes.size(); ++i) delete meshes[i];
      meshes.clear();
    }

    void timer_push(std::string operation)
    {
      timer_ops = timer->time_since_birth();
      opName = operation;
    }

    void timer_pop()
    {
      std::cout << "\n[LOG] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
      opName.clear();
    }

    void ParseCLOptions(int argc, char* argv[]) {
        ProgOptions opts;
        int imeshType=0;
        std::string expectedFName="output.exo";

        opts.addOpt<int>("res,r", "Resolution of the mesh (default=5)", &blockSize);
        opts.addOpt<int>("type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP=3, OVERLAP_V2=4])", &imeshType);
        opts.addOpt<std::string>("file,f", "Output mesh filename (default=output.exo)", &outFilename);
        opts.addOpt<void>("dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual);
        opts.addOpt<void>("weights,w", "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)", &computeDual);
        opts.addOpt<std::string>("load,l", "Input mesh filenames (a source and target mesh)", &expectedFName);

        opts.parseCommandLine(argc, argv);

        switch(imeshType) {
            case 1:
                meshType = RLL;
                break;
            case 2:
                meshType = ICO;
                break;
            case 3:
                meshType = OVERLAP;
                break;
            case 4:
                meshType = OVERLAP_V2;
                break;
            case 0:
            default:
                meshType = CS;
                break;
        }

        if (meshType == OVERLAP || meshType == OVERLAP_V2) {
          opts.getOptAllArgs("load,l", inFilenames);
          assert(inFilenames.size() == 2);
        }
        expectedFName.clear();
    }
  private:
    moab::CpuTimer *timer;
    double timer_ops;
    std::string opName;
};

// Forward declare some methods
moab::ErrorCode LoadTempestMesh(std::string inputFilename, Mesh** tempest_mesh, bool meshValidate=false, bool constructEdgeMap=false);
moab::ErrorCode CreateTempestMesh(ToolContext& ctx, Mesh** tempest_mesh);
moab::ErrorCode ConvertTempestMeshToMOAB(ToolContext& ctx, Mesh* mesh, moab::Interface* mb, moab::EntityHandle& mesh_set);
moab::ErrorCode ConvertMOABMeshToTempest(ToolContext& ctx, moab::Interface* mb, Mesh* mesh, moab::EntityHandle mesh_set);

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;
  NcError error(NcError::verbose_nonfatal);

#ifdef MOAB_HAVE_MPI
  int proc_id = 0, size = 1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  ToolContext ctx;
  ctx.ParseCLOptions(argc, argv);

  Mesh* tempest_mesh=NULL;
  ctx.timer_push("create Tempest mesh");
  rval = CreateTempestMesh(ctx, &tempest_mesh);MB_CHK_ERR(rval);
  ctx.timer_pop();

  moab::Core* mbCore = new (std::nothrow) moab::Core;
  if (NULL == mbCore) return 1;

  // First convert the loaded Tempest mesh to MOAB
#if 0
  {
      moab::EntityHandle mesh_set;
      ctx.timer_push("convert Tempest mesh to MOAB");
      rval = ConvertTempestMeshToMOAB(ctx, tempest_mesh, mbCore, mesh_set);MB_CHK_ERR(rval);
      ctx.timer_pop();
//      mbCore->print_database();
//      mbCore->write_file("test.vtk");

      // Now let us re-convert the MOAB mesh back to Tempest representation
      Mesh tempest_mesh_copy;
      ctx.timer_push("re-convert MOAB mesh back to Tempest mesh");
      rval = ConvertMOABMeshToTempest(ctx, mbCore, &tempest_mesh_copy, mesh_set);MB_CHK_ERR(rval);
      ctx.timer_pop();
//      tempest_mesh_copy.Write("test.exo");
  }
#endif

  if (ctx.meshType == OVERLAP_V2)
  { // Compute intersections with MOAB

    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    moab::EntityHandle red, blue, tempestintx, intxset;

    assert(ctx.meshes.size() == 3);

    // Load the meshes and validate
    rval = ConvertTempestMeshToMOAB(ctx, ctx.meshes[0], mbCore, red);MB_CHK_ERR(rval);
    rval = ConvertTempestMeshToMOAB(ctx, ctx.meshes[1], mbCore, blue);MB_CHK_ERR(rval);
    rval = ConvertTempestMeshToMOAB(ctx, ctx.meshes[2], mbCore, tempestintx);MB_CHK_ERR(rval);
    rval = mbCore->write_mesh("tempest_intersection.h5m",&tempestintx,1);MB_CHK_ERR(rval);

    // print verbosely about the problem setting
    {
        moab::Range rintxverts, rintxelems;
        rval = mbCore->get_entities_by_dimension(red, 0, rintxverts);MB_CHK_ERR(rval);
        rval = mbCore->get_entities_by_dimension(red, 2, rintxelems);MB_CHK_ERR(rval);
        std::cout << "The red set contains " << rintxverts.size() << " vertices and " << rintxelems.size() << " elements \n";

        moab::Range bintxverts, bintxelems;
        rval = mbCore->get_entities_by_dimension(blue, 0, bintxverts);MB_CHK_ERR(rval);
        rval = mbCore->get_entities_by_dimension(blue, 2, bintxelems);MB_CHK_ERR(rval);
        std::cout << "The blue set contains " << bintxverts.size() << " vertices and " << bintxelems.size() << " elements \n";
    }

    const double epsrel=1.e-12;
    const double radius=1.0;
    moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(mbCore);
    mbintx->SetErrorTolerance(epsrel);
    mbintx->SetRadius(radius);

    rval = mbCore->create_meshset(moab::MESHSET_SET, intxset);MB_CHK_ERR(rval);
    ctx.timer_push("compute intersections with MOAB");
    rval = mbintx->intersect_meshes(red, blue, intxset);MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
    ctx.timer_pop();

    delete mbintx;
    {
        moab::Range intxelems;
        rval = mbCore->get_entities_by_dimension(intxset, 2, intxelems);MB_CHK_ERR(rval);
        std::cout << "\nThe intersection set contains " << intxelems.size() << " elements \n";

        double initial_area = area_on_sphere_lHuiller(mbCore, red, radius);
        double area_method1 = area_on_sphere_lHuiller(mbCore, intxset, radius);
        double area_method2 = area_on_sphere(mbCore, intxset, radius);

        std::cout << "initial area: " << initial_area << "\n";
        std::cout<< " area with l'Huiller: " << area_method1 << " with Girard: " << area_method2<< "\n";
        std::cout << " relative difference areas " << fabs(area_method1-area_method2)/area_method1 << "\n";
        std::cout << " relative error " << fabs(area_method1-initial_area)/area_method1 << "\n";
    }

    // Write out our computed intersection file
    rval = mbCore->write_mesh("moab_intersection.h5m",&intxset,1);MB_CHK_ERR(rval);

    if (ctx.computeWeights) {
      // Now let us re-convert the MOAB mesh back to Tempest representation
      Mesh overlap_tempest_mesh_copy;
      rval = ConvertMOABMeshToTempest(ctx, mbCore, &overlap_tempest_mesh_copy, intxset);MB_CHK_ERR(rval);

      ctx.timer_push("compute weights with the Tempest meshes");
      // Call to generate an offline map with the tempest meshes
      // OfflineMap* weightMap = GenerateOfflineMapWithMeshes(  NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
      OfflineMap* weightMap = GenerateOfflineMapWithMeshes(  NULL, *ctx.meshes[0], *ctx.meshes[1], overlap_tempest_mesh_copy,
                                                            "", "",     // std::string strInputMeta, std::string strOutputMeta,
                                                            "fv", "fv", // std::string strInputType, std::string strOutputType,
                                                            4, 4       // int nPin=4, int nPout=4,
  //                                                           bool fBubble=false, int fMonotoneTypeID=0,
  //                                                           bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
  //                                                           std::string strVariables="", std::string strOutputMap="",
  //                                                           std::string strInputData="", std::string strOutputData="",
  //                                                           std::string strNColName="", bool fOutputDouble=false,
  //                                                           std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0
                                                            );
      ctx.timer_pop();
      weightMap->Write("outWeights.nc");
      delete weightMap;
    }
  }

  // Clean up
  ctx.clear();
  delete mbCore;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(0);
}

moab::ErrorCode LoadTempestMesh(std::string inputFilename, Mesh** tempest_mesh, bool meshValidate, bool constructEdgeMap)
{
    std::cout << "\nLoading TempestRemap Mesh object from file = " << inputFilename << " ...\n";
    if (tempest_mesh) {
        NcError error(NcError::silent_nonfatal);

        Mesh* mesh;

        try {
            // Load input mesh
            std::cout << "Loading mesh ...\n";
            mesh = new Mesh(inputFilename);
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

        } catch(Exception & e) {
            std::cout << "TempestRemap ERROR: " << e.ToString() << "\n";
            return moab::MB_FAILURE;

        } catch(...) {
            return moab::MB_FAILURE;
        }

        *tempest_mesh = mesh;
    }
    return moab::MB_SUCCESS;
}

moab::ErrorCode CreateTempestMesh(ToolContext& ctx, Mesh** tempest_mesh)
{
    moab::ErrorCode rval = moab::MB_SUCCESS;
    std::cout << "\nCreating TempestRemap Mesh object ...\n";
    switch(ctx.meshType) {
      case OVERLAP:
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        *tempest_mesh = GenerateOverlapMesh(ctx.inFilenames[0], ctx.inFilenames[1], ctx.outFilename, "exact", true);
        break;
      case OVERLAP_V2:
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        Mesh *src, *dest;
        // Load the meshes and validate
        rval = LoadTempestMesh(ctx.inFilenames[0], &src, true, true);MB_CHK_ERR(rval);
        rval = LoadTempestMesh(ctx.inFilenames[1], &dest, true, true);MB_CHK_ERR(rval);
        ctx.meshes.push_back(src);
        ctx.meshes.push_back(dest);
        // Now let us construct the overlap mesh
        *tempest_mesh = GenerateOverlapWithMeshes(*src, *dest, "" /*ctx.outFilename*/, "exact", false);
        //*tempest_mesh = src;
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
    ctx.meshes.push_back(*tempest_mesh);

    if (!*tempest_mesh) {
      std::cout << "Tempest Mesh is not a complete object; Quitting...";
      exit(-1);
    }

    return rval;
}

moab::ErrorCode ConvertTempestMeshToMOAB(ToolContext& ctx, Mesh* mesh, moab::Interface* mb, moab::EntityHandle& mesh_set)
{
  moab::ErrorCode rval;

  std::cout << "\nConverting TempestRemap Mesh object to MOAB representation ...\n";
  const NodeVector& nodes = mesh->nodes;
  const FaceVector& faces = mesh->faces;

  rval = mb->create_meshset(moab::MESHSET_SET, mesh_set);MB_CHK_SET_ERR(rval, "Can't create new set");

  moab::ReadUtilIface* iface;
  rval = mb->query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

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
  moab::Range mbverts(startv, startv + nodes.size() - 1);
  mb->add_entities(mesh_set, mbverts);

  // We will assume all elements are of the same type - for now;
  // need a better way to categorize without doing a full pass first
  const unsigned lnum_v_per_elem = faces[0].edges.size(); // Linear elements: nedges = nverts ?
  if ((ctx.meshType != OVERLAP && ctx.meshType != OVERLAP_V2) && lnum_v_per_elem <= 4) {
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
          moab::Range mbcells(starte, starte + nPolys[iType] - 1);
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
          moab::Range edges;
          rval = mb->get_adjacencies(mbcells, 1, true, edges,
                                     moab::Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get edges");
      }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode ConvertMOABMeshToTempest(ToolContext& ctx, moab::Interface* mb, Mesh* mesh, moab::EntityHandle mesh_set)
{
  moab::ErrorCode rval;

  std::cout << "\nConverting MOAB Mesh object to TempestRemap Mesh representation ...\n";
  NodeVector& nodes = mesh->nodes;
  FaceVector& faces = mesh->faces;

  moab::Range verts;
  rval = mb->get_entities_by_dimension(mesh_set, 0, verts);MB_CHK_ERR(rval);
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
  rval = mb->get_entities_by_dimension(mesh_set, 2, elems);MB_CHK_ERR(rval);
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

  if (ctx.meshType == OVERLAP_V2) {
      // We need to preserve the source and target FaceID maps so that the offline map can be generated cleanly
      // vecSourceFaceIx, vecTargetFaceIx
  }

  mesh->RemoveZeroEdges();
  return moab::MB_SUCCESS;
}
