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
#include <sstream>
#include <cassert>

#include "moab/Core.hpp"
#include "moab/TempestRemapper.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

struct ToolContext {
  int blockSize;
  std::vector<std::string> inFilenames;
  std::vector<Mesh*> meshes;
  std::vector<moab::EntityHandle> meshsets;
  std::string outFilename;
  moab::TempestRemapper::TempestMeshType meshType;
  const int proc_id, n_procs;
  bool computeDual;
  bool computeWeights;

  ToolContext(const int procid, const int nprocs) : 
    blockSize(5), outFilename("output.exo"), meshType(moab::TempestRemapper::DEFAULT), 
    proc_id(procid), n_procs(nprocs),
    computeDual(false), computeWeights(false)
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
    for (unsigned i = 0; i < meshes.size(); ++i) delete meshes[i];
    meshes.clear();
  }

  void timer_push(std::string operation)
  {
    timer_ops = timer->time_since_birth();
    opName = operation;
  }

  void timer_pop()
  {
    if (!proc_id) std::cout << "\n[LOG] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
    // std::cout << "\n[LOG" << proc_id << "] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
    opName.clear();
  }

  void ParseCLOptions(int argc, char* argv[]) {
    ProgOptions opts;
    int imeshType = 0;
    std::string expectedFName = "output.exo";

    opts.addOpt<int>("res,r", "Resolution of the mesh (default=5)", &blockSize);
    opts.addOpt<int>("type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP=3, OVERLAP_V2=4, OVERLAP_MOAB=5])", &imeshType);
    opts.addOpt<std::string>("file,f", "Output mesh filename (default=output.exo)", &outFilename);
    opts.addOpt<void>("dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual);
    opts.addOpt<void>("weights,w", "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)", &computeWeights);
    opts.addOpt<std::string>("load,l", "Input mesh filenames (a source and target mesh)", &expectedFName);

    opts.parseCommandLine(argc, argv);

    switch (imeshType) {
    case 0:
      meshType = moab::TempestRemapper::CS;
      break;
    case 1:
      meshType = moab::TempestRemapper::RLL;
      break;
    case 2:
      meshType = moab::TempestRemapper::ICO;
      break;
    case 3:
      meshType = moab::TempestRemapper::OVERLAP;
      break;
    case 4:
      meshType = moab::TempestRemapper::OVERLAP_V2;
      break;
    case 5:
      meshType = moab::TempestRemapper::OVERLAP_MOAB;
      break;
    default:
      meshType = moab::TempestRemapper::DEFAULT;
      break;
    }

    if (meshType > moab::TempestRemapper::ICO) {
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
// moab::ErrorCode LoadTempestMesh(std::string inputFilename, Mesh** tempest_mesh, bool meshValidate=false, bool constructEdgeMap=false);
moab::ErrorCode CreateTempestMesh(ToolContext&, moab::TempestRemapper& remapper, Mesh** );

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;
  NcError error(NcError::verbose_nonfatal);
  std::stringstream sstr;

  int proc_id = 0, nprocs = 1;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
#endif

  ToolContext ctx(proc_id, nprocs);
  ctx.ParseCLOptions(argc, argv);

  moab::Interface* mbCore = new (std::nothrow) moab::Core;
  if (NULL == mbCore) return 1;

  moab::ParallelComm* pcomm = new moab::ParallelComm(mbCore, MPI_COMM_WORLD, 0);

  moab::TempestRemapper remapper (mbCore, pcomm);
  remapper.meshValidate = true;
  remapper.constructEdgeMap = true;
  remapper.initialize();

  Mesh* tempest_mesh = NULL;
  ctx.timer_push("create Tempest mesh");
  rval = CreateTempestMesh(ctx, remapper, &tempest_mesh); MB_CHK_ERR(rval);
  ctx.timer_pop();

  if (ctx.meshType == moab::TempestRemapper::OVERLAP_V2)
  {
    // Compute intersections with MOAB
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    // moab::EntityHandle red, blue, tempestintx;

    std::cout << "In overlap_v2, n(meshes) = " << ctx.meshes.size() << " and meshsets = " << ctx.meshsets.size() << "\n";

    assert(ctx.meshes.size() == 3);
 
    rval = pcomm->check_all_shared_handles();MB_CHK_ERR(rval);

    // Load the meshes and validate
    rval = remapper.ConvertTempestMesh(moab::Remapper::SourceMesh); MB_CHK_ERR(rval);
    rval = remapper.ConvertTempestMesh(moab::Remapper::TargetMesh); MB_CHK_ERR(rval);
    rval = remapper.ConvertTempestMesh(moab::Remapper::IntersectedMesh); MB_CHK_ERR(rval);
    rval = mbCore->write_mesh("tempest_intersection.h5m", &ctx.meshsets[2], 1); MB_CHK_ERR(rval);

    // print verbosely about the problem setting
    {
      moab::Range rintxverts, rintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 0, rintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 2, rintxelems); MB_CHK_ERR(rval);
      std::cout << "The red set contains " << rintxverts.size() << " vertices and " << rintxelems.size() << " elements \n";

      moab::Range bintxverts, bintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 0, bintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 2, bintxelems); MB_CHK_ERR(rval);
      std::cout << "The blue set contains " << bintxverts.size() << " vertices and " << bintxelems.size() << " elements \n";
    }

    const double epsrel = 1.e-8;
    const double radius = 1.0 /*2.0*acos(-1.0)*/;
    const double boxeps = 0.1;
    moab::EntityHandle intxset; // == remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // Compute intersections with MOAB
    {
      // Create the intersection on the sphere object
      ctx.timer_push("setup the intersector");
      
      moab::Range local_verts;
      moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(mbCore);
      mbintx->SetErrorTolerance(epsrel);
      mbintx->set_box_error(boxeps);
      mbintx->SetRadius(radius);

      rval = mbintx->FindMaxEdges(ctx.meshsets[0], ctx.meshsets[1]);MB_CHK_ERR(rval);

      rval = mbintx->build_processor_euler_boxes(ctx.meshsets[1], local_verts); MB_CHK_ERR(rval);
      
      ctx.timer_pop();
      
      moab::EntityHandle covering_set;
      ctx.timer_push("communicate the mesh");
      rval = mbintx->construct_covering_set(ctx.meshsets[0], covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
      ctx.timer_pop();

      // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
      ctx.timer_push("compute intersections with MOAB");
      rval = mbCore->create_meshset(moab::MESHSET_SET, intxset);MB_CHK_SET_ERR(rval, "Can't create new set");
      rval = mbintx->intersect_meshes(covering_set, ctx.meshsets[1], intxset); MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
      ctx.timer_pop();

      rval = fix_degenerate_quads(mbCore, ctx.meshsets[2]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[2], radius);MB_CHK_ERR(rval);

      // free the memory
      delete mbintx;
    }

    {
      moab::Range intxelems, intxverts;
      rval = mbCore->get_entities_by_dimension(intxset, 2, intxelems); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(intxset, 0, intxverts,true); MB_CHK_ERR(rval);
      if (!proc_id) std::cout << "\nThe intersection set contains " << intxelems.size() << " elements and " << intxverts.size() << " vertices\n";

      double initial_area = area_on_sphere_lHuiller(mbCore, ctx.meshsets[0], radius);
      double area_method1 = area_on_sphere_lHuiller(mbCore, intxset, radius);
      double area_method2 = area_on_sphere(mbCore, intxset, radius);

      std::cout << "initial area: " << initial_area << "\n";
      std::cout << " area with l'Huiller: " << area_method1 << " with Girard: " << area_method2 << "\n";
      std::cout << " relative difference areas " << fabs(area_method1 - area_method2) / area_method1 << "\n";
      std::cout << " relative error " << fabs(area_method1 - initial_area) / area_method1 << "\n";
    }

    // Write out our computed intersection file
    rval = mbCore->write_mesh("moab_intersection.h5m", &intxset, 1); MB_CHK_ERR(rval);

    if (ctx.computeWeights) {
      // Now let us re-convert the MOAB mesh back to Tempest representation
      // Mesh overlap_tempest_mesh_copy;
      // rval = remapper.ConvertMeshToTempest(&overlap_tempest_mesh_copy, intxset); MB_CHK_ERR(rval);

      ctx.timer_push("compute weights with the Tempest meshes");
      // Call to generate an offline map with the tempest meshes
      OfflineMap* weightMap = GenerateOfflineMapWithMeshes(  NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
      // OfflineMap* weightMap = GenerateOfflineMapWithMeshes(  NULL, *ctx.meshes[0], *ctx.meshes[1], overlap_tempest_mesh_copy,
                              "", "",     // std::string strInputMeta, std::string strOutputMeta,
                              "fv", "fv", // std::string strInputType, std::string strOutputType,
                              nprocs, nprocs       // int nPin=4, int nPout=4,
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
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB)
  {
    const double epsrel = 1.e-8;
    const double radius = 1.0 /*2.0*acos(-1.0)*/;
    const double boxeps = 0.1;
    // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m
    // moab::ParallelComm* pcomm = moab::ParallelComm::get_pcomm(mbCore, 0);

    rval = pcomm->check_all_shared_handles();MB_CHK_ERR(rval);

    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    // const unsigned rank = pcomm->proc_config().proc_rank();
    // const unsigned nprocs = pcomm->proc_config().proc_size();
    // print verbosely about the problem setting
    {
      moab::Range rintxverts, rintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 0, rintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 2, rintxelems); MB_CHK_ERR(rval);
      rval = fix_degenerate_quads(mbCore, ctx.meshsets[0]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[0], radius);MB_CHK_ERR(rval);
      if (!proc_id) std::cout << "The source set contains " << rintxverts.size() << " vertices and " << rintxelems.size() << " elements \n";

      moab::Range bintxverts, bintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 0, bintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 2, bintxelems); MB_CHK_ERR(rval);
      rval = fix_degenerate_quads(mbCore, ctx.meshsets[1]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[1], radius);MB_CHK_ERR(rval);
      if (!proc_id) std::cout << "The target set contains " << bintxverts.size() << " vertices and " << bintxelems.size() << " elements \n";
    }

    // Compute intersections with MOAB
    {
      // Create the intersection on the sphere object
      ctx.timer_push("setup the intersector");
      
      moab::Range local_verts;
      moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(mbCore);
      mbintx->SetErrorTolerance(epsrel);
      mbintx->set_box_error(boxeps);
      mbintx->SetRadius(radius);

      rval = mbintx->FindMaxEdges(ctx.meshsets[0], ctx.meshsets[1]);MB_CHK_ERR(rval);

      rval = mbintx->build_processor_euler_boxes(ctx.meshsets[1], local_verts); MB_CHK_ERR(rval);
      
      ctx.timer_pop();
      
      moab::EntityHandle covering_set;
      ctx.timer_push("communicate the mesh");
      rval = mbintx->construct_covering_set(ctx.meshsets[0], covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
      ctx.timer_pop();

      // // rval = mbCore->add_entities(ctx.meshsets[2], &covering_set, 1);MB_CHK_ERR(rval);
      // rval = mbCore->add_entities(ctx.meshsets[2], &ctx.meshsets[0], 2);MB_CHK_ERR(rval);
      // // rval = mbCore->add_entities(ctx.meshsets[2], local_verts);MB_CHK_ERR(rval);

      // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
      ctx.timer_push("compute intersections with MOAB");
      rval = mbintx->intersect_meshes(covering_set, ctx.meshsets[1], ctx.meshsets[2]); MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
      ctx.timer_pop();

      // moab::Range allVerts;
      // rval = mbCore->get_entities_by_dimension(ctx.meshsets[2], 0, allVerts, true); MB_CHK_ERR(rval);
      // std::cout << "[1] Total number of vertices in the sets (0, 1) = " << allVerts.size() << std::endl;
      // allVerts.clear();
      // rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 0, allVerts, true); MB_CHK_ERR(rval);
      // rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 0, allVerts, true); MB_CHK_ERR(rval);
      // std::cout << "[2] Total number of vertices in the sets (0, 1) = " << allVerts.size() << std::endl;
      // rval = mbCore->add_entities(ctx.meshsets[2], allVerts);MB_CHK_ERR(rval);
      rval = mbCore->add_entities(ctx.meshsets[2], &ctx.meshsets[0], 2);MB_CHK_ERR(rval);

      rval = fix_degenerate_quads(mbCore, ctx.meshsets[2]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[2], radius);MB_CHK_ERR(rval);

      // free the memory
      delete mbintx;
    }

    // std::cout << "MeshSets::: Source = " << ctx.meshsets[0] << " Target = " << ctx.meshsets[1] << " Overlap = " << ctx.meshsets[2] << std::endl;

    {
      // Now let us re-convert the MOAB mesh back to Tempest representation
      rval = remapper.ConvertMeshToTempest(moab::Remapper::IntersectedMesh);MB_CHK_ERR(rval);
      ctx.meshes[2] = remapper.GetMesh(moab::Remapper::IntersectedMesh);

      moab::Tag redPtag,bluePtag;
      rval = mbCore->tag_get_handle("RedParent", redPtag);MB_CHK_ERR(rval);
      rval = mbCore->tag_get_handle("BlueParent", bluePtag);MB_CHK_ERR(rval);

      moab::Range overlapEls, overlapVerts;
      int locsize[2], globsize[2];
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[2], 2, overlapEls); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[2], 0, overlapVerts, true); MB_CHK_ERR(rval);
      rval = pcomm->filter_pstatus(overlapVerts, PSTATUS_NOT_OWNED, PSTATUS_NOT);MB_CHK_ERR(rval);
      locsize[0] = overlapEls.size(); locsize[1] = overlapVerts.size();
      if (!proc_id) std::cout << "-- Local:  Intersection set contains " << locsize[0] << " elements and " << locsize[1] << " vertices\n";
      MPI_Reduce(&locsize, &globsize, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
      if (!proc_id) std::cout << "-- Global: Intersection set contains " << globsize[0] << " elements and " << globsize[1] << " vertices\n";

      // Overlap mesh: mesh[2]
      ctx.meshes[2]->vecSourceFaceIx.resize(overlapEls.size());
      ctx.meshes[2]->vecTargetFaceIx.resize(overlapEls.size());
      ctx.meshes[2]->ConstructEdgeMap();

      rval = mbCore->tag_get_data(redPtag,  overlapEls, &ctx.meshes[2]->vecSourceFaceIx[0]); MB_CHK_ERR(rval);
      rval = mbCore->tag_get_data(bluePtag, overlapEls, &ctx.meshes[2]->vecTargetFaceIx[0]); MB_CHK_ERR(rval);
    }

    // Write out our computed intersection file
    // rval = mbCore->write_mesh("moab_intersection.h5m", &ctx.meshsets[2], 1); MB_CHK_ERR(rval);
    rval = mbCore->write_file("moab_intersection.h5m", NULL, "PARALLEL=WRITE_PART", &ctx.meshsets[2], 1); MB_CHK_ERR(rval);

    {
      double local_areas[3], global_areas[3]; // Array for Initial area, and through Method 1 and Method 2
      local_areas[0] = area_on_sphere_lHuiller(mbCore, ctx.meshsets[1], radius);
      local_areas[1] = area_on_sphere_lHuiller(mbCore, ctx.meshsets[2], radius);
      local_areas[2] = area_on_sphere(mbCore, ctx.meshsets[2], radius);

      MPI_Allreduce(&local_areas, &global_areas, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if (!proc_id) {
        std::cout << "initial area: " << global_areas[0] << "\n";
        std::cout << " area with l'Huiller: " << global_areas[1] << " with Girard: " << global_areas[2] << "\n";
        std::cout << " relative difference areas " << fabs(global_areas[1] - global_areas[2]) / global_areas[1] << "\n";
        std::cout << " relative error " << fabs(global_areas[1] - global_areas[0]) / global_areas[1] << "\n\n";
      }
    }

    if (ctx.computeWeights) {

      // std::cout << "Source = " << ctx.meshes[0] << " Target = " << ctx.meshes[1] << " Overlap = " << ctx.meshes[2] << std::endl;

      ctx.timer_push("compute weights with the Tempest meshes");
      // Call to generate an offline map with the tempest meshes
      OfflineMap* weightMap = new OfflineMap();
      weightMap = GenerateOfflineMapWithMeshes( NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
                              "", "",              // std::string strInputMeta, std::string strOutputMeta,
                              "fv", "fv",          // std::string strInputType, std::string strOutputType,
                              nprocs, nprocs,  // int nPin=4, int nPout=4,
                              false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
                              false, false, (nprocs>1) // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                              //                                                           std::string strVariables="", std::string strOutputMap="",
                              //                                                           std::string strInputData="", std::string strOutputData="",
                              //                                                           std::string strNColName="", bool fOutputDouble=false,
                              //                                                           std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0
                                                          );

      // weightMap->m_vecSourceDimSizes.resize(ctx.meshes[0]->faces.size());
      // weightMap->m_vecTargetDimSizes.resize(ctx.meshes[1]->faces.size());

      rval = remapper.ExchangeGhostWeights(weightMap);MB_CHK_ERR(rval);
      ctx.timer_pop();
      sstr.str("");
      sstr << "outWeights_" << proc_id << ".nc";
      weightMap->Write(sstr.str().c_str());
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


moab::ErrorCode CreateTempestMesh(ToolContext& ctx, moab::TempestRemapper& remapper, Mesh** tempest_mesh)
{
  moab::ErrorCode rval = moab::MB_SUCCESS;

  std::cout << "\nCreating TempestRemap Mesh object ...\n";
  if (ctx.meshType == moab::TempestRemapper::OVERLAP) {
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    *tempest_mesh = GenerateOverlapMesh(ctx.inFilenames[0], ctx.inFilenames[1], ctx.outFilename, "exact", true);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_V2) {
    // Load the meshes and validate
    ctx.meshsets.resize(3);
    ctx.meshes.resize(3);
    ctx.meshsets[0] = remapper.GetMeshSet(moab::Remapper::SourceMesh);
    ctx.meshsets[1] = remapper.GetMeshSet(moab::Remapper::TargetMesh);
    ctx.meshsets[2] = remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // First the source
    rval = remapper.LoadMesh(moab::Remapper::SourceMesh, ctx.inFilenames[0], moab::TempestRemapper::DEFAULT); MB_CHK_ERR(rval);
    ctx.meshes[0] = remapper.GetMesh(moab::Remapper::SourceMesh);
    
    // Next the target
    rval = remapper.LoadMesh(moab::Remapper::TargetMesh, ctx.inFilenames[1], moab::TempestRemapper::DEFAULT); MB_CHK_ERR(rval);
    ctx.meshes[1] = remapper.GetMesh(moab::Remapper::TargetMesh);

    // Now let us construct the overlap mesh, by calling TempestRemap interface directly
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    *tempest_mesh = GenerateOverlapWithMeshes(*ctx.meshes[0], *ctx.meshes[1], "" /*ctx.outFilename*/, "exact", false);
    remapper.SetMesh(moab::Remapper::IntersectedMesh, *tempest_mesh);
    // ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB) {
    ctx.meshsets.resize(3);
    ctx.meshes.resize(3);
    ctx.meshsets[0] = remapper.GetMeshSet(moab::Remapper::SourceMesh);
    ctx.meshsets[1] = remapper.GetMeshSet(moab::Remapper::TargetMesh);
    ctx.meshsets[2] = remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // Load the source mesh and validate
    std::cout << "Loading source mesh\n";
    rval = remapper.LoadNativeMesh(ctx.inFilenames[0], ctx.meshsets[0], 0); MB_CHK_ERR(rval);
    std::cout << "Loading source MOAB mesh to Tempest\n";
    rval = remapper.ConvertMeshToTempest(moab::Remapper::SourceMesh); MB_CHK_ERR(rval);
    ctx.meshes[0] = remapper.GetMesh(moab::Remapper::SourceMesh);
    std::cout << "All done.\n";

    // Load the target mesh and validate 
    rval = remapper.LoadNativeMesh(ctx.inFilenames[1], ctx.meshsets[1], 0); MB_CHK_ERR(rval);
    rval = remapper.ConvertMeshToTempest(moab::Remapper::TargetMesh); MB_CHK_ERR(rval);
    ctx.meshes[1] = remapper.GetMesh(moab::Remapper::TargetMesh);

    // Set the references for the overlap mesh
    // ctx.meshes[2] = remapper.GetMesh(moab::Remapper::IntersectedMesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::ICO) {
    *tempest_mesh = GenerateICOMesh(ctx.blockSize, ctx.computeDual, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::RLL) {
    *tempest_mesh = GenerateRLLMesh(ctx.blockSize * 2, ctx.blockSize, 0.0, 360.0, -90.0, 90.0, false, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else { // default
    *tempest_mesh = GenerateCSMesh(ctx.blockSize, false, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }

  if (ctx.meshType != moab::TempestRemapper::OVERLAP_MOAB && !*tempest_mesh) {
    std::cout << "Tempest Mesh is not a complete object; Quitting...";
    exit(-1);
  }

  return rval;
}

