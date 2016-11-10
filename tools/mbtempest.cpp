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
    std::string outFilename;
    moab::TempestMeshType meshType;
    bool computeDual;
    bool computeWeights;

    ToolContext() : blockSize(5), outFilename("output.exo"), meshType(moab::CS), computeDual(false), computeWeights(false)
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
                meshType = moab::RLL;
                break;
            case 2:
                meshType = moab::ICO;
                break;
            case 3:
                meshType = moab::OVERLAP;
                break;
            case 4:
                meshType = moab::OVERLAP_V2;
                break;
            case 0:
            default:
                meshType = moab::CS;
                break;
        }

        if (meshType == moab::OVERLAP || meshType == moab::OVERLAP_V2) {
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
moab::ErrorCode CreateTempestMesh(ToolContext& ctx, Mesh** tempest_mesh);

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

  moab::Interface* mbCore = new (std::nothrow) moab::Core;
  if (NULL == mbCore) return 1;

  // First convert the loaded Tempest mesh to MOAB
#if 0
  {
      moab::EntityHandle mesh_set;
      ctx.timer_push("convert Tempest mesh to MOAB");
      rval = TempestRemapper::TempestRemapper::ConvertTempestMeshToMOAB(ctx, tempest_mesh, mbCore, mesh_set);MB_CHK_ERR(rval);
      ctx.timer_pop();TempestRemapper::
//      mbCore->print_databaseTempestRemapper::();
//      mbCore->write_file("test.vtk");

      // Now let us re-convert the MOAB mesh back to Tempest representation
      Mesh tempest_mesh_copy;
      ctx.timer_push("re-convert MOAB mesh back to Tempest mesh");
      rval = ConvertMOABMeshToTempest(ctx, mbCore, &tempest_mesh_copy, mesh_set);MB_CHK_ERR(rval);
      ctx.timer_pop();
//      tempest_mesh_copy.Write("test.exo");
  }
#endif

  if (ctx.meshType == moab::OVERLAP_V2)
  { // Compute intersections with MOAB

    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    moab::EntityHandle red, blue, tempestintx, intxset;

    assert(ctx.meshes.size() == 3);

    // Load the meshes and validate
    rval = moab::TempestRemapper::ConvertTempestMeshToMOAB(ctx.meshType, mbCore, ctx.meshes[0], red);MB_CHK_ERR(rval);
    rval = moab::TempestRemapper::ConvertTempestMeshToMOAB(ctx.meshType, mbCore, ctx.meshes[1], blue);MB_CHK_ERR(rval);
    rval = moab::TempestRemapper::ConvertTempestMeshToMOAB(ctx.meshType, mbCore, ctx.meshes[2], tempestintx);MB_CHK_ERR(rval);
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
      rval = moab::TempestRemapper::ConvertMOABMeshToTempest(mbCore, &overlap_tempest_mesh_copy, intxset);MB_CHK_ERR(rval);

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


moab::ErrorCode CreateTempestMesh(ToolContext& ctx, Mesh** tempest_mesh)
{
    moab::ErrorCode rval = moab::MB_SUCCESS;
    std::cout << "\nCreating TempestRemap Mesh object ...\n";
    switch(ctx.meshType) {
      case moab::OVERLAP:
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        *tempest_mesh = GenerateOverlapMesh(ctx.inFilenames[0], ctx.inFilenames[1], ctx.outFilename, "exact", true);
        break;
      case moab::OVERLAP_V2:
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        Mesh *src, *dest;
        // Load the meshes and validate
        rval = moab::TempestRemapper::LoadTempestMesh(ctx.inFilenames[0], &src, true, true);MB_CHK_ERR(rval);
        rval = moab::TempestRemapper::LoadTempestMesh(ctx.inFilenames[1], &dest, true, true);MB_CHK_ERR(rval);
        ctx.meshes.push_back(src);
        ctx.meshes.push_back(dest);
        // Now let us construct the overlap mesh
        *tempest_mesh = GenerateOverlapWithMeshes(*src, *dest, "" /*ctx.outFilename*/, "exact", false);
        //*tempest_mesh = src;
        break;
      case moab::ICO:
        *tempest_mesh = GenerateICOMesh(ctx.blockSize, ctx.computeDual, ctx.outFilename);
        break;
      case moab::RLL:
        *tempest_mesh = GenerateRLLMesh(ctx.blockSize*2, ctx.blockSize, 0.0, 360.0, -90.0, 90.0, false, ctx.outFilename);
        break;
      case moab::CS:
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

