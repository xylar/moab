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
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOfflineMap.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "DebugOutput.hpp"

#ifndef MOAB_HAVE_MPI
    #error mbtempest tool requires MPI configuration
#endif

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"


struct ToolContext
{
        int blockSize;
        std::vector<std::string> inFilenames;
        std::vector<Mesh*> meshes;
        std::vector<moab::EntityHandle> meshsets;
        std::vector<int> disc_orders;
        std::vector<std::string> disc_methods;
        std::string outFilename;
        moab::TempestRemapper::TempestMeshType meshType;
        const int proc_id, n_procs;
        bool computeDual;
        bool computeWeights;
        bool fNoConservation;
        bool fVolumetric;
        moab::DebugOutput outStream;

        ToolContext ( int procid, int nprocs ) :
            blockSize ( 5 ), outFilename ( "output.exo" ), meshType ( moab::TempestRemapper::DEFAULT ),
            proc_id ( procid ), n_procs ( nprocs ),
            computeDual ( false ), computeWeights ( false ), fNoConservation ( false ), fVolumetric ( false ),
            outStream ( std::cout, procid )
        {
            inFilenames.resize ( 2 );
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
            meshes.clear();
        }

        void timer_push ( std::string operation )
        {
            timer_ops = timer->time_since_birth();
            opName = operation;
        }

        void timer_pop()
        {
            outStream.printf ( 0, "[LOG] Time taken to %s = %f\n", opName.c_str(), timer->time_since_birth() - timer_ops );
            // std::cout << "\n[LOG" << proc_id << "] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
            opName.clear();
        }

        void ParseCLOptions ( int argc, char* argv[] )
        {
            ProgOptions opts;
            int imeshType = 0;
            std::string expectedFName = "output.exo";
            std::string expectedMethod = "fv";
            int expectedOrder = 1;

            opts.addOpt<int> ( "res,r", "Resolution of the mesh (default=5)", &blockSize );
            opts.addOpt<int> ( "type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP_FILES=3, OVERLAP_MEMORY=4, OVERLAP_MOAB=5])", &imeshType );
            opts.addOpt<std::string> ( "file,f", "Output mesh filename (default=output.exo)", &outFilename );
            opts.addOpt<void> ( "dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual );
            opts.addOpt<void> ( "weights,w", "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)", &computeWeights );
            opts.addOpt<void> ( "noconserve,c", "Do not apply conservation to the resultant weights (relevant only when computing weights)", &fNoConservation );
            opts.addOpt<void> ( "volumetric,v", "Apply a volumetric projection to compute the weights (relevant only when computing weights)", &fVolumetric );
            opts.addOpt<std::string> ( "load,l", "Input mesh filenames (a source and target mesh)", &expectedFName );
            opts.addOpt<int> ( "order,o", "Discretization orders for the source and target solution fields", &expectedOrder );
            opts.addOpt<std::string> ( "method,m", "Discretization method for the source and target solution fields", &expectedMethod );

            opts.parseCommandLine ( argc, argv );

            switch ( imeshType )
            {
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
                    meshType = moab::TempestRemapper::OVERLAP_FILES;
                    break;

                case 4:
                    meshType = moab::TempestRemapper::OVERLAP_MEMORY;
                    break;

                case 5:
                    meshType = moab::TempestRemapper::OVERLAP_MOAB;
                    break;

                default:
                    meshType = moab::TempestRemapper::DEFAULT;
                    break;
            }

            if ( meshType > moab::TempestRemapper::ICO )
            {
                opts.getOptAllArgs ( "load,l", inFilenames );
                opts.getOptAllArgs ( "order,o", disc_orders );
                opts.getOptAllArgs ( "method,m", disc_methods );

                if ( disc_orders.size() == 0 )
                { disc_orders.resize ( 2, 1 ); }

                if ( disc_orders.size() == 1 )
                { disc_orders.push_back ( 1 ); }

                if ( disc_methods.size() == 0 )
                { disc_methods.resize ( 2, "fv" ); }

                if ( disc_methods.size() == 1 )
                { disc_methods.push_back ( "fv" ); }

                assert ( inFilenames.size() == 2 );
                assert ( disc_orders.size() == 2 );
                assert ( disc_methods.size() == 2 );
            }

            expectedFName.clear();
        }
    private:
        moab::CpuTimer* timer;
        double timer_ops;
        std::string opName;
};

// Forward declare some methods
moab::ErrorCode CreateTempestMesh ( ToolContext&, moab::TempestRemapper& remapper, Mesh* );

int main ( int argc, char* argv[] )
{
    moab::ErrorCode rval;
    NcError error ( NcError::verbose_nonfatal );
    std::stringstream sstr;

    int proc_id = 0, nprocs = 1;
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &proc_id );
    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );

    ToolContext ctx ( proc_id, nprocs );
    ctx.ParseCLOptions ( argc, argv );

    moab::Interface* mbCore = new ( std::nothrow ) moab::Core;

    if ( NULL == mbCore ) { return 1; }

    moab::ParallelComm* pcomm = new moab::ParallelComm ( mbCore, MPI_COMM_WORLD, 0 );

    moab::TempestRemapper remapper ( mbCore, pcomm );
    remapper.meshValidate = true;
    remapper.constructEdgeMap = true;
    remapper.initialize();

    Mesh* tempest_mesh = new Mesh();
    ctx.timer_push ( "create Tempest mesh" );
    rval = CreateTempestMesh ( ctx, remapper, tempest_mesh ); MB_CHK_ERR ( rval );
    ctx.timer_pop();

    // Some constant parameters
    const double epsrel = 1.e-8;
    const double radius_src = 1.0 /*2.0*acos(-1.0)*/;
    const double radius_dest = 1.0 /*2.0*acos(-1.0)*/;
    const double boxeps = 0.1;

    // Rescale the radius of both to compute the intersection
    if (ctx.meshsets.size() > 0) {
        rval = ScaleToRadius(mbCore, ctx.meshsets[0], radius_src);MB_CHK_ERR ( rval );
    }
    if (ctx.meshsets.size() > 1) {
        rval = ScaleToRadius(mbCore, ctx.meshsets[1], radius_dest);MB_CHK_ERR ( rval );
    }

    if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Compute intersections with MOAB
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert ( ctx.meshes.size() == 3 );

        rval = pcomm->check_all_shared_handles(); MB_CHK_ERR ( rval );

        // Load the meshes and validate
        rval = remapper.ConvertTempestMesh ( moab::Remapper::SourceMesh ); MB_CHK_ERR ( rval );
        rval = remapper.ConvertTempestMesh ( moab::Remapper::TargetMesh ); MB_CHK_ERR ( rval );
        remapper.SetMeshType ( moab::Remapper::IntersectedMesh, moab::TempestRemapper::OVERLAP_FILES );
        rval = remapper.ConvertTempestMesh ( moab::Remapper::IntersectedMesh ); MB_CHK_ERR ( rval );
        rval = mbCore->write_mesh ( "tempest_intersection.h5m", &ctx.meshsets[2], 1 ); MB_CHK_ERR ( rval );

        // print verbosely about the problem setting
        {
            moab::Range rintxverts, rintxelems;
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[0], 0, rintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[0], 2, rintxelems ); MB_CHK_ERR ( rval );
            ctx.outStream.printf ( 0, "The red set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size() );

            moab::Range bintxverts, bintxelems;
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[1], 0, bintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[1], 2, bintxelems ); MB_CHK_ERR ( rval );
            ctx.outStream.printf ( 0, "The blue set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size() );
        }

        moab::EntityHandle intxset; // == remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

        // Compute intersections with MOAB
        {
            // Create the intersection on the sphere object
            ctx.timer_push ( "setup the intersector" );

            moab::Intx2MeshOnSphere* mbintx = new moab::Intx2MeshOnSphere ( mbCore );
            mbintx->set_error_tolerance ( epsrel );
            mbintx->set_box_error ( boxeps );
            mbintx->set_radius_source_mesh ( radius_src );
            mbintx->set_radius_destination_mesh ( radius_dest );
            mbintx->set_parallel_comm ( pcomm );

            rval = mbintx->FindMaxEdges ( ctx.meshsets[0], ctx.meshsets[1] ); MB_CHK_ERR ( rval );

            moab::Range local_verts;
            rval = mbintx->build_processor_euler_boxes ( ctx.meshsets[1], local_verts ); MB_CHK_ERR ( rval );

            ctx.timer_pop();

            moab::EntityHandle covering_set;
            ctx.timer_push ( "communicate the mesh" );
            rval = mbintx->construct_covering_set ( ctx.meshsets[0], covering_set ); MB_CHK_ERR ( rval ); // lots of communication if mesh is distributed very differently
            ctx.timer_pop();

            // Now let's invoke the MOAB intersection algorithm in parallel with a
            // source and target mesh set representing two different decompositions
            ctx.timer_push ( "compute intersections with MOAB" );
            rval = mbCore->create_meshset ( moab::MESHSET_SET, intxset ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            rval = mbintx->intersect_meshes ( covering_set, ctx.meshsets[1], intxset ); MB_CHK_SET_ERR ( rval, "Can't compute the intersection of meshes on the sphere" );
            ctx.timer_pop();

            // free the memory
            delete mbintx;
        }

        {
            moab::Range intxelems, intxverts;
            rval = mbCore->get_entities_by_dimension ( intxset, 2, intxelems ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( intxset, 0, intxverts, true ); MB_CHK_ERR ( rval );
            ctx.outStream.printf ( 0, "The intersection set contains %lu elements and %lu vertices \n", intxelems.size(), intxverts.size() );

            double initial_area = area_on_sphere_lHuiller ( mbCore, ctx.meshsets[0], radius_src );
            double area_method1 = area_on_sphere_lHuiller ( mbCore, intxset, radius_src );
            double area_method2 = area_on_sphere ( mbCore, intxset, radius_src );

            ctx.outStream.printf ( 0, "initial area: %12.10f\n", initial_area );
            ctx.outStream.printf ( 0, " area with l'Huiller: %12.10f with Girard: %12.10f\n", area_method1, area_method2 );
            ctx.outStream.printf ( 0, " relative difference areas = %12.10e\n", fabs ( area_method1 - area_method2 ) / area_method1 );
            ctx.outStream.printf ( 0, " relative error = %12.10e\n", fabs ( area_method1 - initial_area ) / area_method1 );
        }

        // Write out our computed intersection file
        rval = mbCore->write_mesh ( "moab_intersection.h5m", &intxset, 1 ); MB_CHK_ERR ( rval );

        if ( ctx.computeWeights )
        {
            ctx.timer_push ( "compute weights with the Tempest meshes" );
            // Call to generate an offline map with the tempest meshes
            OfflineMap weightMap;
            int err = GenerateOfflineMapWithMeshes (  weightMap, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
                      "", "",     // std::string strInputMeta, std::string strOutputMeta,
                      ctx.disc_methods[0], ctx.disc_methods[1], // std::string strInputType, std::string strOutputType,
                      ctx.disc_orders[0], ctx.disc_orders[1]  // int nPin=4, int nPout=4,
                                                   );
            ctx.timer_pop();

            if ( err ) { rval = moab::MB_FAILURE; }
            else { weightMap.Write ( "outWeights.nc" ); }
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m

        rval = pcomm->check_all_shared_handles(); MB_CHK_ERR ( rval );

        // print verbosely about the problem setting
        {
            moab::Range rintxverts, rintxelems;
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[0], 0, rintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[0], 2, rintxelems ); MB_CHK_ERR ( rval );
            rval = fix_degenerate_quads ( mbCore, ctx.meshsets[0] ); MB_CHK_ERR ( rval );
            rval = positive_orientation ( mbCore, ctx.meshsets[0], radius_src ); MB_CHK_ERR ( rval );
            ctx.outStream.printf ( 0, "The red set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size() );

            moab::Range bintxverts, bintxelems;
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[1], 0, bintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( ctx.meshsets[1], 2, bintxelems ); MB_CHK_ERR ( rval );
            rval = fix_degenerate_quads ( mbCore, ctx.meshsets[1] ); MB_CHK_ERR ( rval );
            rval = positive_orientation ( mbCore, ctx.meshsets[1], radius_dest ); MB_CHK_ERR ( rval );
            ctx.outStream.printf ( 0, "The blue set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size() );
        }

        // Compute intersections with MOAB
        ctx.timer_push ( "setup and compute mesh intersections" );
        rval = remapper.ComputeOverlapMesh ( epsrel ); MB_CHK_ERR ( rval );
        ctx.timer_pop();

        // Write out our computed intersection file
        if ( true )
        {
            ctx.outStream.printf ( 0, "writing out the intersection mesh file to %s\n", "moab_intersection.h5m" );
            // rval = mbCore->add_entities ( ctx.meshsets[2], &ctx.meshsets[0], 2 ); MB_CHK_ERR ( rval );
            rval = mbCore->write_file ( "moab_intersection.h5m", NULL, "PARALLEL=WRITE_PART", &ctx.meshsets[2], 1 ); MB_CHK_ERR ( rval );
        }

        {
            double local_areas[3], global_areas[3]; // Array for Initial area, and through Method 1 and Method 2
            local_areas[0] = area_on_sphere_lHuiller ( mbCore, ctx.meshsets[1], radius_src );
            local_areas[1] = area_on_sphere_lHuiller ( mbCore, ctx.meshsets[2], radius_src );
            local_areas[2] = area_on_sphere ( mbCore, ctx.meshsets[2], radius_src );

            MPI_Allreduce ( &local_areas, &global_areas, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            if ( !proc_id )
            {
                ctx.outStream.printf ( 0, "initial area: %12.10f\n", global_areas[0] );
                ctx.outStream.printf ( 0, " area with l'Huiller: %12.10f with Girard: %12.10f\n", global_areas[1], global_areas[2] );
                ctx.outStream.printf ( 0, " relative difference areas = %12.10e\n", fabs ( global_areas[1] - global_areas[2] ) / global_areas[1] );
                ctx.outStream.printf ( 0, " relative error = %12.10e\n", fabs ( global_areas[1] - global_areas[0] ) / global_areas[1] );
            }
        }

        if ( ctx.computeWeights )
        {

            ctx.timer_push ( "in-memory transform overlap mesh (MOAB->Tempest)" );
            {
                // Now let us re-convert the MOAB mesh back to Tempest representation
                // rval = remapper.ConvertMeshToTempest(moab::Remapper::IntersectedMesh);MB_CHK_ERR(rval);
                rval = remapper.AssociateSrcTargetInOverlap(); MB_CHK_ERR ( rval );
                rval = remapper.ConvertMOABMesh_WithSortedEntitiesBySource(); MB_CHK_ERR ( rval );

                ctx.meshes[2] = remapper.GetMesh ( moab::Remapper::IntersectedMesh );
            }
            ctx.timer_pop();

            ctx.timer_push ( "setup computation of weights" );
            // Call to generate an offline map with the tempest meshes
            moab::TempestOfflineMap* weightMap = new moab::TempestOfflineMap ( &remapper );
            ctx.timer_pop();

            ctx.timer_push ( "compute weights with TempestRemap" );
            rval = weightMap->GenerateOfflineMap ( ctx.disc_methods[0], ctx.disc_methods[1],        // std::string strInputType, std::string strOutputType,
                                                   ctx.disc_orders[0],  ctx.disc_orders[1],  // int nPin=4, int nPout=4,
                                                   false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
                                                   ctx.fVolumetric, ctx.fNoConservation, false, // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                                   "", //"",   // std::string strVariables="", std::string strOutputMap="",
                                                   "", "",   // std::string strInputData="", std::string strOutputData="",
                                                   "", false,  // std::string strNColName="", bool fOutputDouble=false,
                                                   "", false, 0.0,   // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0,
                                                   false, false   // bool fInputConcave = false, bool fOutputConcave = false
                                                 );
            ctx.timer_pop();

#if 0
            // OfflineMap* weightMap;
            // weightMap = GenerateOfflineMapWithMeshes( NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
            //                         "", "",              // std::string strInputMeta, std::string strOutputMeta,
            //                         ctx.disc_methods[0], ctx.disc_methods[1],// std::string strInputType, std::string strOutputType,
            //                         ctx.disc_orders[0], ctx.disc_orders[1],  // int nPin=4, int nPout=4,
            //                         false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
            //                         false, false, (nprocs>1) // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
            //                         // std::string strVariables="", std::string strOutputMap="",
            //                         // std::string strInputData="", std::string strOutputData="",
            //                         // std::string strNColName="", bool fOutputDouble=false,
            //                         // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0
            //                                                     );

            // weightMap->m_vecSourceDimSizes.resize(ctx.meshes[0]->faces.size());
            // weightMap->m_vecTargetDimSizes.resize(ctx.meshes[1]->faces.size());

            // sstr.str("");
            // sstr << "outWeights_" << proc_id << ".nc";
            // weightMap->Write(sstr.str().c_str());
#endif

            delete weightMap;
        }
    }

    // Clean up
    ctx.clear();
    delete mbCore;

    MPI_Finalize();
    exit ( 0 );
}


moab::ErrorCode CreateTempestMesh ( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    moab::ErrorCode rval = moab::MB_SUCCESS;
    int err;

    if ( !ctx.proc_id ) { ctx.outStream.printf ( 0, "Creating TempestRemap Mesh object ...\n" ); }

    if ( ctx.meshType == moab::TempestRemapper::OVERLAP_FILES )
    {
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapMesh ( ctx.inFilenames[0], ctx.inFilenames[1], *tempest_mesh, ctx.outFilename, "exact", true );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Load the meshes and validate
        ctx.meshsets.resize ( 3 );
        ctx.meshes.resize ( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet ( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet ( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet ( moab::Remapper::IntersectedMesh );

        // First the source
        rval = remapper.LoadMesh ( moab::Remapper::SourceMesh, ctx.inFilenames[0], moab::TempestRemapper::DEFAULT ); MB_CHK_ERR ( rval );
        ctx.meshes[0] = remapper.GetMesh ( moab::Remapper::SourceMesh );

        // Next the target
        rval = remapper.LoadMesh ( moab::Remapper::TargetMesh, ctx.inFilenames[1], moab::TempestRemapper::DEFAULT ); MB_CHK_ERR ( rval );
        ctx.meshes[1] = remapper.GetMesh ( moab::Remapper::TargetMesh );

        // Now let us construct the overlap mesh, by calling TempestRemap interface directly
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapWithMeshes ( *ctx.meshes[0], *ctx.meshes[1], *tempest_mesh, "" /*ctx.outFilename*/, "exact", false );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            remapper.SetMesh ( moab::Remapper::IntersectedMesh, tempest_mesh );
            ctx.meshes[2] = remapper.GetMesh ( moab::Remapper::IntersectedMesh );
            // ctx.meshes.push_back(*tempest_mesh);
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        ctx.meshsets.resize ( 3 );
        ctx.meshes.resize ( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet ( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet ( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet ( moab::Remapper::IntersectedMesh );

        // Load the source mesh and validate
        rval = remapper.LoadNativeMesh ( ctx.inFilenames[0], ctx.meshsets[0], 0 ); MB_CHK_ERR ( rval );
        rval = remapper.ConvertMeshToTempest ( moab::Remapper::SourceMesh ); MB_CHK_ERR ( rval );
        ctx.meshes[0] = remapper.GetMesh ( moab::Remapper::SourceMesh );

        // Load the target mesh and validate
        rval = remapper.LoadNativeMesh ( ctx.inFilenames[1], ctx.meshsets[1], 0 ); MB_CHK_ERR ( rval );
        rval = remapper.ConvertMeshToTempest ( moab::Remapper::TargetMesh ); MB_CHK_ERR ( rval );
        ctx.meshes[1] = remapper.GetMesh ( moab::Remapper::TargetMesh );

        // Set the references for the overlap mesh
        // ctx.meshes[2] = remapper.GetMesh(moab::Remapper::IntersectedMesh);
    }
    else if ( ctx.meshType == moab::TempestRemapper::ICO )
    {
        err = GenerateICOMesh ( *tempest_mesh, ctx.blockSize, ctx.computeDual, ctx.outFilename );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::RLL )
    {
        err = GenerateRLLMesh ( *tempest_mesh,                  // Mesh& meshOut,
                                ctx.blockSize * 2, ctx.blockSize, // int nLongitudes, int nLatitudes,
                                0.0, 360.0,                       // double dLonBegin, double dLonEnd,
                                -90.0, 90.0,                      // double dLatBegin, double dLatEnd,
                                false, false,                     // bool fFlipLatLon, bool fForceGlobal,
                                "" /*ctx.inFilename*/, ctx.outFilename,  // std::string strInputFile, std::string strOutputFile,
                                true                              // bool fVerbose
                              );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else   // default
    {
        err = GenerateCSMesh ( *tempest_mesh, ctx.blockSize, false, ctx.outFilename );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }

    if ( ctx.meshType != moab::TempestRemapper::OVERLAP_MOAB && !tempest_mesh )
    {
        std::cout << "Tempest Mesh is not a complete object; Quitting...";
        exit ( -1 );
    }

    return rval;
}

