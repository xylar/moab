/*This unit test is for the convergence study of high order reconstruction under uniform refinement*/
#include <iostream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/NestedRefine.hpp"
#include "moab/DiscreteGeometry/DGMSolver.hpp"
#include "moab/DiscreteGeometry/HiReconstruction.hpp"
#include "TestUtil.hpp"
#include "geomObject.cpp"
#include <math.h>
#include <stdlib.h>

#ifdef MOAB_HAVE_MPI
    #include "moab/ParallelComm.hpp"
    #include "MBParallelConventions.h"
    #include "ReadParallel.hpp"
    #include "moab/FileOptions.hpp"
    #include "MBTagConventions.hpp"
    #include "moab_mpi.h"
#endif

using namespace moab;

#define nsamples 10

#ifdef MOAB_HAVE_MPI
    std::string read_options;
#endif

ErrorCode test_closedsurface_mesh ( const char* filename, int* level_degrees, int num_levels, int degree, bool interp, int dim, geomObject* obj )
{
    Core moab;
    Interface* mbImpl = &moab;
    ParallelComm* pc = NULL;
    EntityHandle fileset;
    ErrorCode error;
    error = mbImpl->create_meshset ( moab::MESHSET_SET, fileset ); MB_CHK_ERR ( error );
    //load mesh from file
#ifdef MOAB_HAVE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    EntityHandle partnset;
    error  = mbImpl->create_meshset ( moab::MESHSET_SET, partnset ); MB_CHK_ERR ( error );
    pc = moab::ParallelComm::get_pcomm ( mbImpl, partnset, &comm );
    int procs = 1;
    MPI_Comm_size ( comm, &procs );

    if ( procs > 1 )
    {
        read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        error = mbImpl->load_file ( filename, &fileset, read_options.c_str() ); MB_CHK_ERR ( error );
    }
    else
    {
        error = mbImpl->load_file ( filename, &fileset ); MB_CHK_ERR ( error );
    }

#else
    error = mbImpl->load_file ( filename, &fileset ); MB_CHK_ERR ( error );
#endif

    //Generate hierarchy
    //error = refine_entities(&moab, pc, fileset, level_degrees, num_levels, true);  CHECK_ERR(error);
    NestedRefine uref ( &moab, pc, fileset );
    std::vector<EntityHandle> meshes;
    error = uref.generate_mesh_hierarchy ( num_levels, level_degrees, meshes ); MB_CHK_ERR ( error );

    //Perform high order reconstruction on level 0 mesh
    HiReconstruction hirec ( &moab, pc, fileset );
    assert ( dim == 2 );
    error = hirec.reconstruct3D_surf_geom ( degree, interp, false ); MB_CHK_ERR ( error );

    //High order projection
    Range verts, vorig;
    error = mbImpl->get_entities_by_dimension ( meshes.back(), 0, verts ); MB_CHK_ERR ( error );
    error = mbImpl->get_entities_by_dimension ( meshes.front(), 0, vorig ); MB_CHK_ERR ( error );
    int nvorig = vorig.size();
    double l1err = 0, l2err = 0, linferr = 0;

    for ( Range::iterator ivert = verts.begin() + nvorig; ivert != verts.end(); ++ivert )
    {
        //locate the element in level 0 mesh, on which *ivert is lying
        EntityHandle currvert = *ivert;

        std::vector<EntityHandle> parentEntities;
        error = uref.get_adjacencies(*ivert, 2, parentEntities ); MB_CHK_ERR ( error );
        assert ( parentEntities.size() );

        EntityHandle rootelem;
        error = uref.child_to_parent ( parentEntities[0], num_levels , 0, &rootelem ); MB_CHK_ERR ( error );

        //compute the natural coordinates of *ivert in this element
        assert ( TYPE_FROM_HANDLE ( rootelem ) == MBTRI );

        const EntityHandle* conn; int nvpe = 3;
        error = mbImpl->get_connectivity ( rootelem, conn, nvpe ); MB_CHK_ERR ( error );

        std::vector<double> cornercoords ( 3 * nvpe ), currcoords ( 3 );
        error = mbImpl->get_coords ( conn, nvpe, & ( cornercoords[0] ) ); MB_CHK_ERR ( error );
        error = mbImpl->get_coords ( &currvert, 1, & ( currcoords[0] ) ); MB_CHK_ERR ( error );

        std::vector<double> naturalcoords2fit ( 3 );
        DGMSolver::get_tri_natural_coords ( 3, & ( cornercoords[0] ), 1, & ( currcoords[0] ), & ( naturalcoords2fit[0] ) );

        //project onto the estimated geometry
        double hicoords[3];
        error = hirec.hiproj_walf_in_element ( rootelem, nvpe, 1, & ( naturalcoords2fit[0] ), hicoords ); MB_CHK_ERR ( error );

        //estimate error
        double err = obj->compute_projecterror ( 3, hicoords );
        l1err += err; l2err += err * err; linferr = std::max ( linferr, err );
    }

    l1err /= verts.size() - nvorig; l2err = sqrt ( l2err / ( verts.size() - nvorig ) );
    std::cout << "L1 error " << l1err << " L2 error " << l2err << " Linf error " << linferr << std::endl;
    return error;
}

ErrorCode closedsurface_uref_hirec_convergence_study ( const char* filename, int* level_degrees, int num_levels, std::vector<int>& degs2fit, bool interp, geomObject* obj )
{
    Core moab;
    Interface* mbImpl = &moab;
    ParallelComm* pc = NULL;
    EntityHandle fileset;
    ErrorCode error;
    error = mbImpl->create_meshset ( moab::MESHSET_SET, fileset ); MB_CHK_ERR ( error );
    //load mesh from file
#ifdef MOAB_HAVE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    EntityHandle partnset;
    error  = mbImpl->create_meshset ( moab::MESHSET_SET, partnset ); MB_CHK_ERR ( error );
    pc = moab::ParallelComm::get_pcomm ( mbImpl, partnset, &comm );
    int procs = 1;
    MPI_Comm_size ( comm, &procs );

    if ( procs > 1 )
    {
        read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        error = mbImpl->load_file ( filename, &fileset, read_options.c_str() ); MB_CHK_ERR ( error );
    }
    else
    {
        error = mbImpl->load_file ( filename, &fileset ); MB_CHK_ERR ( error );
    }

#else
    error = mbImpl->load_file ( filename, &fileset ); MB_CHK_ERR ( error );
#endif
    //Generate hierarchy
    NestedRefine uref ( &moab, pc, fileset );
    std::vector<EntityHandle> meshes;
    std::vector<int> meshsizes;
    error = uref.generate_mesh_hierarchy ( num_levels, level_degrees, meshes ); MB_CHK_ERR ( error );
    std::vector< std::vector<double> > geoml1errs ( 1 + degs2fit.size() ), geoml2errs ( 1 + degs2fit.size() ), geomlinferrs ( 1 + degs2fit.size() );

    //Perform high order reconstruction on each level of mesh (projected onto exact geometry) and estimate geometric error with various degrees of fitting
    for ( size_t i = 0; i < meshes.size(); ++i )
    {
        EntityHandle& mesh = meshes[i];
        //project onto exact geometry since each level with uref has only linear coordinates
        Range verts;
        error = mbImpl->get_entities_by_dimension ( mesh, 0, verts ); MB_CHK_ERR ( error );

        for ( Range::iterator ivert = verts.begin(); ivert != verts.end(); ++ivert )
        {
            EntityHandle currvert = *ivert;
            double currcoords[3], exactcoords[3];
            error = mbImpl->get_coords ( &currvert, 1, currcoords );
            obj->project_points2geom ( 3, currcoords, exactcoords, NULL );
            error = mbImpl->set_coords ( &currvert, 1, exactcoords );
        }

        //generate random points on each elements, assument 3D coordinates
        Range elems;
        error = mbImpl->get_entities_by_dimension ( mesh, 2, elems ); MB_CHK_ERR ( error );
        meshsizes.push_back ( elems.size() );
        int nvpe = TYPE_FROM_HANDLE ( *elems.begin() ) == MBTRI ? 3 : 4;
        std::vector<double> testpnts, testnaturalcoords;
        testpnts.reserve ( 3 * elems.size() *nsamples ); testnaturalcoords.reserve ( nvpe * elems.size() *nsamples );

        for ( Range::iterator ielem = elems.begin(); ielem != elems.end(); ++ielem )
        {
            EntityHandle currelem = *ielem;
            std::vector<EntityHandle> conn;
            error = mbImpl->get_connectivity ( &currelem, 1, conn ); MB_CHK_ERR ( error );
            nvpe = conn.size();
            std::vector<double> elemcoords ( 3 * conn.size() );
            error = mbImpl->get_coords ( & ( conn[0] ), conn.size(), & ( elemcoords[0] ) ); MB_CHK_ERR ( error );
            EntityType type = TYPE_FROM_HANDLE ( currelem );

            for ( int s = 0; s < nsamples; ++s )
            {
                if ( type == MBTRI )
                {
                    double a = ( double ) rand() / RAND_MAX, b = ( double ) rand() / RAND_MAX, c = ( double ) rand() / RAND_MAX, sum;
                    sum = a + b + c;

                    if ( sum < 1e-12 )
                    {
                        --s; continue;
                    }
                    else
                    {
                        a /= sum, b /= sum, c /= sum;
                    }

                    testpnts.push_back ( a * elemcoords[0] + b * elemcoords[3] + c * elemcoords[6] );
                    testpnts.push_back ( a * elemcoords[1] + b * elemcoords[4] + c * elemcoords[7] );
                    testpnts.push_back ( a * elemcoords[2] + b * elemcoords[5] + c * elemcoords[8] );
                    testnaturalcoords.push_back ( a );
                    testnaturalcoords.push_back ( b );
                    testnaturalcoords.push_back ( c );
                }
                else if ( type == MBQUAD )
                {
                    double xi = ( double ) rand() / RAND_MAX, eta = ( double ) rand() / RAND_MAX, N[4];
                    xi = 2 * xi - 1; eta = 2 * eta - 1;
                    N[0] = ( 1 - xi ) * ( 1 - eta ) / 4, N[1] = ( 1 + xi ) * ( 1 - eta ) / 4, N[2] = ( 1 + xi ) * ( 1 + eta ) / 4, N[3] = ( 1 - xi ) * ( 1 + eta ) / 4;
                    testpnts.push_back ( N[0]*elemcoords[0] + N[1]*elemcoords[3] + N[2]*elemcoords[6] + N[3]*elemcoords[9] );
                    testpnts.push_back ( N[0]*elemcoords[1] + N[1]*elemcoords[4] + N[2]*elemcoords[7] + N[3]*elemcoords[10] );
                    testpnts.push_back ( N[0]*elemcoords[2] + N[1]*elemcoords[5] + N[2]*elemcoords[8] + N[3]*elemcoords[11] );
                    testnaturalcoords.push_back ( N[0] );
                    testnaturalcoords.push_back ( N[1] );
                    testnaturalcoords.push_back ( N[2] );
                    testnaturalcoords.push_back ( N[3] );
                }
                else
                {
                    throw std::invalid_argument ( "NOT SUPPORTED ELEMENT TYPE" );
                }
            }
        }

        //Compute error of linear interpolation
        double l1err, l2err, linferr;
        obj->compute_projecterror ( 3, elems.size() *nsamples, & ( testpnts[0] ), l1err, l2err, linferr );
        geoml1errs[0].push_back ( l1err ); geoml2errs[0].push_back ( l2err ); geomlinferrs[0].push_back ( linferr );
        //Perform high order projection and compute error
        HiReconstruction hirec ( &moab, pc, mesh );

        for ( size_t ideg = 0; ideg < degs2fit.size(); ++ideg )
        {
            //High order reconstruction
            error = hirec.reconstruct3D_surf_geom ( degs2fit[ideg], interp, false, true ); MB_CHK_ERR ( error );
            int index = 0;

            for ( Range::iterator ielem = elems.begin(); ielem != elems.end(); ++ielem, ++index )
            {
                //Projection
                error = hirec.hiproj_walf_in_element ( *ielem, nvpe, nsamples, & ( testnaturalcoords[nvpe * nsamples * index] ), & ( testpnts[3 * nsamples * index] ) ); MB_CHK_ERR ( error );
            }

            //Estimate error
            obj->compute_projecterror ( 3, elems.size() *nsamples, & ( testpnts[0] ), l1err, l2err, linferr );
            geoml1errs[ideg + 1].push_back ( l1err ); geoml2errs[ideg + 1].push_back ( l2err ); geomlinferrs[ideg + 1].push_back ( linferr );
        }
    }

    std::cout << "Mesh Size: ";

    for ( size_t i = 0; i < meshsizes.size(); ++i ) { std::cout << meshsizes[i] << " "; }

    std::cout << std::endl;
    std::cout << "Degrees: 0 ";

    for ( size_t ideg = 0; ideg < degs2fit.size(); ++ideg ) { std::cout << degs2fit[ideg] << " "; }

    std::cout << std::endl;
    std::cout << "L1-norm error: \n";

    for ( size_t i = 0; i < geoml1errs.size(); ++i )
    {
        for ( size_t j = 0; j < geoml1errs[i].size(); ++j )
        {
            std::cout << geoml1errs[i][j] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << "L2-norm error: \n";

    for ( size_t i = 0; i < geoml2errs.size(); ++i )
    {
        for ( size_t j = 0; j < geoml2errs[i].size(); ++j )
        {
            std::cout << geoml2errs[i][j] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << "Linf-norm error: \n";

    for ( size_t i = 0; i < geomlinferrs.size(); ++i )
    {
        for ( size_t j = 0; j < geomlinferrs[i].size(); ++j )
        {
            std::cout << geomlinferrs[i][j] << " ";
        }

        std::cout << std::endl;
    }

    return error;
}


void usage()
{
    std::cout << "usage: ./uref_hirec_test <mesh file> -degree <degree> -interp <0=least square, 1=interpolation> -dim <mesh dimension> -geom <s=sphere,t=torus>" << std::endl;
    std::cout << "Example: ./uref_hirec_test $MOAB_SOURCE_DIR/MeshFiles/unittest/sphere_tris_5.vtk -degree 2 -interp 0 -dim 2 -geom s" << std::endl;
}

int main ( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
    int nprocs, rank;
    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
#endif

    std::string infile = TestDir + "/sphere_tris_5.vtk";

    int degree = 2, dim = 2, geom = 0;
    bool interp = false;
    ErrorCode error;

    if (argc == 10)
    {
        infile = std::string(argv[1]); bool hasdim = false;

        for ( int i = 2; i < argc; ++i )
        {
            if ( i + 1 != argc )
            {
                if ( std::string ( argv[i] ) == "-degree" )
                {
                    degree = atoi ( argv[++i] );
                }
                else if ( std::string ( argv[i] ) == "-interp" )
                {
                    interp = atoi ( argv[++i] );
                }
                else if ( std::string ( argv[i] ) == "-dim" )
                {
                    dim = atoi ( argv[++i] ); hasdim = true;
                }
                else if ( std::string ( argv[i] ) == "-geom" )
                {
                    geom = std::string ( argv[++i] ) == "s" ? 0 : 1;
                }
                else
                {
#ifdef MOAB_HAVE_MPI

                    if ( 0 == rank )
                    {
                        usage();
                    }
                    MPI_Finalize();

#else
                    usage();
#endif
                    return -1;
                }
            }
        }

        if ( !hasdim )
        {
#ifdef MOAB_HAVE_MPI

            if ( 0 == rank )
            {
                std::cout << "Dimension of input mesh should be provided, positive and less than 3" << std::endl;
            }

#else
            std::cout << "Dimension of input mesh should be provided, positive and less than 3" << std::endl;
#endif
            return -1;
        }

        if ( degree <= 0 || dim > 2 || dim <= 0 )
        {
#ifdef MOAB_HAVE_MPI

            if ( 0 == rank )
            {
                std::cout << "Input degree should be positive number;\n";
                std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
            }

#else
            std::cout << "Input degree should be positive number;\n";
            std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
#endif
            return -1;
        }

#ifdef MOAB_HAVE_MPI

        if ( 0 == rank )
        {
            std::cout << "Testing on " << infile << " with dimension " << dim << "\n";
            std::string opts = interp ? "interpolation" : "least square fitting";
            std::cout << "High order reconstruction with degree " << degree << " " << opts << std::endl;
        }

#else
        std::cout << "Testing on " << infile << " with dimension " << dim << "\n";
        std::string opts = interp ? "interpolation" : "least square fitting";
        std::cout << "High order reconstruction with degree " << degree << " " << opts << std::endl;
#endif
    }
    else
    {
        if (argc > 1)
        {
            usage();
            return -1;
        }
    }

    {
        int level_degrees[3] = {2, 2, 2};
        int num_levels = 3;

        // create the geometry object
        geomObject* obj;
        if ( geom )
          obj = new torus();
        else
          obj = new sphere();

         error = test_closedsurface_mesh ( infile.c_str(), level_degrees, num_levels, degree, interp, dim, obj ); MB_CHK_ERR ( error );

        std::vector<int> degs2fit ( 6 );
        for ( int d = 1; d <= 6; ++d ) { degs2fit[d-1] = d; }

        // Call the higher order reconstruction routines and compute error convergence
        error = closedsurface_uref_hirec_convergence_study ( infile.c_str(), level_degrees, num_levels, degs2fit, interp, obj ); MB_CHK_ERR ( error );
        // cleanup memory
        delete obj;
    }

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

