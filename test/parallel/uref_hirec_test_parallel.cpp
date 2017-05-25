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
#include "moab/DiscreteGeometry/HiReconstruction.hpp"
#include "TestUtil.hpp"
#include "geomObject.cpp"
#include <math.h>
#include <stdlib.h>
#include <assert.h>

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

ErrorCode load_meshset_hirec ( const char* infile, Interface* mbimpl, EntityHandle& meshset, ParallelComm*& pc, const int degree, const int dim )
{
    ErrorCode error;
    error = mbimpl->create_meshset ( moab::MESHSET_SET, meshset ); MB_CHK_ERR ( error );
#ifdef MOAB_HAVE_MPI
    int nprocs, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size ( comm, &nprocs );
    MPI_Comm_rank ( comm, &rank );
    EntityHandle partnset;
    error = mbimpl->create_meshset ( moab::MESHSET_SET, partnset ); MB_CHK_ERR ( error );

    if ( nprocs > 1 )
    { pc = moab::ParallelComm::get_pcomm ( mbimpl, partnset, &comm ); }

    if ( nprocs > 1 )
    {
        int nghlayers = degree > 0 ? HiReconstruction::estimate_num_ghost_layers ( degree, true ) : 0; assert ( nghlayers < 10 );

        if ( nghlayers )
        {
            //get ghost layers
            if ( dim == 2 )
            {
                read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.";
            }
            else if ( dim == 1 )
            {
                read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=1.0.";
            }
            else
            {
                MB_SET_ERR ( MB_FAILURE, "unsupported dimension" );
            }

            read_options += ( char ) ( '0' + nghlayers );
            //std::cout << "On processor " << rank << " Degree=" << degree << " "<< read_options << std::endl;
        }
        else
        {
            read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        }

        /*for debug*/
        /*std::string outfile = std::string(infile); std::size_t dotpos = outfile.find_last_of(".");
        std::ostringstream convert; convert << rank;
        std::string localfile = outfile.substr(0,dotpos)+convert.str()+".vtk";
        //write local mesh
        EntityHandle local;
        error = mbimpl->create_meshset(moab::MESHSET_SET,local);CHECK_ERR(error);
        std::string local_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        error = mbimpl->load_file(infile,&local,local_options.c_str()); MB_CHK_ERR(error);
        error = mbimpl->write_file(localfile.c_str(),0,0,&local,1);*/

        error = mbimpl->load_file ( infile, &meshset, read_options.c_str() ); MB_CHK_ERR ( error );
        /*for debug
        //write local mesh with ghost layers
        localfile = outfile.substr(0,dotpos)+convert.str()+"_ghost.vtk";
        error = mbimpl->write_file(localfile.c_str(),0,0,&meshset,1);*/

    }
    else
    {
        error = mbimpl->load_file ( infile, &meshset ); MB_CHK_ERR ( error );
    }

#else
    assert ( !pc && degree && dim );
    error = mbimpl->load_file ( infile, &meshset ); MB_CHK_ERR ( error );
#endif
    return error;
}

ErrorCode closedsurface_uref_hirec_convergence_study ( const char* infile, std::vector<int>& degs2fit, bool interp, int dim, geomObject* obj, int& ntestverts, std::vector<double>& geoml1errs, std::vector<double>& geoml2errs, std::vector<double>& geomlinferrs )
{
    Core moab;
    Interface* mbImpl = &moab;
    ParallelComm* pc = NULL;
    EntityHandle meshset;

#ifdef MOAB_HAVE_MPI
    int nprocs, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size ( comm, &nprocs );
    MPI_Comm_rank ( comm, &rank );
#endif

    ErrorCode error;
    //mesh will be loaded and communicator pc will be updated
    int mxdeg = 1; for ( size_t i = 0; i < degs2fit.size(); ++i ) { mxdeg = std::max ( degs2fit[i], mxdeg ); }
    error = load_meshset_hirec ( infile, mbImpl, meshset, pc, mxdeg, dim ); MB_CHK_ERR ( error );

    Range elems, elems_owned;
    error = mbImpl->get_entities_by_dimension ( meshset, dim, elems ); MB_CHK_ERR ( error );

#ifdef MOAB_HAVE_MPI

    if ( pc )
    {
        error = pc->filter_pstatus ( elems, PSTATUS_GHOST, PSTATUS_NOT, -1, &elems_owned ); MB_CHK_ERR ( error );
    }
    else
    {
        elems_owned = elems;
    }

#endif

#ifdef MOAB_HAVE_MPI
    //std::cout << "Mesh has " << nelems << " elements on Processor " << rank << " in total;";
    //std::cout << elems_owned.size() << " of which are locally owned elements" << std::endl;
#else
    //std::cout << "Mesh has " << nelems << " elements" << std::endl;
#endif

    /************************
    *   convergence study   *
    *************************/
    //project onto exact geometry since each level with uref has only linear coordinates
    Range verts;
    error = mbImpl->get_entities_by_dimension ( meshset, 0, verts ); MB_CHK_ERR ( error );

    for ( Range::iterator ivert = verts.begin(); ivert != verts.end(); ++ivert )
    {
        EntityHandle currvert = *ivert;
        double currcoords[3], exactcoords[3];
        error = mbImpl->get_coords ( &currvert, 1, currcoords ); MB_CHK_ERR ( error );
        obj->project_points2geom ( 3, currcoords, exactcoords, NULL );

        error = mbImpl->set_coords ( &currvert, 1, exactcoords ); MB_CHK_ERR ( error );
        //for debug
        /*error = mbImpl->get_coords(&currvert,1,currcoords); MB_CHK_ERR(error);
        assert(currcoords[0]==exactcoords[0]&&currcoords[1]==exactcoords[1]&&currcoords[2]==exactcoords[2]);*/
    }

    //test ahf on mesh with ghost layers
    /*Range verts_owned;
    #ifdef MOAB_HAVE_MPI
        if(pc){
            error = pc->filter_pstatus(verts,PSTATUS_GHOST,PSTATUS_NOT,-1,&verts_owned); MB_CHK_ERR(error);
        }else{
            verts_owned = verts;
        }

        HalfFacetRep *ahf = new HalfFacetRep(&moab,pc,meshset);
        error = ahf->initialize(); MB_CHK_ERR(error);

        for(int i=0;i<nprocs;++i){
            if(rank==i){
                std::cout << "Processor " << rank << " Local elements: \n";
                for(Range::iterator ielem=elems.begin();ielem!=elems.end();++ielem){
                    EntityHandle currelem = *ielem;
                    std::vector<EntityHandle> conn;
                    error = mbImpl->get_connectivity(&currelem,1,conn); MB_CHK_ERR(error);
                    std::cout << *ielem << ": ";
                    for(size_t k=0;k<conn.size();++k) std::cout << conn[k] << " ";
                    std::cout << std::endl;
                }

                std::cout << "Processor " << rank << " Local owned elements: \n";
                for(Range::iterator ielem=elems_owned.begin();ielem!=elems_owned.end();++ielem){
                    EntityHandle currelem = *ielem;
                    std::vector<EntityHandle> conn;
                    error = mbImpl->get_connectivity(&currelem,1,conn); MB_CHK_ERR(error);
                    std::cout << *ielem << ": ";
                    for(size_t k=0;k<conn.size();++k) std::cout << conn[k] << " ";
                    std::cout << std::endl;
                }


                std::cout << "Processor " << rank << " Local vertices: ";
                for(Range::iterator ivert=verts.begin();ivert!=verts.end();++ivert){
                    std::cout << *ivert << " ";
                }
                std::cout << std::endl;
                std::cout << "Processor " << rank << " Local owned vertices: ";
                for(Range::iterator ivert=verts_owned.begin();ivert!=verts_owned.end();++ivert){
                    std::cout << *ivert << " ";
                }
                std::cout << std::endl;

                for(Range::iterator ivert=verts_owned.begin();ivert!=verts_owned.end();++ivert){
                    EntityHandle currvid = *ivert;
                    std::cout << "Processor " << rank << " local verts: " << *ivert << " has adjfaces: ";
                    std::vector<EntityHandle> adjfaces;
                    //error = ahf->get_up_adjacencies(currvid,2,adjfaces); MB_CHK_ERR(error);
                    error = mbImpl->get_adjacencies(&currvid,1,2,false,adjfaces); MB_CHK_ERR(error);
                    for(size_t k=0;k<adjfaces.size();++k){
                        std::cout << adjfaces[k] << " ";
                    }
                    std::cout << std::endl;
                }

                for(Range::iterator ivert=verts.begin();ivert!=verts.end();++ivert){
                    EntityHandle currvid = *ivert;
                    std::cout << "Processor " << rank << " all verts: " << *ivert << " has adjfaces: ";
                    std::vector<EntityHandle> adjfaces;
                    //error = ahf->get_up_adjacencies(*ivert,2,adjfaces); MB_CHK_ERR(error);
                    error = mbImpl->get_adjacencies(&currvid,1,2,false,adjfaces); MB_CHK_ERR(error);
                    for(size_t k=0;k<adjfaces.size();++k){
                        std::cout << adjfaces[k] << " ";
                    }
                    std::cout << std::endl;
                }

            }else{
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        exit(0);
    #endif*/

    //generate random points on each elements, assument 3D coordinates
    int nvpe = TYPE_FROM_HANDLE ( *elems.begin() ) == MBTRI ? 3 : 4;
    std::vector<double> testpnts, testnaturalcoords; ntestverts = elems_owned.size() * nsamples;
    testpnts.reserve ( 3 * elems_owned.size() *nsamples ); testnaturalcoords.reserve ( nvpe * elems_owned.size() *nsamples );

    for ( Range::iterator ielem = elems_owned.begin(); ielem != elems_owned.end(); ++ielem )
    {
        EntityHandle currelem = *ielem;
        std::vector<EntityHandle> conn;
        error = mbImpl->get_connectivity ( &currelem, 1, conn ); MB_CHK_ERR ( error );
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

                assert ( a > 0 && b > 0 && c > 0 && fabs ( a + b + c - 1 ) < 1e-12 );
                testpnts.push_back ( a * elemcoords[0] + b * elemcoords[3] + c * elemcoords[6] );
                testpnts.push_back ( a * elemcoords[1] + b * elemcoords[4] + c * elemcoords[7] );
                testpnts.push_back ( a * elemcoords[2] + b * elemcoords[5] + c * elemcoords[8] );
                testnaturalcoords.push_back ( a ), testnaturalcoords.push_back ( b ), testnaturalcoords.push_back ( c );
            }
            else if ( type == MBQUAD )
            {
                double xi = ( double ) rand() / RAND_MAX, eta = ( double ) rand() / RAND_MAX, N[4];
                xi = 2 * xi - 1; eta = 2 * eta - 1;
                N[0] = ( 1 - xi ) * ( 1 - eta ) / 4, N[1] = ( 1 + xi ) * ( 1 - eta ) / 4, N[2] = ( 1 + xi ) * ( 1 + eta ) / 4, N[3] = ( 1 - xi ) * ( 1 + eta ) / 4;
                testpnts.push_back ( N[0]*elemcoords[0] + N[1]*elemcoords[3] + N[2]*elemcoords[6] + N[3]*elemcoords[9] );
                testpnts.push_back ( N[0]*elemcoords[1] + N[1]*elemcoords[4] + N[2]*elemcoords[7] + N[3]*elemcoords[10] );
                testpnts.push_back ( N[0]*elemcoords[2] + N[1]*elemcoords[5] + N[2]*elemcoords[8] + N[3]*elemcoords[11] );
                testnaturalcoords.push_back ( N[0] ), testnaturalcoords.push_back ( N[1] ), testnaturalcoords.push_back ( N[2] ), testnaturalcoords.push_back ( N[3] );
            }
            else
            {
                throw std::invalid_argument ( "NOT SUPPORTED ELEMENT TYPE" );
            }
        }
    }

    //Compute error of linear interpolation in PARALLEL
    double l1err, l2err, linferr;
    geoml1errs.clear(); geoml2errs.clear(); geomlinferrs.clear();
    obj->compute_projecterror ( 3, elems_owned.size() *nsamples, & ( testpnts[0] ), l1err, l2err, linferr );
    geoml1errs.push_back ( l1err ); geoml2errs.push_back ( l2err ); geomlinferrs.push_back ( linferr );
    /*Perform high order projection and compute error*/
    //initialize
    HiReconstruction hirec ( &moab, pc, meshset );

    for ( size_t ideg = 0; ideg < degs2fit.size(); ++ideg )
    {
        //High order reconstruction
        error = hirec.reconstruct3D_surf_geom ( degs2fit[ideg], interp, true, true ); MB_CHK_ERR ( error );

        int index = 0;
        //for debug
        /*double maxlinferr=0,elemlinferr=0,eleml1err=0,eleml2err;
        EntityHandle elem_maxerr;*/

        for ( Range::iterator ielem = elems_owned.begin(); ielem != elems_owned.end(); ++ielem, ++index )
        {
            //Projection
            error = hirec.hiproj_walf_in_element ( *ielem, nvpe, nsamples, & ( testnaturalcoords[nvpe * nsamples * index] ), & ( testpnts[3 * nsamples * index] ) ); MB_CHK_ERR ( error );
            //for debug
            /*obj->compute_projecterror(3,nsamples,&(testpnts[3*nsamples*index]),eleml1err,eleml2err,elemlinferr);
            if(elemlinferr>maxlinferr){
                maxlinferr = elemlinferr; elem_maxerr = *ielem;
            }*/
        }

        assert ( ( size_t ) index == elems_owned.size() );
        //Estimate error
        obj->compute_projecterror ( 3, elems_owned.size() *nsamples, & ( testpnts[0] ), l1err, l2err, linferr );
        geoml1errs.push_back ( l1err ); geoml2errs.push_back ( l2err ); geomlinferrs.push_back ( linferr );

        //for debug
        /*std::cout << "Degree " << degs2fit[ideg] << " has L-inf error " << maxlinferr << " at element " << elems_owned.index(elem_maxerr) << ":";
        std::vector<EntityHandle> conn;
        error = mbImpl->get_connectivity(&elem_maxerr,1,conn); MB_CHK_ERR(error);
        for(size_t ii=0;ii<conn.size();++ii) std::cout << verts.index(conn[ii]) << " "; std::cout << std::endl;

        for(size_t ii=0;ii<conn.size();++ii){
            //assume triangle
            std::vector<double> vertexcoords(3,0);
            error = mbImpl->get_coords(&(conn[ii]),1,&(vertexcoords[0])); MB_CHK_ERR(error);
            std::cout << verts.index(conn[ii]) << ":" << vertexcoords[0] << " "<< vertexcoords[1] << " "<< vertexcoords[2] << "\n";
            GEOMTYPE geomtype;
            std::vector<double> local_coords_system, local_coeffs;
            int deg_out; bool local_interp;
            bool hasfit = hirec.get_fittings_data(conn[ii],geomtype,local_coords_system,deg_out,local_coeffs,local_interp); assert(hasfit);
            std::cout << "Local fitting result: " << std::endl;
            std::cout << "Geom Type: " << geomtype << std::endl;
            std::cout << "Local coordinates system: " << std::endl;
            std::cout << local_coords_system[0] << "\t" << local_coords_system[1] << "\t" << local_coords_system[2] << std::endl;
            std::cout << local_coords_system[3] << "\t" << local_coords_system[4] << "\t" << local_coords_system[5] << std::endl;
            std::cout << local_coords_system[6] << "\t" << local_coords_system[7] << "\t" << local_coords_system[8] << std::endl;
            std::cout << "Fitting degree is " << deg_out << " interp is " << local_interp << std::endl;
            std::cout << "Fitting coefficients are :" << std::endl;
            for(size_t kk=0;kk<local_coeffs.size();++kk) std::cout << local_coeffs[kk] << "\t";
            std::cout << std::endl;
        }*/
    }

    return error;
}


void usage()
{
    std::cout << "usage: mpirun -np <number of processors> ./uref_hirec_test_parallel <prefix of files> -str <first number> -end <lastnumber> -suffix <suffix> -interp <0=least square, 1=interpolation> -dim <mesh dimension> -geom <s=sphere,t=torus>" << std::endl;
}

int main ( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
    int nprocs, rank;
    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
#endif
    std::string prefix, suffix, istr, iend;
    int dim = 2, geom = 0;
    bool interp = false;
    ErrorCode error;

    if ( argc == 1 )
    {
        usage();
        prefix = TestDir + "/sphere_quads_5_ML_";
        suffix = ".h5m";
        istr = "0";
        iend = "1";
        std::cout << "Using defaults: MeshFiles/unittest/sphere_quads_5_ML_*.h5m -str 0 -end 1 -suffix \".h5m\" -interp 0 -dim 2 -geom s" << std::endl;
    }
    else
    {
        prefix = std::string ( argv[1] ); bool hasdim = false;

        for ( int i = 2; i < argc; ++i )
        {
            if ( i + 1 != argc )
            {
                if ( std::string ( argv[i] ) == "-suffix" )
                {
                    suffix = std::string ( argv[++i] );
                }
                else if ( std::string ( argv[i] ) == "-str" )
                {
                    istr = std::string ( argv[++i] );
                }
                else if ( std::string ( argv[i] ) == "-end" )
                {
                    iend = std::string ( argv[++i] );
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
                    return 0;
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
            return 0;
        }

        if ( dim > 2 || dim <= 0 )
        {
#ifdef MOAB_HAVE_MPI

            if ( 0 == rank )
            {
                std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
            }

#else
            std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
#endif
            return 0;
        }

#ifdef MOAB_HAVE_MPI

        if ( 0 == rank )
        {
            std::cout << "Testing on " << prefix + istr + "-" + iend + suffix << " with dimension " << dim << "\n";
            std::string opts = interp ? "interpolation" : "least square fitting";
            std::cout << "High order reconstruction in " << opts << std::endl;
        }

#else
        std::cout << "Testing on " << prefix + istr + "-" + iend + suffix << " with dimension " << dim << "\n";
        std::string opts = interp ? "interpolation" : "least square fitting";
        std::cout << "High order reconstruction in " << opts << std::endl;
#endif
    }

    geomObject* obj;

    if ( geom == 0 )
    {
        obj = new sphere();
    }
    else
    {
        obj = new torus();
    }

    std::vector<int> degs2fit; for ( int d = 1; d <= 6; ++d ) { degs2fit.push_back ( d ); }

    int begin = atoi ( istr.c_str() ), end = atoi ( iend.c_str() ) + 1;
#ifdef MOAB_HAVE_MPI
    std::vector< std::vector<double> > geoml1errs_global, geoml2errs_global, geomlinferrs_global;

    if ( rank == 0 )
    {
        geoml1errs_global = std::vector< std::vector<double> > ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
        geoml2errs_global = std::vector< std::vector<double> > ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
        geomlinferrs_global = std::vector< std::vector<double> > ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
    }

#else
    std::vector< std::vector<double> > geoml1errs_global ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
    std::vector< std::vector<double> > geoml2errs_global ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
    std::vector< std::vector<double> > geomlinferrs_global ( 1 + degs2fit.size(), std::vector<double> ( end - begin, 0 ) );
#endif

    for ( int i = begin; i < end; ++i )
    {
        std::ostringstream convert; convert << i;
        std::string infile = prefix + convert.str() + suffix;
        int ntestverts; std::vector<double> geoml1errs, geoml2errs, geomlinferrs;

#ifdef MOAB_HAVE_MPI
        std::cout << "Processor " << rank << " is working on file " << infile << std::endl;
#endif
        error = closedsurface_uref_hirec_convergence_study ( infile.c_str(), degs2fit, interp, dim, obj, ntestverts, geoml1errs, geoml2errs, geomlinferrs ); MB_CHK_ERR ( error );
        assert ( geoml1errs.size() == 1 + degs2fit.size() && geoml2errs.size() == 1 + degs2fit.size() && geomlinferrs.size() == 1 + degs2fit.size() );
#ifdef MOAB_HAVE_MPI

        if ( nprocs > 1 )
        {
            int ntestverts_global = 0;
            MPI_Reduce ( &ntestverts, &ntestverts_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

            for ( size_t d = 0; d < degs2fit.size() + 1; ++d )
            {
                double local_l1err = ntestverts * geoml1errs[d], local_l2err = geoml2errs[d] * ( geoml2errs[d] * ntestverts ), local_linferr = geomlinferrs[d];
                std::cout << "On Processor " << rank << " with mesh " << i << " Degree = " << ( d == 0 ? 0 : degs2fit[d - 1] ) << " L1:" << geoml1errs[d] << " L2:" << geoml2errs[d] << " Li:" << geomlinferrs[d] << std::endl;
                double global_l1err = 0, global_l2err = 0, global_linferr = 0;
                MPI_Reduce ( &local_l1err, &global_l1err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
                MPI_Reduce ( &local_l2err, &global_l2err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
                MPI_Reduce ( &local_linferr, &global_linferr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

                if ( rank == 0 )
                {
                    geoml1errs_global[d][i - begin] = global_l1err / ntestverts_global;
                    geoml2errs_global[d][i - begin] = sqrt ( global_l2err / ntestverts_global );
                    geomlinferrs_global[d][i - begin] = global_linferr;
                }
            }
        }
        else
        {
            for ( size_t d = 0; d < degs2fit.size() + 1; ++d )
            {
                geoml1errs_global[d][i - begin] = geoml1errs[d];
                geoml2errs_global[d][i - begin] = geoml2errs[d];
                geomlinferrs_global[d][i - begin] = geomlinferrs[d];
            }
        }

#else

        for ( size_t d = 0; d < degs2fit.size() + 1; ++d )
        {
            geoml1errs_global[d][i - begin] = geoml1errs[d];
            geoml2errs_global[d][i - begin] = geoml2errs[d];
            geomlinferrs_global[d][i - begin] = geomlinferrs[d];
        }

#endif
    }

#ifdef MOAB_HAVE_MPI

    if ( rank == 0 )
    {
        std::cout << "Degrees: 0 ";

        for ( size_t ideg = 0; ideg < degs2fit.size(); ++ideg ) { std::cout << degs2fit[ideg] << " "; }

        std::cout << std::endl;
        std::cout << "L1-norm error: \n";

        for ( size_t i = 0; i < geoml1errs_global.size(); ++i )
        {
            for ( size_t j = 0; j < geoml1errs_global[i].size(); ++j )
            {
                std::cout << geoml1errs_global[i][j] << " ";
            }

            std::cout << std::endl;
        }

        std::cout << "L2-norm error: \n";

        for ( size_t i = 0; i < geoml2errs_global.size(); ++i )
        {
            for ( size_t j = 0; j < geoml2errs_global[i].size(); ++j )
            {
                std::cout << geoml2errs_global[i][j] << " ";
            }

            std::cout << std::endl;
        }

        std::cout << "Linf-norm error: \n";

        for ( size_t i = 0; i < geomlinferrs_global.size(); ++i )
        {
            for ( size_t j = 0; j < geomlinferrs_global[i].size(); ++j )
            {
                std::cout << geomlinferrs_global[i][j] << " ";
            }

            std::cout << std::endl;
        }
    }

#else
    std::cout << "Degrees: 0 ";

    for ( size_t ideg = 0; ideg < degs2fit.size(); ++ideg ) { std::cout << degs2fit[ideg] << " "; }

    std::cout << std::endl;
    std::cout << "L1-norm error: \n";

    for ( size_t i = 0; i < geoml1errs_global.size(); ++i )
    {
        for ( size_t j = 0; j < geoml1errs_global[i].size(); ++j )
        {
            std::cout << geoml1errs_global[i][j] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << "L2-norm error: \n";

    for ( size_t i = 0; i < geoml2errs_global.size(); ++i )
    {
        for ( size_t j = 0; j < geoml2errs_global[i].size(); ++j )
        {
            std::cout << geoml2errs_global[i][j] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << "Linf-norm error: \n";

    for ( size_t i = 0; i < geomlinferrs_global.size(); ++i )
    {
        for ( size_t j = 0; j < geomlinferrs_global[i].size(); ++j )
        {
            std::cout << geomlinferrs_global[i][j] << " ";
        }

        std::cout << std::endl;
    }

#endif
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
}

