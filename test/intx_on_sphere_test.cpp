/*
 * intx_on_sphere_test.cpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "moab/Intx2MeshOnSphere.hpp"
#include "moab/IntxUtils.hpp"
#include "TestUtil.hpp"
#include <math.h>

using namespace moab;

int main(int argc, char* argv[])
{
  int rank=0, size=1;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  bool output_fraction = true;
  bool write_files_rank = false;
  // check command line arg second grid is red, arrival, first mesh is blue, departure
  // will will keep the
  const char *filename_mesh1 = STRINGIFY(MESHDIR) "/mbcslam/lagrangeHomme.vtk";
  const char *filename_mesh2 = STRINGIFY(MESHDIR) "/mbcslam/eulerHomme.vtk";
  double R = 6. * sqrt(3.) / 2; // input
  double epsrel=1.e-8;
  double boxeps=0.1;
  const char *newFile = "intx.h5m";
  if (argc == 8)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    R = atof(argv[3]);
    epsrel = atof(argv[4]);
    newFile = argv[5];
    output_fraction = (bool) atoi(argv[6]); // could be 0 or 1, default true
    write_files_rank = (bool) atoi(argv[7]); // could be 0 or 1, false or true (default false)
  }
  else
  {
    printf("Usage: mpiexec -np <N> %s <mesh_filename1> <mesh_filename2> <radius> <epsrel> <newFile> <output_fraction (bool)> <write_files (bool)> \n",
        argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %f %f %s\n",
        filename_mesh1, filename_mesh2, R, epsrel, newFile);
  }

  std::string opts = (size == 1 ? "" : std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
        std::string(";PARALLEL_RESOLVE_SHARED_ENTS"));

  // read meshes in 2 file sets
  ErrorCode rval;
  Core moab;
  Interface * mb = &moab; // global
  EntityHandle sf1, sf2, outputSet;
   
  // create meshsets and load files
  if (0==rank)
    std::cout << "Creating mesh sets\n";
  rval = mb->create_meshset(MESHSET_SET, sf1);MB_CHK_ERR(rval);
  rval = mb->create_meshset(MESHSET_SET, sf2);MB_CHK_ERR(rval);
  if (0==rank)
    std::cout << "Loading mesh file 1\n";
  rval = mb->load_file(filename_mesh1, &sf1, opts.c_str());MB_CHK_ERR(rval);
  if (0==rank)
    std::cout << "Loading mesh file 2\n";
  rval = mb->load_file(filename_mesh2, &sf2, opts.c_str());MB_CHK_ERR(rval);

  rval = mb->create_meshset(MESHSET_SET, outputSet);MB_CHK_ERR(rval);

  // std::cout << "Fix orientation etc ..\n";
  //IntxUtils; those calls do nothing for a good mesh
  rval = fix_degenerate_quads(mb, sf1);MB_CHK_ERR(rval);
  rval = fix_degenerate_quads(mb, sf2);MB_CHK_ERR(rval);

  rval = positive_orientation(mb, sf1, R);MB_CHK_ERR(rval);
  rval = positive_orientation(mb, sf2, R);MB_CHK_ERR(rval);

  //ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
 
  Intx2MeshOnSphere  worker(mb);

  worker.SetErrorTolerance(R*epsrel);
  worker.set_box_error(boxeps);
  //worker.SetEntityType(moab::MBQUAD);
  worker.SetRadius(R);
  //worker.enable_debug();


  if (size>1)
  {
    Range local_verts;
    rval = worker.build_processor_euler_boxes(sf2, local_verts); MB_CHK_ERR(rval);// output also the local_verts
    if (write_files_rank)
    {
      std::stringstream outf;
      outf<<"second_mesh" << rank<<".h5m";
      rval = mb->write_file(outf.str().c_str(), 0, 0, &sf2, 1); MB_CHK_ERR(rval);
    }
  }
  EntityHandle covering_set;
  if (size>1)
  {
#ifdef MOAB_HAVE_MPI
    double elapsed = MPI_Wtime();
    rval = worker.construct_covering_set(sf1, covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
    elapsed = MPI_Wtime() - elapsed;
    if (0==rank)
      std::cout << "\nTime to communicate the mesh = " << elapsed << std::endl;
    // area fraction of the covering set that needed to be communicated from other processors
    // number of elements in the covering set communicated, compared to total number of elements in the covering set
    if (output_fraction)
    {
      EntityHandle comm_set; // set with elements communicated from other tasks
      rval = mb->create_meshset(MESHSET_SET, comm_set); MB_CHK_ERR(rval);
      // see how much more different is compared to sf1
      rval = mb->unite_meshset(comm_set, covering_set) ;  MB_CHK_ERR(rval); // will have to subtract from covering set, initial set
      // subtract
      rval = mb->subtract_meshset(comm_set, sf1); MB_CHK_ERR(rval);
      // compute fractions
      double area_cov_set = area_on_sphere(mb, covering_set, R);
      assert(area_cov_set>0);
      double comm_area = area_on_sphere(mb, comm_set, R);
      // more important is actually the number of elements communicated
      int num_cov_cells, num_comm_cells;
      rval = mb->get_number_entities_by_dimension(covering_set, 2, num_cov_cells); MB_CHK_ERR(rval);
      rval = mb->get_number_entities_by_dimension(comm_set, 2, num_comm_cells); MB_CHK_ERR(rval);
      double fraction_area = comm_area/area_cov_set;
      double fraction_num_cells = (double) num_comm_cells/num_cov_cells; // determine min, max, average of these fractions

      double max_fraction_area, max_fraction_num_cells, min_fraction_area, min_fraction_num_cells;
      double average_fraction_area, average_fraction_num_cells;
      MPI_Reduce(&fraction_area, &max_fraction_area, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &max_fraction_num_cells, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_area, &min_fraction_area, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &min_fraction_num_cells, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_area, &average_fraction_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &average_fraction_num_cells, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      average_fraction_area /= size;
      average_fraction_num_cells /= size;
      if (rank==0)
      {
        std::cout << " fraction area:      min: " << min_fraction_area << " max: " << max_fraction_area << " average :" << average_fraction_area << " \n";
        std::cout << " fraction num cells: min: " << min_fraction_num_cells << " max: " << max_fraction_num_cells << " average: " <<average_fraction_num_cells <<  " \n";
      }

    }
    if (write_files_rank)
    {
      std::stringstream outf, cof;
      outf<<"first_mesh" << rank<<".h5m";
      rval = mb->write_file(outf.str().c_str(), 0, 0, &sf1, 1); MB_CHK_ERR(rval);
      cof<<"covering_mesh" << rank<<".h5m";
      rval = mb->write_file(cof.str().c_str(), 0, 0, &covering_set, 1); MB_CHK_ERR(rval);
    }


#endif
  }
  else
    covering_set = sf1;

  if (0==rank)
    std::cout << "Computing intersections ..\n";
#ifdef MOAB_HAVE_MPI
  double elapsed = MPI_Wtime();
#endif
  rval = worker.intersect_meshes(covering_set, sf2, outputSet);MB_CHK_SET_ERR(rval,"failed to intersect meshes");
#ifdef MOAB_HAVE_MPI
  elapsed = MPI_Wtime() - elapsed;
  if (0==rank)
    std::cout << "\nTime to compute the intersection between meshes = " << elapsed << std::endl;
#endif
  // the output set does not have the intx vertices on the boundary shared, so they will be duplicated right now
  // we write this file just for checking it looks OK

  //compute total area with 2 methods
  //double initial_area = area_on_sphere_lHuiller(mb, sf1, R);
  //double area_method1 = area_on_sphere_lHuiller(mb, outputSet, R);
  //double area_method2 = area_on_sphere(mb, outputSet, R);

  //std::cout << "initial area: " << initial_area << "\n";
  //std::cout<< " area with l'Huiller: " << area_method1 << " with Girard: " << area_method2<< "\n";
  //std::cout << " relative difference areas " << fabs(area_method1-area_method2)/area_method1 << "\n";
  //std::cout << " relative error " << fabs(area_method1-initial_area)/area_method1 << "\n";

  if (write_files_rank)
  {
    std::stringstream outf;
    outf<<"intersect" << rank<<".h5m";
    rval = mb->write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  }
  double intx_area = area_on_sphere(mb, outputSet, R);
  double arrival_area = area_on_sphere(mb, sf2, R) ;
  std::cout<< "On rank : " << rank << " arrival area: " << arrival_area<<
      "  intersection area:" << intx_area << " rel error: " << fabs((intx_area-arrival_area)/arrival_area) << "\n";

 // rval = mb->write_file(newFile, 0, "PARALLEL=WRITE_PART", &outputSet, 1);MB_CHK_SET_ERR(rval,"failed to write intx file");

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;

}


