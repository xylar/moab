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
#include "moab/ParallelComm.hpp"
#include "moab/Intx2MeshOnSphere.hpp"
#include "moab/IntxUtils.hpp"
#include "TestUtil.hpp"
#include <math.h>

using namespace moab;

int main(int argc, char* argv[])
{
  // check command line arg// Euler grid is red, arrival, Lagrangian is blue, departure
  // will will keep the
  const char *filename_mesh1 = STRINGIFY(MESHDIR) "/mbcslam/lagrangeHomme.vtk";
  const char *filename_mesh2 = STRINGIFY(MESHDIR) "/mbcslam/eulerHomme.vtk";
  double R = 6. * sqrt(3.) / 2; // input
  double epsrel=1.e-8;
  double boxeps=0.1;
  const char *newFile = "intx.h5m";
  if (argc == 6)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    R = atof(argv[3]);
    epsrel = atof(argv[4]);
    newFile = argv[5];
  }
  else
  {
    printf("Usage: %s <mesh_filename1> <mesh_filename2> <radius> <epsrel> <newFile>\n",
        argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %f %f %s\n",
        filename_mesh1, filename_mesh2, R, epsrel, newFile);
  }

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::string opts = (size == 1 ? "" : std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
        std::string(";PARALLEL_RESOLVE_SHARED_ENTS"));

  // read meshes in 2 file sets
  ErrorCode rval;
  Core moab;
  Interface * mb = &moab; // global
  EntityHandle sf1, sf2, outputSet;
   
  // create meshsets and load files
  std::cout << "Creating mesh sets\n";
  rval = mb->create_meshset(MESHSET_SET, sf1);MB_CHK_ERR(rval);
  rval = mb->create_meshset(MESHSET_SET, sf2);MB_CHK_ERR(rval);
  std::cout << "Loading mesh file 1\n";
  rval = mb->load_file(filename_mesh1, &sf1, opts.c_str());MB_CHK_ERR(rval);
  std::cout << "Loading mesh file 2\n";
  rval = mb->load_file(filename_mesh2, &sf2, opts.c_str());MB_CHK_ERR(rval);

  rval = mb->create_meshset(MESHSET_SET, outputSet);MB_CHK_ERR(rval);

  std::cout << "Fix orientation etc ..\n";
  // IntxUtils; those calls do nothing for a good mesh
  //rval = fix_degenerate_quads(mb, sf1);MB_CHK_ERR(rval);
  //rval = fix_degenerate_quads(mb, sf2);MB_CHK_ERR(rval);

  //rval = positive_orientation(mb, sf1, R);MB_CHK_ERR(rval);
  //rval = positive_orientation(mb, sf2, R);MB_CHK_ERR(rval);

  //ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
 
  Intx2MeshOnSphere  worker(mb);

  worker.SetErrorTolerance(R*epsrel);
  worker.set_box_error(boxeps);
  //worker.SetEntityType(moab::MBQUAD);
  worker.SetRadius(R);
  //worker.enable_debug();

  Range local_verts;
  if (size>1)
  {
    rval = worker.build_processor_euler_boxes(sf1, local_verts); MB_CHK_ERR(rval);// output also the local_verts
  }
  EntityHandle covering_set;
  if (size>1)
  {
    rval = worker.create_departure_mesh_3rd_alg(sf2, covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
  }
  else
    covering_set = sf2;

  std::cout << "Computing intersections ..\n";
  double elapsed = MPI_Wtime();
  rval = worker.intersect_meshes(sf1, covering_set, outputSet);MB_CHK_SET_ERR(rval,"failed to intersect meshes");
  elapsed = MPI_Wtime() - elapsed;
  std::cout << "\nTime to compute the intersection between meshes = " << elapsed << std::endl; 

  // the output set does not have the intx vertices on the boundary shared, so they will be duplicated right now
  // we write this file just for checking it looks OK
  rval = mb->write_file(newFile, 0, "PARALLEL=WRITE_PART", &outputSet, 1);MB_CHK_SET_ERR(rval,"failed to write intx file");


  //compute total area with 2 methods
  //double initial_area = area_on_sphere_lHuiller(mb, sf1, R);
  //double area_method1 = area_on_sphere_lHuiller(mb, outputSet, R);
  //double area_method2 = area_on_sphere(mb, outputSet, R);

  //std::cout << "initial area: " << initial_area << "\n";
  //std::cout<< " area with l'Huiller: " << area_method1 << " with Girard: " << area_method2<< "\n";
  //std::cout << " relative difference areas " << fabs(area_method1-area_method2)/area_method1 << "\n";
  //std::cout << " relative error " << fabs(area_method1-initial_area)/area_method1 << "\n";

  std::stringstream outf;
  outf<<"intersect" << rank<<".h5m";
  rval = mb->write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  double intx_area = area_on_sphere(mb, outputSet, R);
  double arrival_area = area_on_sphere(mb, sf1, R) ;
  std::cout<< "On rank : " << rank << " arrival area: " << arrival_area<<
      "  intersection area:" << intx_area << " rel error: " << fabs((intx_area-arrival_area)/arrival_area) << "\n";
  CHECK_ERR(rval);

  MPI_Finalize();
  return 0;

}


