#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "TestUtil.hpp"
#include "MBTagConventions.hpp"
#include <iostream>
#include <sstream>

// a file with 4 quads, in line, partitioned in 4 parts
std::string filename = TestDir + "/io/ln4.h5m";

using namespace moab;

void test_correct_ghost()
{
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;

  ErrorCode rval = MB_SUCCESS;

  // Get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);


  char read_opts[]="PARALLEL=READ_PART;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;PARALLEL_GHOSTS=2.0.1";
  rval = mb->load_file(filename.c_str() , 0, read_opts);CHECK_ERR(rval);

  if (nproc >= 3)
  {
    rval = pcomm-> correct_thin_ghost_layers();
    CHECK_ERR(rval);
  }

  rval = pcomm->exchange_ghost_cells(2, 0, 2, 0, true); CHECK_ERR(rval);// true to store remote handles

  // write in serial the database , on each rank
  std::ostringstream outfile;
  outfile <<"testReadThin_n" <<nproc<<"."<< rank<<".h5m";

  rval = mb->write_file(outfile.str().c_str()); // everything on local root
  CHECK_ERR(rval);
  delete mb;
}

void test_read_with_thin_ghost_layer()
{
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  moab::Core *mb = new moab::Core();

  ErrorCode rval = MB_SUCCESS;
  // Get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);

  char read_opts[]="PARALLEL=READ_PART;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;PARALLEL_GHOSTS=2.0.1;PARALLEL_THIN_GHOST_LAYER;";
  rval = mb->load_file(filename.c_str() , 0, read_opts);CHECK_ERR(rval);

  rval = pcomm->exchange_ghost_cells(2, 0, 2, 0, true); CHECK_ERR(rval);// true to store remote handles

  // write in serial the database , on each rank
  std::ostringstream outfile;
  outfile <<"testReadGhost_n" <<nproc<<"."<< rank<<".h5m";

  rval = mb->write_file(outfile.str().c_str()); // everything on local root
  CHECK_ERR(rval);
  delete mb;
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (nproc <= 3 && rank == 0)
  {
    std::cout << " launch it on at least 4 processes. \n";
    return 0;
  }

  int result = 0;
  if (argc>=2)
    filename = argv[1]; // to be able to test other files too

  result += RUN_TEST(test_read_with_thin_ghost_layer);
  result += RUN_TEST(test_correct_ghost);


  MPI_Finalize();
  return 0;
}
