/*
 * migrate_nontrivial.cpp
 *
 *  migrate_nontrivial will contain tests for migrating meshes in parallel environments, with iMOAB api,
 *  using nontrivial partitions; for starters, we will use zoltan to compute partition in parallel,
 *    and then use the result to migrate the mesh
 *  trivial partition gave terrible results for mpas type meshes, that were numbered
 *  like a fractal; the resulting migrations were like Swiss cheese, full of holes
 *  a mesh is read on senders tasks
 *  we will use graph like methods or geometry methods, and will modify the ZoltanPartitioner
 *   to add needed methods
 *
 *  mesh will be sent to receivers tasks, with nonblocking MPI_Isend calls, and then received
 *    with blocking MPI_Recv calls;
 *
 *  we will not modify the GLOBAL_ID tag, we assume it was set correctly before we started
 */

#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/ProgOptions.hpp"

#define RUN_TEST_ARG2(A, B) run_test( &A, #A, B)

using namespace moab;

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s\n", message); return MB_FAILURE;}

int is_any_proc_error( int is_my_error )
{
  int result = 0;
  int err = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  return err || result;
}

int run_test( ErrorCode (*func)(const char*),
              const char* func_name,
              const char* file_name)
{
  ErrorCode result = (*func)(file_name);
  int is_err = is_any_proc_error( (MB_SUCCESS != result) );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (rank == 0) {
    if (is_err)
      std::cout << func_name << " : FAILED!!" << std::endl;
    else
      std::cout << func_name << " : success" << std::endl;
  }

  return is_err;
}

ErrorCode migrate_graph( const char* filename );
ErrorCode migrate_geom( const char* filename );

// some global variables, used by all tests
int rank, size, ierr;

int compid1, compid2;  // component ids are unique over all pes, and established in advance;
int nghlay; // number of ghost layers for loading the file
std::vector<int> groupTasks; // at most 4 tasks
int startG1, startG2, endG1, endG2;

MPI_Comm jcomm; // will be a copy of the global
MPI_Group jgroup;

ErrorCode migrate_smart(const char*filename, const char * outfile, int partMethod)
{
  // first create MPI groups

  std::string filen(filename);
  MPI_Group group1, group2;
  groupTasks.resize(endG1-startG1+1);
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, &groupTasks[0], &group1);
  CHECKRC(ierr, "can't create group1")

  groupTasks.resize(endG2-startG2+1);
  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, &groupTasks[0], &group2);
  CHECKRC(ierr, "can't create group2")

  // create 2 communicators, one for each group
  int tagcomm1 = 1, tagcomm2 = 2;
  MPI_Comm comm1, comm2;
  ierr = MPI_Comm_create_group(jcomm, group1, tagcomm1, &comm1);
  CHECKRC(ierr, "can't create comm1")

  ierr = MPI_Comm_create_group(jcomm, group2, tagcomm2, &comm2);
  CHECKRC(ierr, "can't create comm2")

  ierr = iMOAB_Initialize(0, 0); // not really needed anything from argc, argv, yet; maybe we should
  CHECKRC(ierr, "can't initialize iMOAB")

  // give some dummy values to component ids, just to differentiate between them
  // the par comm graph is unique between components
  compid1 = 4;
  compid2 = 7;

  int appID1;
  iMOAB_AppID pid1=&appID1;
  int appID2;
  iMOAB_AppID pid2=&appID2;

  if (comm1 != MPI_COMM_NULL) {
    ierr = iMOAB_RegisterApplication("APP1", &comm1, &compid1, pid1);
    CHECKRC(ierr, "can't register app1 ")
  }
  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_RegisterApplication("APP2", &comm2, &compid2, pid2);
    CHECKRC(ierr, "can't register app2 ")
  }

  if (comm1 != MPI_COMM_NULL) {

      std::string   readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");

      nghlay = 0;

      ierr = iMOAB_LoadMesh(pid1, filen.c_str(), readopts.c_str(), &nghlay, filen.length(), strlen(readopts.c_str()) );
      CHECKRC(ierr, "can't load mesh ")
      ierr = iMOAB_SendMesh(pid1, &jcomm, &group2, &compid2, &partMethod); // send to component 2
      CHECKRC(ierr, "cannot send elements" )
  }

  if (comm2 != MPI_COMM_NULL) {
     ierr = iMOAB_ReceiveMesh(pid2, &jcomm, &group1, &compid1); // receive from component 1
     CHECKRC(ierr, "cannot receive elements")
     std::string wopts;
     wopts   = "PARALLEL=WRITE_PART;";
     ierr = iMOAB_WriteMesh(pid2, (char*)outfile , (char*)wopts.c_str(), strlen(outfile), strlen(wopts.c_str()) );
     CHECKRC(ierr, "cannot write received mesh" )
  }

  MPI_Barrier(jcomm);

  // we can now free the sender buffers
  if (comm1 != MPI_COMM_NULL)
     ierr = iMOAB_FreeSenderBuffers(pid1, &jcomm, &compid2);

  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(pid2);
    CHECKRC(ierr, "cannot deregister app 2 receiver" )
  }

  if (comm1 != MPI_COMM_NULL) {
     ierr = iMOAB_DeregisterApplication(pid1);
     CHECKRC(ierr, "cannot deregister app 1 sender" )
  }


  ierr = iMOAB_Finalize();
  CHECKRC(ierr, "did not finalize iMOAB" )


  if (MPI_COMM_NULL != comm1) MPI_Comm_free(&comm1);
  if (MPI_COMM_NULL != comm2) MPI_Comm_free(&comm2);

  MPI_Group_free(&group1);
  MPI_Group_free(&group2);
  return MB_SUCCESS;
}
// migrate from 2 tasks to 3 tasks
ErrorCode migrate_graph( const char* filename )
{
  return migrate_smart(filename, "migrate_graph.h5m", 1);
}

ErrorCode migrate_geom( const char* filename )
{
  return migrate_smart(filename, "migrate_geom.h5m", 2);
}

int main( int argc, char* argv[] )
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );


  MPI_Comm_dup(MPI_COMM_WORLD, &jcomm);
  MPI_Comm_group(jcomm, &jgroup);

  ProgOptions opts;

  //std::string inputfile, outfile("out.h5m"), netcdfFile, variable_name, sefile_name;
  std::string filename;
  filename = TestDir + "/field1.h5m";
  startG1=0;
  startG2=0;
  endG1 = 0;
  endG2 = 1;
  opts.addOpt<std::string>("file,f","source file", &filename);

  opts.addOpt<int>("startSender,a", "start task for source layout", &startG1);
  opts.addOpt<int>("endSender,b", "end task for source layout", &endG1);
  opts.addOpt<int>("startRecv,c", "start task for receiver layout", &startG2);
  opts.addOpt<int>("endRecv,d", "end task for receiver layout", &endG2);

  opts.parseCommandLine(argc, argv);

  if (rank == 0)
  {
    std::cout << " input file : " << filename << "\n";
    std::cout << " sender   on tasks: " << startG1 << ":" << endG1 <<  "\n";
    std::cout << " receiver on tasks: " << startG2 << ":" << endG2 <<  "\n";
  }

  int num_errors = 0;


  num_errors += RUN_TEST_ARG2( migrate_graph, filename.c_str() );
  num_errors += RUN_TEST_ARG2( migrate_geom, filename.c_str() );

  if (rank == 0) {
    if (!num_errors)
      std::cout << "All tests passed" << std::endl;
    else
      std::cout << num_errors << " TESTS FAILED!" << std::endl;
  }

  MPI_Group_free(&jgroup);
  MPI_Comm_free(&jcomm);
  MPI_Finalize();
  return num_errors;
}





