/*
 * migrate_test will contain tests for migrating meshes in parallel environments, with iMOAB api
 * these  methods are tested also in the example MigrateMesh.F90, with variable
 * numbers of processes; migrate_test is launched usually on 2 processes, and it tests
 * various cases
 * a mesh is read on senders tasks, sent to receivers tasks, and then written out for verification
 * It depends on hdf5 parallel for reading and writing in parallel
 */

#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"

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

ErrorCode migrate_1_1( const char* filename );
ErrorCode migrate_1_2( const char* filename );
ErrorCode migrate_2_1( const char* filename );
ErrorCode migrate_2_2( const char* filename );

// some global variables, used by all tests
int rank, size, ierr;

int compid1, compid2;  // component ids are unique over all pes, and established in advance;
int nghlay; // number of ghost layers for loading the file
int groupTasks[2]; // at most 2 tasks
int startG1, startG2, endG1, endG2;

MPI_Comm jcomm; // will be a copy of the global
MPI_Group jgroup;

int main( int argc, char* argv[] )
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  if (size!=2)
  {
    std::cout <<"run on 2 tasks only.\n";
    MPI_Finalize();
    return 1;
  }

  MPI_Comm_dup(MPI_COMM_WORLD, &jcomm);
  MPI_Comm_group(jcomm, &jgroup);

	std::string filename;
	filename = TestDir + "/Homme_2pt.h5m";
  if (argc>1)
  {
    filename = argv[1];
  }
  int num_errors = 0;
  num_errors += RUN_TEST_ARG2( migrate_1_1, filename.c_str() );
  num_errors += RUN_TEST_ARG2( migrate_1_2, filename.c_str() );
  num_errors += RUN_TEST_ARG2( migrate_2_1, filename.c_str() );
  num_errors += RUN_TEST_ARG2( migrate_2_2, filename.c_str() );

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

ErrorCode migrate(const char*filename)
{
  // first create MPI groups

  std::string filen(filename);
  MPI_Group group1, group2;
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, groupTasks, &group1);
  CHECKRC(ierr, "can't create group1")

  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, groupTasks, &group2);
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
      ierr = iMOAB_SendMesh(pid1, &jcomm, &group2, &compid2); // send to component 2
      CHECKRC(ierr, "cannot send elements" )
  }

  if (comm2 != MPI_COMM_NULL) {
     ierr = iMOAB_ReceiveMesh(pid2, &jcomm, &group1, &compid1); // receive from component 1
     CHECKRC(ierr, "cannot receive elements")
     std::string wopts;
     if (appID2 == 0)
       wopts   = "PARALLEL=WRITE_PART;PARALLEL_COMM=0";
     else if (appID2 == 1)
       wopts   = "PARALLEL=WRITE_PART;PARALLEL_COMM=1";
     else
       CHECKRC(1, "wrong pid2" )

     char wfile[] = "recvMesh.h5m";
     ierr = iMOAB_WriteMesh(pid2, wfile , (char*)wopts.c_str(), strlen(wfile), strlen(wopts.c_str()) );
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
// migrate from task 0 to task 1, non overlapping
ErrorCode migrate_1_1( const char* filename )
{
  startG1= endG1 = 0;
  startG2 = endG2 = 1;
  return migrate(filename);
}
// migrate from task 0 to 2 tasks (0 and 1)
ErrorCode migrate_1_2( const char* filename )
{
  startG1= endG1 = startG2 = 0;
  endG2 = 1;
  return migrate(filename);
}

// migrate from 2 tasks (0, 1) to 1 task (0)
ErrorCode migrate_2_1( const char* filename )
{
  startG1= endG2 = startG2 = 0;
  endG1 = 1;
  return migrate(filename);
}

// migrate from 2 tasks to 2 tasks (overkill)
ErrorCode migrate_2_2( const char* filename )
{
  startG1 = startG2 = 0;
  endG1 = endG2 = 1;
  return migrate(filename);
}

