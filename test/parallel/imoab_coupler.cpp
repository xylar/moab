/*
 * This imoab_coupler test will simulate coupling between 2 components
 * 2 meshes will be loaded from 2 files (atm and ocean), and they will be migrated to
 * 2 processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 *
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
    #error mbtempest tool requires MPI configuration
#endif

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include <iostream>
#include <sstream>

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s\n", message); return 1;}
#define PUSH_TIMER(operation)  { timer_ops = timer.time_since_birth(); opName = operation;}
#define POP_TIMER() { \
  double locElapsed=timer.time_since_birth() - timer_ops, minElapsed=0, maxElapsed=0; \
  MPI_Reduce(&locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); \
  MPI_Reduce(&locElapsed, &minElapsed, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); \
  if (!rank) std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed << ", avg = " << (maxElapsed+minElapsed)/2 << "\n"; \
  opName.clear(); \
}

using namespace moab;

// #define VERBOSE

int main( int argc, char* argv[] )
{
  int rank, size, ierr;
  MPI_Comm jcomm; // will be a copy of the global
  MPI_Group jgroup;
  std::string readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");

  // Timer data
  moab::CpuTimer timer;
  double timer_ops;
  std::string opName;

  int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
  repartitioner_scheme = 1; // use the graph partitioner in that caseS
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  MPI_Comm_dup(MPI_COMM_WORLD, &jcomm);
  MPI_Comm_group(jcomm, &jgroup);

  // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ; intx atm/ocn is not in e3sm yet, give a number
  //   6 * 100+ 18 = 618
  int compid1=5, compid2=17, compid3=6, compid4=18, intxid=618;  // component ids are unique over all pes, and established in advance;
  // compid1 is for atm on atm pes
  // compid2 is for ocean, on ocean pe
  // compid3 is for atm on coupler pes
  // compid4 is for ocean on coupelr pes
  // intxid is for intx atm / ocn on coupler pes
  int nghlay=0; // number of ghost layers for loading the file
  std::vector<int> groupTasks; // at most 2 tasks
  int startG1=0, startG2=0, endG1=size/2-1, endG2=size-1; // Support launch of imoab_coupler test on any combo of 2*x processes
  /* COMBOS THAT WORK */
  // int startG1=0, startG2=0, endG1=0, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=1, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=size-1, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  
  /* COMBOS THAT **DO NOT** WORK */
  // int startG1=0, startG2=0, endG1=0, endG2=1; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=1, endG2=1; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=size-1, endG2=size-1; // Support launch of imoab_coupler test on any combo of 2*x processes

  // load atm on 2 proc, ocean on 2, migrate both to 2 procs, then compute intx
  // later, we need to compute weight matrix with tempestremap
  std::string filename1, filename2;
#ifdef MOAB_HAVE_HDF5
  filename1 = TestDir + "/wholeATM_T.h5m";
  filename2 = TestDir + "/recMeshOcn.h5m";
#endif

  if (argc>2) {
    filename1 = std::string(argv[1]);
    filename2 = std::string(argv[2]);
    if (argc>3) repartitioner_scheme = atoi(argv[3]);
    if (argc>4)
    {
      int collocation;
      collocation = atoi(argv[4]); // 0) original (half atm) 1) (half ocn) 2) collocated
      if (1==collocation)
      {
        endG2=size/2-1; endG1=size-1;
      }
      else if (2==collocation)
      {
        endG2=size-1; endG1=size-1;
      }
    }
  }

  if (!rank)
  {
    std::cout << " atm file: " << filename1 << "\n   on tasks : " << startG1 << ":"<<endG1 <<
        "\n ocn file: " << filename2 << "\n     on tasks : " << startG2 << ":"<<endG2 <<
        "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
  }

  // load files on 2 different communicators, groups
  // first groups has task 0, second group tasks 0 and 1
  // coupler will be on joint tasks, will be on a third group (0 and 1, again)

  groupTasks.resize(size, 0);

  MPI_Group group1, group2;
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, &groupTasks[0], &group1);
  CHECKRC(ierr, "can't create group1")

  groupTasks.clear();
  groupTasks.resize(size, 0);
  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, &groupTasks[0], &group2);
  CHECKRC(ierr, "can't create group2")

  // create 2 communicators, one for each group
  int tagcomm1 = 1, tagcomm2 = 2;
  int rankInComm1 = -1, rankInComm2 = -1;
  MPI_Comm comm1, comm2;
  // comm1 is for atmosphere app;
  ierr = MPI_Comm_create_group(jcomm, group1, tagcomm1, &comm1);
  CHECKRC(ierr, "can't create comm1")

  // comm2 is for ocean app
  ierr = MPI_Comm_create_group(jcomm, group2, tagcomm2, &comm2);
  CHECKRC(ierr, "can't create comm2")

  ierr = iMOAB_Initialize(0, 0); // not really needed anything from argc, argv, yet; maybe we should
  CHECKRC(ierr, "can't initialize iMOAB")

  int appID1 =-1;
  iMOAB_AppID pid1=&appID1; // atm
  int appID2 =-1;
  iMOAB_AppID pid2=&appID2; // ocn
  int appID3=-1, appID4=-1, appID5=-1;// -1 means it is not initialized
  iMOAB_AppID pid3=&appID3; // atm on coupler PEs
  iMOAB_AppID pid4=&appID4; // ocn on coupler PEs
  iMOAB_AppID pid5= &appID5; // intx atm -ocn on coupler PEs

  ierr = iMOAB_RegisterApplication("ATMX", &jcomm, &compid3, pid3); // atm on coupler pes
  CHECKRC(ierr, "can't register atm over coupler pes ")

  ierr = iMOAB_RegisterApplication("OCNX", &jcomm, &compid4, pid4); // ocn on coupler pes
  CHECKRC(ierr, "can't register ocn over coupler pes ")

  PUSH_TIMER("Load source mesh")
  if (comm1 != MPI_COMM_NULL) {
    MPI_Comm_rank( comm1, &rankInComm1 );
    double t1= MPI_Wtime();
    ierr = iMOAB_RegisterApplication("ATM1", &comm1, &compid1, pid1);
    CHECKRC(ierr, "can't register app1 ")
    // load first mesh
    ierr = iMOAB_LoadMesh(pid1, filename1.c_str(), readopts.c_str(), &nghlay, filename1.length(), strlen(readopts.c_str()) );
    CHECKRC(ierr, "can't load mesh1 ")
    double t2= MPI_Wtime();
    if (!rankInComm1) std::cout << "[LOG] load atm mesh:" << t2-t1 << "\n";
    // then send mesh to coupler pes, on pid3
    ierr = iMOAB_SendMesh(pid1, &jcomm, &jgroup, &compid3, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKRC(ierr, "cannot send elements" )
    double t3= MPI_Wtime();
    if (!rankInComm1) std::cout << "[LOG] compute partition and send atm mesh: " << t3-t2 << "\n";
  }
  // now, receive meshes, on joint communicator; first mesh 1
  double t4=MPI_Wtime();
  ierr = iMOAB_ReceiveMesh(pid3, &jcomm, &group1, &compid1); // receive from component 1
  CHECKRC(ierr, "cannot receive elements on ATMX app")
  double t5=MPI_Wtime();
  if (!rank) std::cout << "[LOG] receive atm mesh and resolve shared entities: " << t5-t4 << "\n";
  POP_TIMER()

  // we can now free the sender buffers
  if (comm1 != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(pid1, &jcomm, &compid3);
    CHECKRC(ierr, "cannot free buffers used to send atm mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  PUSH_TIMER("Load target mesh")
  if (comm2 != MPI_COMM_NULL) {
    MPI_Comm_rank( comm2, &rankInComm2 );
    ierr = iMOAB_RegisterApplication("OCN1", &comm2, &compid2, pid2);
    CHECKRC(ierr, "can't register app2 ")
    // load second mesh
    double t6 = MPI_Wtime();
    ierr = iMOAB_LoadMesh(pid2, filename2.c_str(), readopts.c_str(), &nghlay, filename2.length(), strlen(readopts.c_str()) );
    CHECKRC(ierr, "can't load mesh2, on pid2 ")
    double t7 = MPI_Wtime();
    if (!rankInComm2) std::cout << "[LOG] load ocn mesh:"<< t7-t6 << "\n";
    // then send mesh to coupler pes, on pid3
    ierr = iMOAB_SendMesh(pid2, &jcomm, &jgroup, &compid4, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKRC(ierr, "cannot send elements" )
    double t8 = MPI_Wtime();
    if (!rankInComm2) std::cout << "[LOG] compute partition and send ocn mesh:"<< t8-t7 << "\n";
  }

  double t9=MPI_Wtime();
  ierr = iMOAB_ReceiveMesh(pid4, &jcomm, &group2, &compid2); // receive from component 2, ocn, on coupler pes
  CHECKRC(ierr, "cannot receive elements on OCNX app")
  double t10=MPI_Wtime();
  if (!rank) std::cout << "[LOG] receive ocn mesh and resolve shared entities: "<< t10-t9 << "\n";
  POP_TIMER()

  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(pid2, &jcomm, &compid4);
    CHECKRC(ierr, "cannot free buffers used to send ocn mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);
  char outputFileTgt3[] = "recvTarget.h5m";
  #ifdef MOAB_HAVE_MPI
    char writeOptions3[] ="PARALLEL=WRITE_PART";
  #else
    char writeOptions3[] ="";
  #endif
  ierr = iMOAB_WriteMesh(pid4, outputFileTgt3, writeOptions3,
    strlen(outputFileTgt3), strlen(writeOptions3) );
  CHECKRC(ierr, "cannot write ocn mesh after receiving")

  double t11=MPI_Wtime();
    if (!rank) std::cout << "[LOG] write received ocn mesh on coupler pes: "<< t11-t10 << "\n";

#ifdef MOAB_HAVE_TEMPESTREMAP
  // now compute intersection between OCNx and ATMx on coupler PEs
  ierr = iMOAB_RegisterApplication("AO", &jcomm, &intxid, pid5);
  CHECKRC(ierr, "can't register ocn_atm intx over coupler pes ")

  PUSH_TIMER("Compute source-target mesh intersection")
  ierr = iMOAB_ComputeMeshIntersectionOnSphere(pid3, pid4, pid5); // coverage mesh was computed here, for pid3, atm on coupler pes
  // basically, atm was redistributed according to target (ocean) partition, to "cover" the ocean partitions
  // check if intx valid, write some h5m intx file
  CHECKRC(ierr, "cannot compute intersection" )
  POP_TIMER()

#ifdef VERBOSE
  std::stringstream outf;
  outf<<"intx_0" << rank<<".h5m";
  std::string intxfile=outf.str(); // write in serial the intx file, for debugging
  char writeOptions[] ="";
  ierr = iMOAB_WriteMesh(pid5, (char*)intxfile.c_str(), writeOptions, (int)intxfile.length(), strlen(writeOptions));
  CHECKRC(ierr, "cannot write intx file result" )
#endif

  // the new graph will be for sending data from atm comp to coverage mesh;
  // it involves initial atm app; pid1; also migrate atm mesh on coupler pes, pid3
  // results are in pid5, intx mesh; remapper also has some info about coverage mesh
  // afteer this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more precise info about
  // what to send ; every time, we will send also the element global id, which should uniquely identify the element
  PUSH_TIMER("Compute coverage graph")
  ierr = iMOAB_CoverageGraph(&jcomm, pid1,  &compid1, pid3,  &compid3, pid5); // it happens over joint communicator
  CHECKRC(ierr, "cannot recompute direct coverage graph" )
  POP_TIMER()

  const char* weights_identifiers[1] = {"scalar"};
  int disc_orders[2] = {4, 1};
  const char* disc_methods[2] = {"cgll", "fv"};
  const char* dof_tag_names[2] = {"GLOBAL_DOFS", "GLOBAL_ID"};
  int fVolumetric=0, fValidate=1, fNoConserve=0;
  PUSH_TIMER("Compute the projection weights with TempestRemap")
  ierr = iMOAB_ComputeScalarProjectionWeights ( pid5, weights_identifiers[0],
                                                disc_methods[0], &disc_orders[0],
                                                disc_methods[1], &disc_orders[1],
                                                &fVolumetric, &fNoConserve, &fValidate,
                                                dof_tag_names[0], dof_tag_names[1],
                                                strlen(disc_methods[0]), strlen(disc_methods[1]),
                                                strlen(dof_tag_names[0]), strlen(dof_tag_names[1]) );
  CHECKRC(ierr, "cannot compute scalar projection weights" )
  POP_TIMER()

  const char* fieldname = "a2oTAG";
  const char* fieldnameT = "a2oTAG_proj";
  int tagIndex[2];
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int num_components1 = disc_orders[0]*disc_orders[0], num_components2 = disc_orders[1]*disc_orders[1];

  ierr = iMOAB_DefineTagStorage(pid3, fieldname, &tagTypes[0], &num_components1, &tagIndex[0],  strlen(fieldname) );
  CHECKRC(ierr, "failed to define the field tag");

  ierr = iMOAB_DefineTagStorage(pid4, fieldnameT, &tagTypes[1], &num_components2, &tagIndex[1],  strlen(fieldnameT) );
  CHECKRC(ierr, "failed to define the field tag");

  // need to make sure that the coverage mesh (created during intx method) received the tag that need to be projected to target
  // so far, the coverage mesh has only the ids and global dofs;
  // need to change the migrate method to accommodate any GLL tag
  // now send a tag from original atmosphere (pid1) towards migrated coverage mesh (pid3), using the new coverage graph communicator

  // make the tag 0, to check we are actually sending needed data
  {
    if (appID3 >= 0)
    {
      int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
        /*
         * Each process in the communicator will have access to a local mesh instance, which will contain the
         * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible blocks,
         * number of sidesets and nodesets boundary conditions will be returned in size 3 arrays, for local, ghost and total
         * numbers.
         */
        ierr = iMOAB_GetMeshInfo(  pid3, nverts, nelem, nblocks, nsbc, ndbc);
        CHECKRC(ierr, "failed to get num primary elems");
        int numAllElem = nelem[2];
        std::vector<double> vals;
        int storLeng = num_components1*numAllElem;
        vals.resize(storLeng);
        for (int k=0; k<storLeng; k++)
          vals[k] = 0.;
        int eetype = 1;
        ierr = iMOAB_SetDoubleTagStorage ( pid3, "a2oTAG", &storLeng, &eetype, &vals[0], strlen("a2oTAG"));
        CHECKRC(ierr, "cannot make tag nul")
        // set the tag to 0
    }
  }
  
  PUSH_TIMER("Send/receive data from component to coupler")
  if (comm1 != MPI_COMM_NULL ){

    // basically, adjust the migration of the tag we want to project; it was sent initially with
    // trivial partitioning, now we need to adjust it for "coverage" mesh
     // as always, use nonblocking sends
     ierr = iMOAB_SendElementTag(pid1, &compid1, &compid3, "a2oTAG", &jcomm, strlen("a2oTAG"));
     CHECKRC(ierr, "cannot send tag values")
  }
  // receive on atm on coupler pes, that was redistributed according to coverage
  ierr = iMOAB_ReceiveElementTag(pid3, &compid3, &compid1, "a2oTAG", &jcomm, strlen("a2oTAG"));
  CHECKRC(ierr, "cannot receive tag values")
  POP_TIMER()

#ifdef VERBOSE
    char outputFileRecvd[] = "recvAtmCoup.h5m";
    ierr = iMOAB_WriteMesh(pid3, outputFileRecvd, writeOptions3,
        strlen(outputFileRecvd), strlen(writeOptions3) );
#endif
    // we can now free the sender buffers
     if (comm1 != MPI_COMM_NULL) {
       ierr = iMOAB_FreeSenderBuffers(pid1, &jcomm, &compid3);
       CHECKRC(ierr, "cannot free buffers used to resend atm mesh tag towards the coverage mesh")
     }

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined 
     on the source mesh and get the projection on the target mesh */
  PUSH_TIMER("Apply Scalar projection weights")
  ierr = iMOAB_ApplyScalarProjectionWeights ( pid5, weights_identifiers[0],
                                            fieldname,
                                            fieldnameT,
                                            strlen(fieldname),
                                            strlen(fieldnameT)
                                            );
  CHECKRC(ierr, "failed to compute projection weight application");
  POP_TIMER()

  char outputFileTgt[] = "fIntxTarget.h5m";
#ifdef MOAB_HAVE_MPI
    char writeOptions2[] ="PARALLEL=WRITE_PART";
#else
    char writeOptions2[] ="";
#endif
  ierr = iMOAB_WriteMesh(pid4, outputFileTgt, writeOptions2,
    strlen(outputFileTgt), strlen(writeOptions2) );

  ierr = iMOAB_DeregisterApplication(pid5);
  CHECKRC(ierr, "cannot deregister app intx AO" )

#endif

  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(pid2);
    CHECKRC(ierr, "cannot deregister app OCN1" )
  }
  if (comm1 != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(pid1);
    CHECKRC(ierr, "cannot deregister app ATM1" )
  }

  ierr = iMOAB_DeregisterApplication(pid4);
  CHECKRC(ierr, "cannot deregister app OCNX" )

  ierr = iMOAB_DeregisterApplication(pid3);
  CHECKRC(ierr, "cannot deregister app OCNX" )

  ierr = iMOAB_Finalize();
  CHECKRC(ierr, "did not finalize iMOAB" )

  if (MPI_COMM_NULL != comm1) MPI_Comm_free(&comm1);
  if (MPI_COMM_NULL != comm2) MPI_Comm_free(&comm2);

  MPI_Group_free(&group1);
  MPI_Group_free(&group2);
  MPI_Group_free(&jgroup);
  MPI_Comm_free(&jcomm);

  return 0;
}


