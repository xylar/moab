/*
 * This imoab_coupler test will simulate coupling between 2 components
 * 2 meshes will be loaded from 2 files (atm and ocean), and they will be migrated to
 * 2 processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 *
 */

#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include <iostream>
#include <sstream>

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s\n", message); return 1;}

using namespace moab;

int main( int argc, char* argv[] )
{
  int rank, size, ierr;
  MPI_Comm jcomm; // will be a copy of the global
  MPI_Group jgroup;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  if (size!=2)
    return 1; // launch it on 2 processors only

  MPI_Comm_dup(MPI_COMM_WORLD, &jcomm);
  MPI_Comm_group(jcomm, &jgroup);

  int compid1=4, compid2=7, compid3=11, compid4=12, intxid=15;  // component ids are unique over all pes, and established in advance;
  int nghlay=0; // number of ghost layers for loading the file
  int groupTasks[2]; // at most 2 tasks
  int startG1=0, startG2=0, endG1=0, endG2=1;

  // load atm on 1 proc, ocean on 2, migrate both to 2 procs, then compute intx
  // later, we need to compute weight matrix with tempestremap
  std::string filename1, filename2;
  filename1 = TestDir + "/atm.h5m";
  filename2 = TestDir + "/mpas.h5m";

  // load files on 2 different communicators, groups
  // first groups has task 0, second group tasks 0 and 1
  // coupler will be on joint tasks, will be on a third group (0 and 1, again)

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

  int appID1;
  iMOAB_AppID pid1=&appID1; // atm
  int appID2;
  iMOAB_AppID pid2=&appID2; // ocn
  int appID3, appID4, appID5;
  iMOAB_AppID pid3=&appID3; // atm on coupler PEs
  iMOAB_AppID pid4=&appID4; // ocn on coupler PEs
  iMOAB_AppID pid5= &appID5; // intx atm -ocn on coupler PEs

  ierr = iMOAB_RegisterApplication("ATMX", &jcomm, &compid3, pid3);
  CHECKRC(ierr, "can't register atm over coupler pes ")

  ierr = iMOAB_RegisterApplication("OCNX", &jcomm, &compid4, pid4);
  CHECKRC(ierr, "can't register ocn over coupler pes ")

  std::string   readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");

  if (comm1 != MPI_COMM_NULL) {
    ierr = iMOAB_RegisterApplication("ATM1", &comm1, &compid1, pid1);
    CHECKRC(ierr, "can't register app1 ")
    // load first mesh
    ierr = iMOAB_LoadMesh(pid1, filename1.c_str(), readopts.c_str(), &nghlay, filename1.length(), strlen(readopts.c_str()) );
    CHECKRC(ierr, "can't load mesh1 ")
    // then send mesh to coupler pes, on pid3
    ierr = iMOAB_SendMesh(pid1, &jcomm, &jgroup, &compid3); // send to component 3, on coupler pes
    CHECKRC(ierr, "cannot send elements" )
  }
  // now, receive meshes, on joint communicator; first mesh 1
  ierr = iMOAB_ReceiveMesh(pid3, &jcomm, &group1, &compid1); // receive from component 1
  CHECKRC(ierr, "cannot receive elements on ATMX app")

  // we can now free the sender buffers
  if (comm1 != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(pid1, &jcomm, &compid3);
    CHECKRC(ierr, "cannot free buffers used to send atm mesh")
  }

  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_RegisterApplication("OCN1", &comm2, &compid2, pid2);
    CHECKRC(ierr, "can't register app2 ")
    // load second mesh
    ierr = iMOAB_LoadMesh(pid2, filename2.c_str(), readopts.c_str(), &nghlay, filename2.length(), strlen(readopts.c_str()) );
    CHECKRC(ierr, "can't load mesh2, on pid2 ")
    // then send mesh to coupler pes, on pid3
    ierr = iMOAB_SendMesh(pid2, &jcomm, &jgroup, &compid4); // send to component 3, on coupler pes
    CHECKRC(ierr, "cannot send elements" )
  }


  ierr = iMOAB_ReceiveMesh(pid4, &jcomm, &group2, &compid2); // receive from component 2
  CHECKRC(ierr, "cannot receive elements on OCNX app")

  if (comm2 != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(pid2, &jcomm, &compid4);
    CHECKRC(ierr, "cannot free buffers used to send ocn mesh")
  }

#ifdef MOAB_HAVE_TEMPESTREMAP
  // now compute intersection between OCNx and ATMx on coupler PEs
  ierr = iMOAB_RegisterApplication("AO", &jcomm, &intxid, pid5);
  CHECKRC(ierr, "can't register ocn_atm intx over coupler pes ")

  ierr = iMOAB_ComputeMeshIntersectionOnSphere(pid3, pid4, pid5);
  // check if intx valid, write some h5m intx file
  CHECKRC(ierr, "cannot compute intersection" )

  std::stringstream outf;
  outf<<"intx_0" << rank<<".h5m";
  std::string intxfile=outf.str(); // write in serial the intx file, for debugging
  char writeOptions[] ="";
  ierr = iMOAB_WriteMesh(pid5, (char*)intxfile.c_str(), writeOptions, (int)intxfile.length(), strlen(writeOptions));
  CHECKRC(ierr, "cannot write intx file result" )
  // the new graph will be for sending data from atm comp to coverage mesh;
  // it involves initial atm app; pid1; also migrate atm mesh on coupler pes, pid3
  // results are in pid5, intx mesh; remapper also has some info about coverage mesh
  ierr = iMOAB_CoverageGraph(&jcomm, pid1,  &compid1, pid3,  &compid3, pid5); // it happens over joint communicator
  CHECKRC(ierr, "cannot recompute direct coverage graph" )
  
  int disc_orders[2] = {4, 1};
  const char* disc_methods[2] = {"CGLL", "fv"};
  const char* dof_tag_names[2] = {"GLOBAL_ID", "GLOBAL_ID"};
  int fVolumetric=0, fValidate=1, fNoConserve=0;
  ierr = iMOAB_ComputeScalarProjectionWeights ( pid5,
                                                disc_methods[0], &disc_orders[0],
                                                disc_methods[1], &disc_orders[1],
                                                &fVolumetric, &fNoConserve, &fValidate,
                                                dof_tag_names[0], dof_tag_names[1],
                                                strlen(disc_methods[0]), strlen(disc_methods[1]),
                                                strlen(dof_tag_names[0]), strlen(dof_tag_names[1]) );
  CHECKRC(ierr, "cannot compute scalar projection weights" )

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


