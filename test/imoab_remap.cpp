#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_MPI
#  include "moab_mpi.h"
#endif

#include "moab/iMOAB.h"

// for malloc, free:
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#ifdef MESHDIR
const char* TestDir = STRINGIFY(MESHDIR);
#else
#error Specify MESHDIR to compile test
#endif

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s", message); return 1;}

// this test will be run in serial only
int main(int argc, char * argv[])
{
  int nprocs=1, rank=0;
  char filen1[200], filen2[200];
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
#endif

#ifdef MOAB_HAVE_HDF5
  strcpy(filen1, TestDir);
  strcpy(filen2, TestDir);
  strcat(filen1, "/wholeATM_T.h5m");
  strcat(filen2, "/recMeshOcn.h5m");
#endif

  if (argc>2) {
    strcpy(filen1, argv[1]);
    strcpy(filen2, argv[2]);
  }

  if (!strlen(filen1) || !strlen(filen2)) {
    printf("Please provide file names as arguments\n");
    return 1;
  }

  /*
   * MOAB needs to be initialized; A MOAB instance will be created, and will be used by each application
   * in this framework. There is no parallel context yet.
   */
  ErrCode rc = iMOAB_Initialize(argc, argv);
  CHECKRC(rc, "failed to initialize MOAB");

  int appID1, appID2, appID3;
  iMOAB_AppID pid1=&appID1;
  iMOAB_AppID pid2=&appID2;
  iMOAB_AppID pid3=&appID3;
  /*
   * Each application has to be registered once. A mesh set and a parallel communicator will be associated
   * with each application. A unique application id will be returned, and will be used for all future
   * mesh operations/queries.
   */
  int compid1 = 10, compid2 = 20, compid3 = 100;
  rc = iMOAB_RegisterApplication( "COUP_APP1",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &compid1,
      pid1);
  CHECKRC(rc, "failed to register application1");

  rc = iMOAB_RegisterApplication( "COUP_APP2",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &compid2,
      pid2);
  CHECKRC(rc, "failed to register application2");

  rc = iMOAB_RegisterApplication( "COUPLER",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &compid3,
      pid3);
  CHECKRC(rc, "failed to register application2");

  const char *read_opts="";
  int num_ghost_layers=0;

  /*
   * Loading the mesh is a parallel IO operation. Ghost layers can be exchanged too, and default MOAB
   * sets are augmented with ghost elements. By convention, blocks correspond to MATERIAL_SET sets,
   * side sets with NEUMANN_SET sets, node sets with DIRICHLET_SET sets. Each element and vertex entity should have
   * a GLOBAL ID tag in the file, which will be available for visible entities
   */
  rc = iMOAB_LoadMesh(  pid1, filen1, read_opts, &num_ghost_layers, strlen(filen1), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");
  rc = iMOAB_LoadMesh(  pid2, filen2, read_opts, &num_ghost_layers, strlen(filen2), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");

  int nverts[3], nelem[3];
  /*
   * Each process in the communicator will have access to a local mesh instance, which will contain the
   * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible blocks,
   * number of sidesets and nodesets boundary conditions will be returned in size 3 arrays, for local, ghost and total
   * numbers.
   */
  rc = iMOAB_GetMeshInfo( pid1, nverts, nelem, 0, 0, 0);
  CHECKRC(rc, "failed to get mesh info");
  printf("Source Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  rc = iMOAB_GetMeshInfo( pid2, nverts, nelem, 0, 0, 0);
  CHECKRC(rc, "failed to get mesh info");
  printf("Destination Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  /*
   * The 2 tags used in this example exist in the file, already.
   * If a user needs a new tag, it can be defined with the same call as this one
   * this method, iMOAB_DefineTagStorage, will return a local index for the tag.
   * The name of the tag is case sensitive.
   * This method is collective.
   */
  // int disc_orders[2] = {1, 1};
  // const char* disc_methods[2] = {"fv", "fv"};
  // const char* dof_tag_names[2] = {"GLOBAL_ID", "GLOBAL_ID"};
  int disc_orders[2] = {4, 1};
  const char* disc_methods[2] = {"cgll", "fv"};
  const char* dof_tag_names[2] = {"GLOBAL_DOFS", "GLOBAL_ID"};
  
  const char* weights_identifiers[2] = {"scalar", "scalar_conservative"};
  int fMonotoneTypeID=0, fVolumetric=0, fValidate=1, fNoConserve=0;

  const char* fieldname = "a2oTbot";
  const char* fieldnameT = "a2oTbot_proj";
  const char* fieldnameTnc = "a2oTbot_projnocons";
  int tagIndex[2];
  // int entTypes[2] = {1, 1}; /* both on elements; */
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int num_components1 = disc_orders[0]*disc_orders[0], num_components2 = disc_orders[1]*disc_orders[1];

  rc = iMOAB_DefineTagStorage(pid1, fieldname, &tagTypes[0], &num_components1, &tagIndex[0],  strlen(fieldname) );
  CHECKRC(rc, "failed to define the field tag");

  rc = iMOAB_DefineTagStorage(pid2, fieldnameT, &tagTypes[1], &num_components2, &tagIndex[1],  strlen(fieldnameT) );
  CHECKRC(rc, "failed to define the field tag");

  rc = iMOAB_DefineTagStorage(pid2, fieldnameTnc, &tagTypes[1], &num_components2, &tagIndex[1],  strlen(fieldnameTnc) );
  CHECKRC(rc, "failed to define the field tag");

  /* Next compute the mesh intersection on the sphere between the source and target meshes */
  rc = iMOAB_ComputeMeshIntersectionOnSphere(pid1, pid2, pid3);
  CHECKRC(rc, "failed to compute mesh intersection");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=1;
  rc = iMOAB_ComputeScalarProjectionWeights ( pid3,
                                              weights_identifiers[0],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[0]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKRC(rc, "failed to compute remapping projection weights for scalar non-conservative field");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=0;
  rc = iMOAB_ComputeScalarProjectionWeights ( pid3,
                                              weights_identifiers[1],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[1]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKRC(rc, "failed to compute remapping projection weights for scalar conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  rc = iMOAB_ApplyScalarProjectionWeights ( pid3,
                                            weights_identifiers[0],
                                            fieldname,
                                            fieldnameTnc,
                                            strlen(weights_identifiers[0]),
                                            strlen(fieldname),
                                            strlen(fieldnameTnc)
                                            );
  CHECKRC(rc, "failed to compute projection weight application for scalar non-conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  rc = iMOAB_ApplyScalarProjectionWeights ( pid3,
                                            weights_identifiers[1],
                                            fieldname,
                                            fieldnameT,
                                            strlen(weights_identifiers[1]),
                                            strlen(fieldname),
                                            strlen(fieldnameT)
                                            );
  CHECKRC(rc, "failed to compute projection weight application for scalar conservative field");

  /*
   * the file can be written in parallel, and it will contain additional tags defined by the user
   * we may extend the method to write only desired tags to the file
   */
  {
    // free allocated data
    char outputFileOv[]  = "fIntxOverlap.h5m";
    char outputFileTgt[] = "fIntxTarget.h5m";
    char writeOptions[] ="";

    rc = iMOAB_WriteMesh(pid2, outputFileTgt, writeOptions,
      strlen(outputFileTgt), strlen(writeOptions) );

    rc = iMOAB_WriteMesh(pid3, outputFileOv, writeOptions,
      strlen(outputFileOv), strlen(writeOptions) );
  }

  /*
   * deregistering application will delete all mesh entities associated with the application and will
   *  free allocated tag storage.
   */
  rc = iMOAB_DeregisterApplication(pid3);
  CHECKRC(rc, "failed to de-register application3");
  rc = iMOAB_DeregisterApplication(pid2);
  CHECKRC(rc, "failed to de-register application2");
  rc = iMOAB_DeregisterApplication(pid1);
  CHECKRC(rc, "failed to de-register application1");

  /*
   * this method will delete MOAB instance
   */
  rc = iMOAB_Finalize();
  CHECKRC(rc, "failed to finalize MOAB");

  return 0;
}

