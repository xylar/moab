
#include "mpi.h"
#include  "../src/moab/imoab.h"
// for malloc, free:
#include <stdlib.h>
#include <string.h>

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s", message); return 1;}
int main(int argc, char * argv[])
{

  MPI_Init(&argc, &argv);
  int nprocs, rank;

  MPI_Comm comm=MPI_COMM_WORLD;

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  char * filen = "p8ex1.h5m";
  if (argc>1)
    filen = argv[1];

  ErrCode rc = iMOABInitialize(argc, argv);

  CHECKRC(rc, "failed to initialize MOAB");
  int num_global_vertices=0, num_global_elements=0, num_dimension=0, num_parts=0;
  rc = ReadHeaderInfo ( filen, &num_global_vertices, &num_global_elements, &num_dimension,
      &num_parts, (int)strlen(filen) );

  CHECKRC(rc, "failed to read header info");

  if (0==rank)
  {
    printf("file %s has %d vertices, %d elements, %d parts in partition\n", filen,
        num_global_vertices, num_global_elements, num_parts);
  }
  int appID;
  iMOAB_AppID pid=&appID;
  rc = RegisterApplication( "PROTEUS", &comm,  pid);
  CHECKRC(rc, "failed to register application");
  char *read_opts="PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
  int num_ghost_layers=1;
  rc = LoadMesh(  pid, filen, read_opts, &num_ghost_layers, strlen(filen), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");

  // nverts[0]: number of local vertices (used to represent initial partition, before ghosting
  //   in moab terminology, nverts[0] has owned and shared vertices on the interfaces
  //   some of the shared vertices could be actually owned by another task
  // nverts[1]: number of ghost vertices, appear after ghosting
  // nverts[2]: total number of vertices, local + ghost

  // we do not have that problem / terminology for elements
  // elements are either owned (local, in the initial partition) or ghosted

  // nblocks as the number of material sets( visible blocks) is different
  // blocks[0] represent number of blocks completely owned by current partition
  // blocks[1] is number of blocks with at least one ghost element in it
  // blocks[2] is the total number of blocks visible on current rank
  // any of the blocks [0] or blocks[1] can be 0
  int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
  rc = GetMeshInfo(  pid, nverts, nelem, nblocks, nsbc, ndbc);
  CHECKRC(rc, "failed to get mesh info");

  iMOAB_GlobalID * vGlobalID = (iMOAB_GlobalID*)malloc(nverts[2]*sizeof(iMOAB_GlobalID)) ;
  rc = GetVertexID(pid, &nverts[2], vGlobalID);
  CHECKRC(rc, "failed to get vertex id info");

  int * vranks = (int*)malloc(nverts[2]*sizeof(int));
  rc =GetVertexOwnership(pid, &nverts[2], vranks );
  CHECKRC(rc, "failed to get vertex ranks");

  double * coords = (double*) malloc(3*nverts[2]*sizeof(double));
  int size_coords= 3*nverts[2];
  rc = GetVisibleVerticesCoordinates( pid, &size_coords, coords);
  CHECKRC(rc, "failed to get coordinates");

  iMOAB_GlobalID * gbIDs = (iMOAB_GlobalID*) malloc(nblocks[2]*sizeof(iMOAB_GlobalID));
  rc = GetBlockID(pid, &nblocks[2], gbIDs);
  CHECKRC(rc, "failed to get block info");

  /* the 2 tags used in this example exist in the file, already  */
  int tagIndex[2];
  int entTypes[2] = {0, 1}; /* first is on vertex, second is on elements; */
  int tagTypes[2] = { DENSE_INTEGER, DENSE_DOUBLE } ;
  int num_components = 1;
  rc = DefineTagStorage(pid, "INTFIELD", &tagTypes[0], &num_components, &tagIndex[0],  strlen("INTFIELD") );
  CHECKRC(rc, "failed to get tag INTFIELD ");

  rc = DefineTagStorage(pid, "DFIELD", &tagTypes[1], &num_components, &tagIndex[1],  strlen("DFIELD") );
  CHECKRC(rc, "failed to get tag DFIELD ");

  // synchronize one of the tags only, just to see what happens
  int num_tags_to_sync=1;
  rc = SynchronizeTags(pid, &num_tags_to_sync, &tagIndex[0], &tagTypes[0] );
  CHECKRC(rc, "failed to sync tag INTFIELD ");

  for (int irank=0; irank<nprocs; irank++)
  {
    if (irank==rank)
    {
      /* printf some of the block info */
      printf("on rank %d, there are \n"
              "  %3d visible vertices of which  %3d local  %3d ghost \n"
              "  %3d visible elements of which  %3d owned  %3d ghost \n"
              "  %3d visible blocks\n"
              "  %3d visible neumann BCs\n"
              "  %3d visible dirichlet BCs\n", rank,
              nverts[2], nverts[0], nverts[1],
              nelem[2], nelem[0], nelem[1],
              nblocks[2], nsbc[2], ndbc[2]);
      /* printf some of the vertex id infos */
      int numToPrint = nverts[2];
      printf("on rank %d vertex info:\n", rank);
      for (int i=0; i<numToPrint; i++)
        printf(" vertex local id: %3d, rank ID:%d  global ID: %3d  coords: %g, %g, %g\n", i, vranks[i], vGlobalID[i],
              coords[3*i], coords[3*i+1], coords[3*i+2]);

      iMOAB_GlobalID * element_global_IDs = (iMOAB_GlobalID*)malloc(nelem[2]*sizeof(iMOAB_GlobalID));
      iMOAB_GlobalID * block_IDs = (iMOAB_GlobalID*)malloc(nelem[2]*sizeof(iMOAB_GlobalID));
      int * ranks = (int*)malloc(nelem[2]*sizeof(int));
      rc = GetVisibleElementsInfo(pid, &nelem[2], element_global_IDs, ranks, block_IDs);
      CHECKRC(rc, "failed to get all elem info");
      for (int i=0; i<nelem[2]; i++)
        printf(" element local id: %3d,  global ID: %3d  rank:%d  block ID: %2d \n", i, element_global_IDs[i],
            ranks[i], block_IDs[i]);
      free(element_global_IDs);
      free(ranks);
      free(block_IDs);

      // get first element connectivity
      int conn[27];
      int nv = 27;
      int eindex = 0;
      rc = GetElementConnectivity(pid, &eindex, &nv, conn);
      CHECKRC(rc, "failed to get first element connectivity");
      printf(" conn for first element: \n");
      for (int i=0; i<nv; i++)
        printf(" %3d", conn[i]);
      printf("\n");


      for (int i=0; i<nblocks[2]; i++)
      {
        printf(" block index: %3d, block ID: %3d \n", i, gbIDs[i] );
        int vertices_per_element, num_elements_in_block;
        rc = GetBlockInfo(pid,  &gbIDs[i] , &vertices_per_element, &num_elements_in_block);
        CHECKRC(rc, "failed to elem block info");
        printf("    has %4d elements with %d vertices per element\n",  num_elements_in_block, vertices_per_element);
        int size_conn= num_elements_in_block*vertices_per_element;
        iMOAB_LocalID * element_connectivity = (iMOAB_LocalID*) malloc (sizeof(iMOAB_LocalID)*size_conn);
        rc = GetBlockElementConnectivities(pid, &gbIDs[i], &size_conn, element_connectivity);
        CHECKRC(rc, "failed to get block elem connectivity");
        int * element_ownership = (int*) malloc (sizeof(int)*num_elements_in_block);

        rc = GetElementOwnership(pid, &gbIDs[i], &num_elements_in_block,  element_ownership);
        CHECKRC(rc, "failed to get block elem ownership");
        iMOAB_GlobalID* global_element_ID = (iMOAB_GlobalID*)malloc(sizeof(iMOAB_GlobalID)*num_elements_in_block);
        iMOAB_LocalID* local_element_ID =(iMOAB_LocalID*)malloc(sizeof(iMOAB_LocalID)*num_elements_in_block);

        rc = GetElementID(pid, &gbIDs[i], &num_elements_in_block, global_element_ID, local_element_ID);
        CHECKRC(rc, "failed to get block elem IDs");
        for (int j=0; j< num_elements_in_block; j++)
        {
          printf("  elem %3d owned by %d gid: %4d lid: %4d  -- ", j, element_ownership[j], global_element_ID[j], local_element_ID[j]);
          for (int k=0; k<vertices_per_element; k++)
            printf( " %5d", element_connectivity[j*vertices_per_element+k]);
          printf("\n");
        }
        free(global_element_ID);
        free(local_element_ID);
        free (element_connectivity);
        free (element_ownership);
      }
      /*
       * query int tag values on vertices
       */
      int * int_tag_vals = (int *) malloc (sizeof(int) * nverts[2]); // for all visible vertices on the rank
      rc = GetIntTagStorage(pid, "INTFIELD", &nverts[2], &entTypes[0],
          int_tag_vals, strlen("INTFIELD"));
      CHECKRC(rc, "failed to get INTFIELD tag");
      printf("INTFIELD tag values:\n");
      for (int i=0; i<nverts[2]; i++)
      {
        printf(" %4d", int_tag_vals[i]);
        if (i%20==19)
          printf("\n");
      }
      printf("\n");
      free(int_tag_vals);

      /*
       * query double tag values on elements
       */
      double * double_tag_vals = (double *) malloc (sizeof(double) * nelem[2]); // for all visible elements on the rank
      rc = GetDoubleTagStorage(pid, "DFIELD", &nelem[2], &entTypes[1],
          double_tag_vals, strlen("DFIELD"));
      CHECKRC(rc, "failed to get DFIELD tag");
      printf("DFIELD tag values: (not exchanged) \n");
      for (int i=0; i<nelem[2]; i++)
      {
        printf(" %f", double_tag_vals[i]);
        if (i%8==7)
          printf("\n");
      }
      printf("\n");
      free(double_tag_vals);

      /* query surface BCs */
      iMOAB_LocalID * surfBC_ID = (iMOAB_LocalID*) malloc (sizeof(iMOAB_LocalID)*nsbc[2]);
      int * ref_surf = (int *)  malloc (sizeof(int)*nsbc[2]);
      int * bc_value = (int *)  malloc (sizeof(int)*nsbc[2]);
      rc = GetPointerToSurfaceBC(pid, &nsbc[2], surfBC_ID, ref_surf, bc_value);
      CHECKRC(rc, "failed to get surf boundary conditions");
      printf(" Surface boundary conditions:\n");
      for (int i=0; i<nsbc[2]; i++)
      {
        printf("  elem_localID %4d  side:%d  BC:%2d\n",surfBC_ID[i] ,ref_surf[i], bc_value[i] );
      }
      free(surfBC_ID);
      free(ref_surf);
      free(bc_value);
      /* query vertex BCs */
      iMOAB_LocalID * vertBC_ID = (iMOAB_LocalID*) malloc (sizeof(iMOAB_LocalID)*ndbc[2]);
      int * vertBC_value = (int *)  malloc (sizeof(int)*ndbc[2]);
      rc = GetPointerToVertexBC(pid, &ndbc[2], vertBC_ID, vertBC_value);
      CHECKRC(rc, "failed to get vertex boundary conditions");
      printf("  Vertex boundary conditions:\n");
      for (int i=0; i<ndbc[2]; i++)
      {
        printf("   vertex %4d   BC:%2d\n",vertBC_ID[i], vertBC_value[i] );
      }
      free(vertBC_ID);
      free(vertBC_value);
    }
    MPI_Barrier(comm); // to avoid printing problems, as we write all procs data
  }

  // free allocated data
  free(coords);
  free (vGlobalID);
  free (vranks);
  char outputFile[] = "fnew.h5m";
  char writeOptions[] ="PARALLEL=WRITE_PART";
  rc = WriteMesh(pid, outputFile, writeOptions,
      strlen(outputFile), strlen(writeOptions) );

  rc = DeregisterApplication(pid);
  CHECKRC(rc, "failed to de-register application");

  rc = iMOABFinalize();
  CHECKRC(rc, "failed to finalize MOAB");

  MPI_Finalize();

  return 0;
}
