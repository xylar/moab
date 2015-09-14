
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

  int nverts, nelem, nblocks, nsbc, ndbc;
  rc = GetMeshInfo(  pid, &nverts, &nelem, &nblocks, &nsbc, &ndbc);
  CHECKRC(rc, "failed to get mesh info");
  if (1==rank)
  {
    printf("on rank %d, there are \n"
        "  %d visible vertices\n"
        "  %d visible elements\n"
        "  %d visible blocks\n"
        "  %d visible neumann BCs\n"
        "  %d visible dirichlet BCs\n", rank, nverts, nelem, nblocks, nsbc, ndbc);
  }

  iMOAB_GlobalID * vGlobalID = (iMOAB_GlobalID*)malloc(nverts*sizeof(iMOAB_GlobalID)) ;
  iMOAB_LocalID * vLocalID = (iMOAB_LocalID*)malloc(nverts*sizeof(iMOAB_GlobalID)) ;
  rc = GetVertexID(pid, nverts, vGlobalID, vLocalID );
  CHECKRC(rc, "failed to get vertex id info");

  int * vranks = (int*)malloc(nverts*sizeof(int));
  rc =GetVertexOwnership(pid, nverts, vranks );
  CHECKRC(rc, "failed to get vertex ranks");

  double * coords = (double*) malloc(3*nverts*sizeof(double));
  rc = GetVisibleVerticesCoordinates( pid, 3*nverts, coords);
  CHECKRC(rc, "failed to get coordinates");

  if (1==rank)
  {
    // printf some of the vertex id infos
    int numToPrint = nverts;
    printf("on rank %d some vertex info:\n", rank);
    for (int i=0; i<numToPrint; i++)
      printf(" vertex local id: %3d, rank ID:%d  global ID: %3d  coords: %g, %g, %g\n",vLocalID[i], vranks[i], vGlobalID[i],
            coords[3*i], coords[3*i+1], coords[3*i+2]);
  }

  iMOAB_GlobalID * gbIDs = (iMOAB_GlobalID*) malloc(nblocks*sizeof(iMOAB_GlobalID));
  iMOAB_LocalID *  lbIDs = (iMOAB_LocalID*)  malloc(nblocks*sizeof(iMOAB_LocalID));
  rc = GetBlockID(pid, nblocks, gbIDs, lbIDs);
  CHECKRC(rc, "failed to get block info");
  if (1==rank)
  {
    // printf some of the block info
    printf("on rank %d some block info:\n", rank);
    for (int i=0; i<nblocks; i++)
    {
      printf(" block index: %3d, block ID: %3d \n", lbIDs[i], gbIDs[i] );
      int vertices_per_element, num_elements_in_block;
      rc = GetBlockInfo(pid,  gbIDs[i] , &vertices_per_element, &num_elements_in_block);
      CHECKRC(rc, "failed to elem block info");
      printf("    has %4d elements with %d vertices per element\n",  num_elements_in_block, vertices_per_element);
      int size_conn= num_elements_in_block*vertices_per_element;
      iMOAB_GlobalID * element_connectivity = (iMOAB_GlobalID*) malloc (sizeof(iMOAB_GlobalID)*size_conn);
      rc = GetElementConnectivity(pid, gbIDs[i], size_conn, element_connectivity);
      CHECKRC(rc, "failed to get block elem connectivity");
      int * element_ownership = (int*) malloc (sizeof(int)*num_elements_in_block);

      GetElementOwnership(pid, gbIDs[i], num_elements_in_block,  element_ownership);
      CHECKRC(rc, "failed to get block elem ownership");
      iMOAB_GlobalID* global_element_ID = (iMOAB_GlobalID*)malloc(sizeof(iMOAB_GlobalID)*num_elements_in_block);
      iMOAB_LocalID* local_element_ID =(iMOAB_LocalID*)malloc(sizeof(iMOAB_LocalID)*num_elements_in_block);

      rc = GetElementID(pid, gbIDs[i], num_elements_in_block, global_element_ID, local_element_ID);
      CHECKRC(rc, "failed to get block elem IDs");
      for (int j=0; j< num_elements_in_block; j++)
      {
        printf("  elem %3d owned by %d gid: %4d -- ", j, element_ownership[j], global_element_ID[j]);
        for (int k=0; k<vertices_per_element; k++)
          printf( " %5d", element_connectivity[j*vertices_per_element+k]);
        printf("\n");
      }
      free(global_element_ID);
      free(local_element_ID);
      free (element_connectivity);
      free (element_ownership);
    }

  }


  // free allocated data
  free(coords);
  free (vGlobalID);
  free (vLocalID);
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
