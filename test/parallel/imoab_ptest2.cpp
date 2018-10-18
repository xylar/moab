/*
 * imoab_ptest_C.cpp
 *
 *  push parallel mesh into moab, using iMoab API, C version
 *      This program shows how to push a mesh into MOAB in parallel using iMoab, with sufficient
 *     information to resolve boundary sharing and exchange a layer of ghost information.
 *
 *      After resolving the sharing, mesh is saved to a file, then read back again, in parallel
 *       the mesh is a quad mesh, 2x2 on each task,
 *      by default, this test is run on 2 processors

 *      Test is similar to itaps/iMeshP_unit_tests.cpp, and with imoab_ptest.F, fortran 77 version
 *
 *
 *      Create a mesh:
 *      Groups of four quads will be arranged into parts as follows:
 *      +------+------+------+------+------+-----
 *      |             |             |
 *      |             |             |
 *      +    Part 0   +    Part 2   +    Part 4
 *      |             |             |
 *      |             |             |
 *      +------+------+------+------+------+-----
 *      |             |             |
 *      |             |             |
 *      +    Part 1   +    Part 3   +    Part 5
 *      |             |             |
 *      |             |             |
 *      +------+------+------+------+------+-----
 *
 *      Vertices will be enumerated as follows:
 *      1------6-----11-----16-----21-----26----- > x  x from 0, 1, 2
 *      |             |             |
 *      |             |             |
 *      2      7     12     17     22     27
 *      |             |             |
 *      |             |             |
 *      3------8-----13-----18-----23-----28-----
 *      |             |             |
 *      |             |             |
 *      4      9     14     19     24     29
 *      |             |             |
 *      |             |             |
 *      5-----10-----15-----20-----25-----30-----
 *      |
 *      y varies from 0 to 4
 *
 *      Processor 0 will have vertices 1, 2, 3, 6, 7, 8, 11, 12, 13
 *        and 4 quads, with ids from 1 to 4, and connectivity
 *         1, 2, 7, 6;  2, 3, 8, 7; 6 ,7, 12, 11; 7, 8, 13, 12
 *      Processor 1 will have vertices 3, 4, 5, 8, 9, 10, 13, 14, 15
 *        and 4 quads, with ids from 5 to 8, and connectivity
 *         3, 4, 8, 9;  4, 5, 10, 9; 8, 9, 14, 13; 9, 10, 15, 14,
 *      and so on
 *
 *      Vertex Global IDs will be used to resolve sharing
 *
 *      Element IDs will be [4*rank+1,4*rank+5]
 *      */

/*
 * #define ERROR(rval) if (0 .ne. rval) call exit(1)
*/

#include "moab/MOABConfig.h"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include <stdio.h>
#include <string.h>
#include <vector>

#define ERROR(rc, msg) if (0 != rc)  { printf ("error %s", msg); return 1;}

int main(int argc, char ** argv)
{
  int my_id, num_procs, ix, iy, numv, nume;
  int dime, lco, mbtype, blockid, npe;
  char appname[10]="IMTEST";
/* coordinates for 9 vertices */
  double  coordinates[27] , deltax, deltay;
  int ids[9]={1, 2, 3, 6, 7, 8, 11, 12, 13};
  int connec[16] = {
      1, 2, 5, 4,
      2, 3, 6, 5,
      4, 5, 8, 7,
      5, 6, 9, 8 };
/*      used for ghosting */
  int dimgh, bridge, num_layers;
  double coords_core[27]= {
    0., 0., 0.,
    0., 1., 0.,
    0., 2., 0.,
    1., 0., 0.,
    1., 1., 0.,
    1., 2., 0.,
    2., 0., 0.,
    2., 1., 0.,
    2., 2., 0.
  };
  int appID, compid;
  iMOAB_AppID pid=&appID;

  MPI_Init( &argc, &argv );
  MPI_Comm comm=MPI_COMM_WORLD;

  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &my_id);

  ErrCode rc = iMOAB_Initialize(argc, argv); ERROR(rc,"can't initialize");
  if (!my_id)
    printf(" I'm process %d out of %d\n", my_id, num_procs);

  appname[7]=0; // make sure it is NULL
  compid = 8; // given number, external application number
  rc = iMOAB_RegisterApplication( appname, &comm, &compid, pid); ERROR(rc, " can't register app");
  /* create first 9 vertices,  in a square 3x3; */
  deltax = (my_id/2) * 2.;
  deltay = (my_id%2) * 2.;
  ix = (my_id/2) * 10;
  iy =my_id%2 * 2;
  for (int i=0; i<9; i++)
  {
    coordinates[ 3*i     ] = coords_core [3*i  ] + deltax;
    coordinates[ 3*i + 1 ] = coords_core [3*i+1] + deltay;
    coordinates[ 3*i + 2 ] = coords_core [3*i+2];

/*       translate the ids too, by multiples of 10 or add 2, depending on my_id*/
    ids[i]= ids[i] + ix+iy;
  }
  numv = 9;
  nume = 4;
  lco = numv*3;
  dime = 3;
  rc = iMOAB_CreateVertices(pid, &lco, &dime, coordinates); ERROR(rc, "can't create vertices");
  /* create now 4 quads with those 9 vertices*/
  mbtype = 3;
  blockid  =100;
  npe = 4;
  rc = iMOAB_CreateElements(pid, &nume, &mbtype, &npe, connec, &blockid); ERROR(rc, "can't create elements");

  rc = iMOAB_ResolveSharedEntities( pid, &numv, ids ); ERROR(rc, "can't resolve shared ents");

  // test iMOAB_ReduceTagsMax on a tag
  int tagType = DENSE_INTEGER;
  int num_components = 1;
  int tagIndex = 0; // output

  rc = iMOAB_DefineTagStorage(pid, "INTFIELD", &tagType, &num_components, &tagIndex,  strlen("INTFIELD") );
  ERROR(rc, "failed to get tag INTFIELD ");
  //set some values
  std::vector<int> valstest(numv);
  for (int k=0; k<numv; k++)
  {
    valstest[k] = my_id+k;
  }
  int num_tag_storage_length = numv*num_components;
  int entType = 0; // vertex
  rc = iMOAB_SetIntTagStorage ( pid, "INTFIELD", &num_tag_storage_length, &entType, &valstest[0], strlen("INTFIELD") );
  ERROR(rc, "failed to set tag INTFIELD ");

  rc = iMOAB_ReduceTagsMax ( pid, &tagIndex, &entType );
  ERROR(rc, "failed reduce tags max ");

  /* see ghost elements */
  dimgh = 2; /* will ghost quads, topological dim 2 */
  bridge = 0 ; /* use vertex as bridge */
  num_layers = 1 ; /* so far, one layer only */
  rc = iMOAB_DetermineGhostEntities( pid, &dimgh, &num_layers, &bridge); ERROR(rc, "can't determine ghosts");

/*     write out the mesh file to disk, in parallel, if h5m*/
#ifdef MOAB_HAVE_HDF5_PARALLEL
  char outfile[32]="whole.h5m";
  char wopts[100]="PARALLEL=WRITE_PART";
  rc = iMOAB_WriteMesh(pid, outfile, wopts, 9, 19); ERROR(rc,"can't write mesh");
#endif
/*     all done. de-register and finalize */
  rc = iMOAB_DeregisterApplication(pid); ERROR(rc, "can't de-register app");

  rc = iMOAB_Finalize(); ERROR(rc, "can't finalize");

  MPI_Finalize();
  return 0;
}

