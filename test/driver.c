
#include "mpi.h"
#include  "../src/moab/imoab.h"
#include <string.h>

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s", message); return 1;}
int main(int argc, char * argv[])
{

  char * filen = "p8ex1.h5m";
  if (argc>1)
    filen = argv[1];

  ErrCode rc = iMOABInitialize(argc, argv);

  CHECKRC(rc, "failed to initialize MOAB");
  int num_global_vertices=0, num_global_elements=0, num_dimension=0, num_parts=0;
  rc = ReadHeaderInfo ( filen, &num_global_vertices, &num_global_elements, &num_dimension,
      &num_parts, (int)strlen(filen) );

  CHECKRC(rc, "failed to read header info");

  printf("file %s has %d vertices, %d elements, %d parts in partition\n", filen,
      num_global_vertices, num_global_elements, num_parts);

  rc = iMOABFinalize();
  CHECKRC(rc, "failed to finalize MOAB");
  return 0;
}
