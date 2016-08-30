#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/ProgOptions.hpp"
#include "../verdict/moab/VerdictWrapper.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

// Tempest includes
#include "TempestRemapAPI.h"

/* Exit values */
#define SUCCESS 0
#define USAGE_ERROR 1
#define TEMPEST_ERROR 2
#define NOT_IMPLEMENTED 3

using namespace moab;

moab::ErrorCode translate_tempest_mesh(Mesh* mesh, moab::Interface* mb);

// static void usage_error( ProgOptions& opts )
// {
//   opts.printUsage();
// #ifdef MOAB_HAVE_MPI
//   MPI_Finalize();
// #endif
//   exit(USAGE_ERROR);
// }

inline char* create_char_array(const char* s) {
  const int len = strlen(s);
  char *a = new char[len+1];
  a[len]=0;
  memcpy(a,s,len);
  return a;
}

template<typename T>
inline char* create_carray(T val) {
  std::stringstream sstr;
  sstr << val;
  return create_char_array (sstr.str().c_str());
}

enum MeshType { CS=0, RLL=1, ICO=2 };

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;
  int blockSize = 5;
  std::string outFilename="output.exo";
  int meshType=0;
  bool computeDual=false;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  ProgOptions opts;

  opts.addOpt<int>("res,r",
      "Resolution of the mesh (default=5)", &blockSize);
  opts.addOpt<int>("type,t",
      "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2])", &meshType);
  opts.addOpt<std::string>("file,f",
      "Output mesh filename (default=output.exo)", &outFilename);

  // opts.addOpt<int>(string("xproc,M"),
  //     std::string("Number of processors in x dir (default=1)"), &M);
  // opts.addOpt<int>(string("yproc,N"),
  //     std::string("Number of processors in y dir (default=1)"), &N);

  opts.addOpt<void>("dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual);

  opts.parseCommandLine(argc, argv);

  std::vector<char*> progargs;
  Mesh* tempest_mesh;
  switch(meshType) {
    case ICO:
      tempest_mesh = GenerateICOMesh(blockSize, computeDual, outFilename);
      break;
    case RLL:
      tempest_mesh = GenerateRLLMesh(blockSize*2, blockSize, 0.0, 360.0, -90.0, 90.0, false, outFilename);
      break;
    case CS:
    default:
      tempest_mesh = GenerateCSMesh(blockSize, false, outFilename);
      break;
  }

  if (!tempest_mesh) {
    std::cout << "Tempest Mesh is not a complete object; Quitting...";
    exit(TEMPEST_ERROR);
  }

  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb) {
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  ErrorCode rval = translate_tempest_mesh(tempest_mesh, mb);MB_CHK_ERR(rval);

  // mb->print_database();


  delete mb;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(SUCCESS);
}

moab::ErrorCode translate_tempest_mesh(Mesh* mesh, moab::Interface* mb)
{
  return moab::MB_SUCCESS;
}
