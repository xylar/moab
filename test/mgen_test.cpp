#include "moab/Core.hpp"
#include "moab/MGen.hpp"
#include "TestUtil.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif
  Core mcore;
  Interface* mb = &mcore;
  EntityHandle fileset;

  ErrorCode rval = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm *pc = new ParallelComm(mb, MPI_COMM_WORLD);
  MGen * mgen = new MGen(mb, pc, fileset);
#else
  MGen * mgen = new MGen(mb);
#endif

  brOpts opts;
  // default options
  opts.A=opts.B=opts.C=2;
  opts.M=opts.N=opts.K=1;
  opts.blockSize = 4;
  opts.xsize = opts.ysize = opts.zsize = 1.;
  opts.ui=CartVect(1.,0,0.);
  opts.uj=CartVect(0.,1.,0.);
  opts.uk=CartVect(0.,0.,1.);

  opts.newMergeMethod =  opts.quadratic =  opts.keep_skins =  opts.tetra = false;
  opts.adjEnts =  opts.parmerge =  opts.nosave = false;

  rval = mgen->BrickInstance(opts); MB_CHK_ERR(rval);

  return 0;
}
