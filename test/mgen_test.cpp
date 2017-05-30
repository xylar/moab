#include "moab/Core.hpp"
#include "moab/MeshGeneration.hpp"
#include "TestUtil.hpp"
#include "moab/ProgOptions.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;
using std::string;

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
  MeshGeneration::BrickOpts opts;
  // default options
  opts.A=opts.B=opts.C=2;
  opts.M=opts.N=opts.K=1;
  opts.blockSize = 4;
  opts.xsize = opts.ysize = opts.zsize = 1.;
  opts.ui=CartVect(1.,0,0.);
  opts.uj=CartVect(0.,1.,0.);
  opts.uk=CartVect(0.,0.,1.);
  opts.newMergeMethod =  opts.quadratic =  opts.keep_skins =  opts.tetra = false;
  opts.adjEnts = opts.parmerge = false;
  opts.GL = 0;

  ProgOptions popts;

  popts.addOpt<int>(string("blockSize,b"),
      string("Block size of mesh (default=4)"), &opts.blockSize);
  popts.addOpt<int>(string("xproc,M"),
      string("Number of processors in x dir (default=1)"), &opts.M);
  popts.addOpt<int>(string("yproc,N"),
      string("Number of processors in y dir (default=1)"), &opts.N);
  popts.addOpt<int>(string("zproc,K"),
      string("Number of processors in z dir (default=1)"), &opts.K);

  popts.addOpt<int>(string("xblocks,A"),
      string("Number of blocks on a task in x dir (default=2)"), &opts.A);
  popts.addOpt<int>(string("yblocks,B"),
      string("Number of blocks on a task in y dir (default=2)"), &opts.B);
  popts.addOpt<int>(string("zblocks,C"),
      string("Number of blocks on a task in x dir (default=2)"), &opts.C);

  popts.addOpt<double>(string("xsize,x"),
      string("Total size in x direction (default=1.)"), &opts.xsize);
  popts.addOpt<double>(string("ysize,y"),
      string("Total size in y direction (default=1.)"), &opts.ysize);
  popts.addOpt<double>(string("zsize,z"),
      string("Total size in z direction (default=1.)"), &opts.zsize);

  popts.addOpt<void>("newMerge,w", "use new merging method", &opts.newMergeMethod);

  popts.addOpt<void>("quadratic,q", "use hex 27 elements", &opts.quadratic);

  popts.addOpt<void>("keep_skins,k", "keep skins with shared entities", &opts.keep_skins);

  popts.addOpt<void>("tetrahedrons,t", "generate tetrahedrons", &opts.tetra);

  popts.addOpt<void>("faces_edges,f", "create all faces and edges", &opts.adjEnts);

  popts.addOpt<int>(string("ghost_layers,g"),
    string("Number of ghost layers (default=0)"), &opts.GL);

  popts.addOpt<void>("parallel_merge,p", "use parallel mesh merge, not vertex ID based merge", &opts.parmerge);

  popts.parseCommandLine(argc, argv);

  ErrorCode rval = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm *pc = new ParallelComm(mb, MPI_COMM_WORLD);
  MeshGeneration * mgen = new MeshGeneration(mb, pc, fileset);
#else
  MeshGeneration * mgen = new MeshGeneration(mb, 0, fileset);
#endif

  rval = mgen->BrickInstance(opts); MB_CHK_ERR(rval);

  return 0;
}
