/** @example GenLargeMesh.cpp \n
 * \brief Create a large structured mesh, partitioned \n
 *
 *  It shows how to create a mesh on the fly, on multiple processors
 *  Each processor will create its version of a block mesh, partitioned
 *  as AxBxC blocks. Each block will be with blockSize^3 hexahedrons, and
 *  will get a different PARALLEL_PARTITION tag
 *
 *  The number of tasks will be MxNxK, and it must match the mpi size
 *  Each task p will generate its mesh at location (m,n,k), and it is
 *  lexicographic ordering: rank = m + n * M + k * M * N
 *
 *  By default M=1, N=1, K=1, so by default it should be launched on 1 proc
 *  By default, blockSize is 4, and A=2, B=2, C=2, so each task will generate locally
 *  blockSize^3 x A x B x C hexahedrons (value = 64x8 = 512 hexas, in 8 partitions)
 *  (if -t, multiple by 6 for total number of cells/tets)
 *  The total number of partitions will be A*B*C*M*N*K (default 8)
 *
 *  Each part in partition will get a proper tag, and the numbering is in
 *  lexicographic ordering; x direction varies first, then y, then z.
 *  The same principle is used for global id numbering of the nodes and cells.
 *  (x varies first)
 *
 *  The vertices will get a proper global id, which will be used to resolve the
 *  shared entities
 *  The output will be written in parallel, and we will try sizes as big as we can
 *  (up to a billion vertices, because we use int for global ids)
 *
 *  Within each partition, the hexas entity handles will be contiguous, and also the
 *  vertices; The global id will be determined easily, for a vertex, but the entity
 *  handle space will be more interesting to handle, within a partition (we want
 *  contiguous handles within a partition). After merging the vertices, some fragmentation
 *  will occur in the vertices handle space, within each partition.
 *
 *  To run: ./GenLargeMesh
 *
 *  When launched on more procs, you have to make sure
 *  num procs = M*N*K
 *
 *  So you can launch with
 *  mpiexec -np 8 ./GenLargeMesh -M 2 -N 2 -K 2
 *
 *  We also added -q option; it works now only for hexa mesh, it will generate
 *  quadratic hex27 elements
 *
 *  -t option will generate tetrahedrons instead of hexahedra. Each hexahedra is
 *  decomposed into 6 tetrahedrons.
 *
 *  -f option will also generate all edges and faces in the model.
 *  -w will use a newer merging method locally. Merging is necessary to merge
 *  vertices on the local task, and the new method does not use a searching tree,
 *  but rather the global id set on the vertices in a consistent manner
 *
 *  -d and -i options can be used to add some artificial tags on the model;
 *  you can have multiple -d and -i options; -i <tag_name> will set an integer
 *  tag with name tag_name on the vertices; -d < tag_name2> will generate
 *  double tags on cells (3d elements). You can have multiple tags, like
 *  -i tag1 -i tag2 -i tag3 -d tag4
 *
 *  -x, -y, -z options will control the geometric dimensions of the final mesh, in
 *  x, y and z directions.
 *
 *  -o <out_file> controls the name of the output file; it needs to have extension h5m,
 *  because the file is written in parallel.
 *
 *  -k will keep the edges and faces that are generated as part of resolving shared entities
 *  (by default these edges and faces are removed); when -f option is used, the -k option is
 *  enabled too (so no faces and edges are deleted)
 *
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/MeshGeneration.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

#include <time.h>
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{

  bool nosave = false;

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  MeshGeneration::BrickOpts bopts;
    // default options
  bopts.A=bopts.B=bopts.C=2;
  bopts.M=bopts.N=bopts.K=1;
  bopts.blockSize = 4;
  bopts.xsize = bopts.ysize = bopts.zsize = 1.;
  bopts.ui=CartVect(1.,0,0.);
  bopts.uj=CartVect(0.,1.,0.);
  bopts.uk=CartVect(0.,0.,1.);
  bopts.newMergeMethod =  bopts.quadratic =  bopts.keep_skins =  bopts.tetra = false;
  bopts.adjEnts =  bopts.parmerge = false;
  bopts.GL = 0;

  ProgOptions opts;

  opts.addOpt<int>(string("blockSize,b"),
      string("Block size of mesh (default=4)"), &bopts.blockSize);
  opts.addOpt<int>(string("xproc,M"),
      string("Number of processors in x dir (default=1)"), &bopts.M);
  opts.addOpt<int>(string("yproc,N"),
      string("Number of processors in y dir (default=1)"), &bopts.N);
  opts.addOpt<int>(string("zproc,K"),
      string("Number of processors in z dir (default=1)"), &bopts.K);

  opts.addOpt<int>(string("xblocks,A"),
      string("Number of blocks on a task in x dir (default=2)"), &bopts.A);
  opts.addOpt<int>(string("yblocks,B"),
      string("Number of blocks on a task in y dir (default=2)"), &bopts.B);
  opts.addOpt<int>(string("zblocks,C"),
      string("Number of blocks on a task in x dir (default=2)"), &bopts.C);

  opts.addOpt<double>(string("xsize,x"),
      string("Total size in x direction (default=1.)"), &bopts.xsize);
  opts.addOpt<double>(string("ysize,y"),
      string("Total size in y direction (default=1.)"), &bopts.ysize);
  opts.addOpt<double>(string("zsize,z"),
      string("Total size in z direction (default=1.)"), &bopts.zsize);

  opts.addOpt<void>("newMerge,w", "use new merging method", &bopts.newMergeMethod);

  opts.addOpt<void>("quadratic,q", "use hex 27 elements", &bopts.quadratic);

  opts.addOpt<void>("keep_skins,k", "keep skins with shared entities", &bopts.keep_skins);

  opts.addOpt<void>("tetrahedrons,t", "generate tetrahedrons", &bopts.tetra);

  opts.addOpt<void>("faces_edges,f", "create all faces and edges", &bopts.adjEnts);

  opts.addOpt<int>(string("ghost_layers,g"),
  string("Number of ghost layers (default=0)"), &bopts.GL);

  vector<string> intTagNames;
  string firstIntTag;
  opts.addOpt<string>("int_tag_vert,i", "add integer tag on vertices", &firstIntTag);

  vector<string> doubleTagNames;
  string firstDoubleTag;
  opts.addOpt<string>("double_tag_cell,d", "add double tag on cells", &firstDoubleTag);

  string outFileName = "GenLargeMesh.h5m";
  opts.addOpt<string>("outFile,o", "Specify the output file name string (default GenLargeMesh.h5m)", &outFileName);

#ifdef MOAB_HAVE_HDF5_PARALLEL
  bool readb = false;
  opts.addOpt<void>("readback,r", "read back the generated mesh", &readb);

  bool readAndGhost = false;
  opts.addOpt<void>("readAndGhost,G", "read back the generated mesh and ghost one layer", &readAndGhost);
#endif

  opts.addOpt<void>("parallel_merge,p", "use parallel mesh merge, not vertex ID based merge", &bopts.parmerge);

  opts.addOpt<void>("no_save,n", "do not save the file", &nosave);

  opts.parseCommandLine(argc, argv);

  opts.getOptAllArgs("int_tag_vert,i", intTagNames);
  opts.getOptAllArgs("double_tag_cell,d", doubleTagNames);

  if (bopts.adjEnts)
    bopts.keep_skins = true; // Do not delete anything

  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb) {
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  int rank=0, size=1;

#ifdef MOAB_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  EntityHandle fileset;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(rval);
#ifdef MOAB_HAVE_MPI
  ParallelComm *pc = new ParallelComm(mb, MPI_COMM_WORLD);
  MeshGeneration * mgen = new MeshGeneration(mb, pc, fileset);
#else
  MeshGeneration * mgen = new MeshGeneration(mb, 0, fileset);
#endif

  clock_t tt = clock();

  rval = mgen->BrickInstance(bopts); MB_CHK_ERR(rval);

  Range all3dcells;
  rval = mb->get_entities_by_dimension(fileset, 3, all3dcells); MB_CHK_ERR(rval);

  if (0 == rank) {
    cout << "generate local mesh: "
         << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
    cout << "number of elements on rank 0: " << all3dcells.size() << endl;
    cout << "Total number of elements " << all3dcells.size()*size << endl;
    cout << "Element type: " << ( bopts.tetra ? "MBTET" : "MBHEX") << " order:" <<
          ( bopts.quadratic? "quadratic" : "linear" ) << endl;
  }
  Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_SET_ERR(rval, "Can't get all vertices");





  if (!nosave){
#ifdef MOAB_HAVE_HDF5_PARALLEL
    rval = mb->write_file(outFileName.c_str(), 0, ";;PARALLEL=WRITE_PART;CPUTIME;", &fileset, 1);MB_CHK_SET_ERR(rval, "Can't write in parallel");
#else
    // should be a vtk file, actually, maybe make sure of that
    rval = mb->write_file(outFileName.c_str(), 0, "", &fileset, 1);MB_CHK_SET_ERR(rval, "Can't write in serial");
#endif
    if (0 == rank) {
      cout << "write file " << outFileName << " in "
           << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
      tt = clock();
    }
  }
  // delete the mesh that we already have in-memory
  size_t nLocalVerts = verts.size();
  size_t nLocalCells = all3dcells.size();

  mb->delete_mesh();

#ifdef MOAB_HAVE_HDF5_PARALLEL
  if (!nosave && readb)
  {
    // now recreate a core instance and load the file we just wrote out to verify
    Core mb2;
    std::string read_opts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;CPUTIME;");
    if (readAndGhost)
      read_opts+="PARALLEL_GHOSTS=3.0.1;";
    rval = mb2.load_file(outFileName.c_str(), 0, read_opts.c_str());MB_CHK_SET_ERR(rval, "Can't read in parallel");
    if (0 == rank) {
      cout << "read back file " << outFileName << " with options: \n" << read_opts <<
          " in "  << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
      tt = clock();
    }
    moab::Range nverts, ncells;
    rval = mb2.get_entities_by_dimension(0, 0, nverts);MB_CHK_SET_ERR(rval, "Can't get all vertices");
    rval = mb2.get_entities_by_dimension(0, 3, ncells);MB_CHK_SET_ERR(rval, "Can't get all 3d cells elements");

    if (readAndGhost && size > 1)
    {
      // filter out the ghost nodes and elements, for comparison with original mesh
      // first get the parallel comm
      ParallelComm* pcomm2 = ParallelComm::get_pcomm(&mb2, 0);
      if (NULL == pcomm2) MB_SET_ERR(MB_FAILURE, "can't get parallel comm.");
      rval = pcomm2->filter_pstatus(nverts, PSTATUS_GHOST, PSTATUS_NOT);MB_CHK_SET_ERR(rval, "Can't filter ghost vertices");
      rval = pcomm2->filter_pstatus(ncells, PSTATUS_GHOST, PSTATUS_NOT);MB_CHK_SET_ERR(rval, "Can't filter ghost cells");
    }
    if (nverts.size() != nLocalVerts && ncells.size() != nLocalCells ) {
      MB_SET_ERR(MB_FAILURE, "Reading back the output file led to inconsistent number of entities.");
    }

    // delete the mesh that we already have in-memory
    mb2.delete_mesh();
  }
#endif

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
