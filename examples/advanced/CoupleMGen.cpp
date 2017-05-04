/*
 * Driver to test coupling online, without using IO hdf5 files
 * Will instantiate 2 different meshes, that cover the same domain, and will
 * call a
 */
// MOAB includes
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include "Coupler.hpp"
#include "moab_mpi.h"
#include "ElemUtil.hpp"
#include "moab/MGen.hpp"
#include "moab/ProgOptions.hpp"


using namespace moab;
using std::string;

double physField(double x, double y, double z){

  double out = sin(x + y + z);

  return out;
}

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  Core mcore;
  Interface* mb = &mcore;
  EntityHandle fileset1, fileset2; // for 2 different meshes
  brOpts opts;
  // default options
  opts.A=opts.B=opts.C=1;
  opts.M=opts.N=opts.K=1;
  opts.blockSize = 4;
  opts.xsize = opts.ysize = opts.zsize = 1.;
  opts.ui=CartVect(1.,0,0.);
  opts.uj=CartVect(0.,1.,0.);
  opts.uk=CartVect(0.,0.,1.);
  opts.newMergeMethod =  opts.quadratic =  opts.keep_skins =  opts.tetra = false;
  opts.adjEnts =  opts.parmerge = false;
  opts.nosave = true; // do not save the files

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

  popts.addOpt<void>("no_save,n", "do not save the file", &opts.nosave);

  Coupler::Method method = Coupler::LINEAR_FE;

  double toler = 1.e-6;
  popts.addOpt<double>(string("eps,e"),
      string("tolerance for coupling, used in locating points"), &toler);

  bool writeMeshes=false;
  popts.addOpt<void>("print,p", "write meshes", &writeMeshes);

  popts.parseCommandLine(argc, argv);

  ErrorCode rval = mb->create_meshset(MESHSET_SET, fileset1);MB_CHK_ERR(rval);
  rval = mb->create_meshset(MESHSET_SET, fileset2);MB_CHK_ERR(rval);

  ParallelComm *pc1 = new ParallelComm(mb, MPI_COMM_WORLD);
  MGen * mgen1 = new MGen(mb, pc1, fileset1);

  rval = mgen1->BrickInstance(opts); MB_CHK_ERR(rval); // this will generate first mesh on fileset1

  // set an interpolation tag on source mesh, from phys field
  std::string interpTag("interp_tag");
  Tag tag;
  rval = mb->tag_get_handle(interpTag.c_str(), 1, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT | MB_TAG_DENSE);MB_CHK_ERR(rval);

  Range src_elems;
  rval = pc1->get_part_entities(src_elems, 3);MB_CHK_ERR(rval);
  Range src_verts;
  rval = mb->get_connectivity(src_elems, src_verts);MB_CHK_ERR(rval);
  for(Range::iterator vit=src_verts.begin(); vit!=src_verts.end(); ++vit){
    EntityHandle vert = *vit; //?

    double vertPos[3];
    mb->get_coords(&vert, 1, vertPos);

    double fieldValue =  physField(vertPos[0], vertPos[1], vertPos[2]);

    rval = mb->tag_set_data(tag, &vert, 1, &fieldValue); MB_CHK_ERR(rval);
  }

  // change some options, so it is a different mesh
  int tmp1=opts.K; opts.K=opts.M; opts.M=tmp1; // swap (opts.K, opts.M)
  opts.tetra = !opts.tetra;
  opts.blockSize++;

  ParallelComm *pc2 = new ParallelComm(mb, MPI_COMM_WORLD);
  MGen * mgen2 = new MGen(mb, pc2, fileset2);

  rval = mgen2->BrickInstance(opts); MB_CHK_ERR(rval); // this will generate second mesh on fileset2

  // test the sets are fine
  if (writeMeshes)
  {
    rval = mb->write_file("mesh1.h5m", 0, ";;PARALLEL=WRITE_PART;CPUTIME;PARALLEL_COMM=0;", &fileset1, 1);MB_CHK_SET_ERR(rval, "Can't write in parallel mesh 1");
    rval = mb->write_file("mesh2.h5m", 0, ";;PARALLEL=WRITE_PART;CPUTIME;PARALLEL_COMM=1;", &fileset2, 1);MB_CHK_SET_ERR(rval, "Can't write in parallel mesh 1");
  }


  // Instantiate a coupler, which also initializes the tree
  Coupler mbc(mb, pc1, src_elems, 0);


  // Get points from the target mesh to interpolate
  // We have to treat differently the case when the target is a spectral mesh
  // In that case, the points of interest are the GL points, not the vertex nodes
  std::vector<double> vpos; // This will have the positions we are interested in
  int numPointsOfInterest = 0;

  Range targ_elems;
  Range targ_verts;

  // First get all vertices adj to partition entities in target mesh
  rval = pc2->get_part_entities(targ_elems, 3);MB_CHK_ERR(rval);

  rval = mb->get_adjacencies(targ_elems, 0, false, targ_verts,
                                     Interface::UNION);MB_CHK_ERR(rval);
  Range tmp_verts;
  // Then get non-owned verts and subtract
  rval = pc2->get_pstatus_entities(0, PSTATUS_NOT_OWNED, tmp_verts);MB_CHK_ERR(rval);
  targ_verts = subtract(targ_verts, tmp_verts);
  // get position of these entities; these are the target points
  numPointsOfInterest = (int)targ_verts.size();
  vpos.resize(3*targ_verts.size());
  rval = mb->get_coords(targ_verts, &vpos[0]);MB_CHK_ERR(rval);
  // Locate those points in the source mesh
  std::cout<<"rank "<< proc_id<< " points of interest: " << numPointsOfInterest << "\n";
  rval = mbc.locate_points(&vpos[0], numPointsOfInterest, 0, toler);MB_CHK_ERR(rval);

  // Now interpolate tag onto target points
  std::vector<double> field(numPointsOfInterest);

  rval = mbc.interpolate(method, interpTag, &field[0]);MB_CHK_ERR(rval);

  // compare with the actual phys field
  double err_max = 0;
  for (int i=0; i<numPointsOfInterest; i++)
  {
    double trval = physField(vpos[3*i], vpos[3*i+1],vpos[3*i+2]);
    double err2=fabs(trval-field[i]);
    if (err2>err_max)
      err_max = err2;
  }

  std::cout<<"err max on proc " << proc_id << " is " << err_max << "\n";


  MPI_Finalize();

  return 0;
}
