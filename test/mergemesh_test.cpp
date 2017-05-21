#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MergeMesh.hpp"
#include <iostream>
#include "TestUtil.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#ifndef MOAB_HAVE_HDF5
#error This test requires HDF5 support
#endif

using namespace moab;

std::string meshfile = TestDir + "/16_unmerged_hex.h5m";
std::string meshfile2 = TestDir + "/merge_with_tag.h5m";
std::string meshfile3 = TestDir + "/triangles.h5m";
std::string outfile = "mm_out.h5m";

void mergesimple_test();
void merge_with_tag_test();
void merge_all_test();

#ifdef MOAB_HAVE_MPI
int main(int argc, char** argv)
#else
int main()
#endif
{
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  int result = 0;

  result += RUN_TEST(mergesimple_test);
  result += RUN_TEST(merge_with_tag_test);
  result += RUN_TEST(merge_all_test);

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return result;
}

void mergesimple_test()
{
  ErrorCode rval;
  Core mb;
  Interface* iface = &mb;
  // can be generalized to load user defined input/output file

  rval = iface->load_mesh(meshfile.c_str());
  CHECK_ERR(rval);
  int dim = 3;
  moab::Range ents;
  iface->get_entities_by_dimension(0, dim, ents);

  MergeMesh mm(iface);
  double merge_tol = 1e-3;

  rval = mm.merge_entities(ents, merge_tol);
  CHECK_ERR(rval);

  // Fixed for now

  rval = iface->write_file( outfile.c_str() );
  CHECK_ERR(rval);

  return ;
}

void merge_with_tag_test()
{
  ErrorCode rval;
  Core mb;
  Interface* iface = &mb;
  // can be generalized to load user defined input/output file

  rval = iface->load_mesh(meshfile2.c_str());
  CHECK_ERR(rval);
  int dim = 0;
  moab::Range verts;
  iface->get_entities_by_dimension(0, dim, verts);
  Tag  tag_for_merge;
  rval = iface->tag_get_handle("IDFTAG", tag_for_merge);
  CHECK_ERR(rval);

  MergeMesh mm(iface);
  rval = mm.merge_using_integer_tag(verts, tag_for_merge);
  CHECK_ERR(rval);
  rval = iface->write_file( outfile.c_str() );
  CHECK_ERR(rval);

  verts.clear();
  iface->get_entities_by_dimension(0, dim, verts);
  CHECK_EQUAL( 405, (int)verts.size()) ;

  return;
}

void merge_all_test()
{
  ErrorCode rval;
  Core mb;
  Interface* iface = &mb;

  rval = iface->load_mesh(meshfile3.c_str());
  CHECK_ERR(rval);

  MergeMesh mm(iface);
  double merge_tol = 1e-3;
  rval = mm.merge_all(0, merge_tol); // root set
  CHECK_ERR(rval);
  rval = iface->write_file( outfile.c_str() );
  CHECK_ERR(rval);

  return;
}

