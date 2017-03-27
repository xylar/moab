#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"

using namespace moab;

GeomTopoTool* GTT;
GeomQueryTool* GQT;
Interface* MBI;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#else
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#endif

void gqt_load_file() 
{
  MBI = new moab::Core();
  ErrorCode rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true);
  GQT = new GeomQueryTool(GTT);
}

void gqt_load_file_dagmc_build_obb() 
{
  /* 1 - Test with external moab, load file in GTT/GQT*/
  // make new moab core
  ErrorCode rval;

  // load a file
  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI);
  GQT = new GeomQueryTool(GTT);
  rval = GQT->initialize();
  CHECK_ERR(rval);

  // delete instances
  delete GTT;
  delete GQT;
  delete MBI;

}

void gqt_load_file_dagmc_via_moab_build_obb() {
  /* 2 - Test with external moab, load file in MOAB*/
  // build implicit complement and obb trees
  ErrorCode rval;

  // load file
  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the file
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI);
  GQT = new GeomQueryTool(GTT);
  // create obb trees
  EntityHandle implicit_complement;
  rval = GQT->gttool()->get_implicit_complement(implicit_complement, true);
  CHECK_ERR(rval);
  rval = GQT->gttool()->construct_obb_trees();
  CHECK_ERR(rval);			      

  delete GTT;
  delete GQT;
  delete MBI;

}

void gqt_load_file_dagmc_internal_build_obb() {
  /* 3 - Test with internal moab, load file in MOAB*/
  // find geomsets, build implicit complement, and construct obb trees
  ErrorCode rval;

  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  EntityHandle implicit_complement;
  rval = GQT->gttool()->get_implicit_complement(implicit_complement, true);
  CHECK_ERR(rval);
  rval = GQT->gttool()->construct_obb_trees();
  CHECK_ERR(rval);			      

  delete GTT;
  delete GQT;
  delete MBI;

}

void gqt_test_obb_retreval() {
    // make new dagmc
  std::cout << "test_obb_retreval" << std::endl;

  ErrorCode rval;
  // load a file
  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  GQT->initialize();
  CHECK_ERR(rval);			      

  //write to file
  rval = MBI->write_file("fcad");
  CHECK_ERR(rval);

  // delete instances
  delete GTT;
  delete GQT;
  delete MBI;

  // re-load written file
  MBI = new moab::Core();
  rval = MBI->load_file("fcad"); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  rval = GQT->initialize();
  CHECK_ERR(rval);			      

  // remove file
  remove("fcad");
  // delete instances
  delete GTT;
  delete GQT;
  delete MBI;
  
}


void gqt_create_impl_compl() {
  EntityHandle implicit_complement;
  ErrorCode rval = GQT->gttool()->get_implicit_complement(implicit_complement, true);
  CHECK_ERR(rval);
}

void gqt_build_obb() 
{  
  ErrorCode rval = GQT->gttool()->construct_obb_trees();
  CHECK_ERR(rval);
}

void gqt_num_vols()
{
  int expect_num_vols = 2;
  int num_vols = GQT->gttool()->num_ents_of_dim(3); 
  CHECK_EQUAL(expect_num_vols, num_vols);

}

void gqt_point_in()
{
  int result = 0;
  int expect_result = 1;
  int vol_idx = 1;
  double xyz[3] = {0.0, 0.0, 0.0};
  EntityHandle vol_h;

  // load file and create instances
  MBI = new moab::Core();
  ErrorCode rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  // initialize for geom queries
  GQT->initialize();

  // get volume
  vol_h = GQT->gttool()->entity_by_id(3,1);
  // test point in volume
  GQT->point_in_volume(vol_h, xyz, result);
  CHECK_ERR(rval);			      
  CHECK_EQUAL(expect_result, result);
}

void gqt_test_obb_retreval_rayfire() {
  // make new dagmc
  std::cout << "test_obb_retreval and ray_fire" << std::endl;
  
  ErrorCode rval;

  // load file and initialize
  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  GQT->initialize();
  CHECK_ERR(rval);			      

  // write to file
  rval = MBI->write_file("fcad");
  CHECK_ERR(rval);

  // delete instances
  delete GTT;
  delete GQT;
  delete MBI;

  // load the file that was written and initialize
  MBI = new moab::Core();
  rval = MBI->load_file("fcad"); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  rval = GQT->initialize();
  CHECK_ERR(rval);			      

  // remove generated test file
  remove("fcad");

  // now perform full ray fire
  double eps = 1.e-6;
  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h;

  // get volume handle to fire on
  vol_h = GQT->gttool()->entity_by_id(3,vol_idx);
  // test ray fire
  rval = GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  CHECK_ERR(rval);
  CHECK_REAL_EQUAL(expect_next_surf_dist, next_surf_dist, eps);

  // delete instances
  delete GTT;
  delete GQT;
  delete MBI;
}

void gqt_rayfire()
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h;
  ErrorCode rval;

  MBI = new moab::Core();
  rval = MBI->load_file(input_file); // open the
  CHECK_ERR(rval);
  GTT = new GeomTopoTool(MBI,true,0);
  GQT = new GeomQueryTool(GTT);
  GQT->initialize();
  CHECK_ERR(rval);			      

  vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  rval = GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  CHECK_ERR(rval);
  CHECK_REAL_EQUAL(expect_next_surf_dist, next_surf_dist, eps);
  
}

void gqt_closest_to()
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {-6.0, 0.0, 0.0};
  double distance; // distance from point to nearest surface
  double expect_distance = 1.0;
  EntityHandle vol_h;
  ErrorCode rval;
  
  // get the volume handle
  vol_h = GQT->gttool()->entity_by_id(3, vol_idx);

  // test closest to location
  rval = GQT->closest_to_location(vol_h, xyz, distance);
  CHECK_ERR(rval);
  
  // distance should be 1.0 cm
  CHECK_REAL_EQUAL(expect_distance, distance, eps);

}

void gqt_test_boundary()
{
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  int surf_idx = 1;
  EntityHandle surf_h = GQT->gttool()->entity_by_id(2, surf_idx);
  
  double xyz[3] = {0.0, 0.0, 5.0};
  double dir[3] = {0.0, 0.0, 1.0};
  int result;
  int expect_result = 0;
  
  // test vol boundary
  ErrorCode rval = GQT->test_volume_boundary(vol_h, surf_h, xyz, dir, result);
  CHECK_ERR(rval);
  // check ray leaving volume
  CHECK_EQUAL(expect_result, result);
}
  
int main(int /* argc */, char** /* argv */)
{
  int result = 0;
  
  result += RUN_TEST(gqt_load_file); // test ray fire
  result += RUN_TEST(gqt_build_obb); // build the obb
  result += RUN_TEST(gqt_create_impl_compl); // build the obb
  result += RUN_TEST(gqt_num_vols); // make sure the num of vols correct
  result += RUN_TEST(gqt_load_file_dagmc_build_obb); //
  result += RUN_TEST(gqt_load_file_dagmc_via_moab_build_obb); //
  result += RUN_TEST(gqt_load_file_dagmc_internal_build_obb); // 
  result += RUN_TEST(gqt_test_obb_retreval); // check that we are retreving loaded obbs
  result += RUN_TEST(gqt_test_obb_retreval_rayfire); // check that we can ray fire on loaded obbs
  result += RUN_TEST(gqt_point_in); // check entity by point
  result += RUN_TEST(gqt_rayfire); // ensure ray fire distance is correct
  result += RUN_TEST(gqt_closest_to); // check the distance to surface nearest point
  result += RUN_TEST(gqt_test_boundary); // check particle entering leaving

  // cleanup instances
  delete GQT;
  delete GTT;
  delete MBI;
  
  return result;
}
