#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

#include "moab/GeomQueryTool.hpp"
#include "moab/GeomTopoTool.hpp"

using namespace moab;

Core* MBI;
GeomTopoTool* GTT;
GeomQueryTool* GQT;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

const std::string input_file = TestDir + "/test_geom.h5m";

void gqt_setup_test()
{
  MBI = new Core();
  ErrorCode rval = MBI->load_file(input_file.c_str());
  CHECK_ERR(rval);

  GTT = new GeomTopoTool(MBI);
  GQT = new GeomQueryTool(GTT);
  rval = GQT->initialize();
  CHECK_ERR(rval);
}

void gqt_point_in()
{
  int result = 0;
  int expected_result = 1;
  double xyz[3] = {0.0, 0.0, 0.0};
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  ErrorCode rval = GQT->point_in_volume(vol_h, xyz, result);
  CHECK_ERR(rval);
  CHECK_EQUAL(expected_result, result);
}

int gqt_point_in_vol_dir(double origin[3], double dir[3], int vol_idx)
{
  int result = 0;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  // normalise the vector
  double dir_norm = (dir[0]*dir[0]) + (dir[1]*dir[1]) + (dir[2]*dir[2]);

  dir[0] = dir[0] / sqrt(dir_norm);
  dir[1] = dir[1] / sqrt(dir_norm);
  dir[2] = dir[2] / sqrt(dir_norm);

  ErrorCode rval = GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  CHECK_ERR(rval);

  xyz[0] = origin[0] + (next_surf_dist*dir[0]);
  xyz[1] = origin[1] + (next_surf_dist*dir[1]);
  xyz[2] = origin[2] + (next_surf_dist*dir[2]);

  std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;

  rval = GQT->point_in_volume(vol_h, xyz, result, dir);
  CHECK_ERR(rval);
  return result;
}

void gqt_point_in_vol_1()
{
  double dir[3] = {-1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};
  int vol_idx = 1;
  int expected_result = 1;

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);
  CHECK_EQUAL(expected_result, result);
}

void gqt_point_in_vol_2()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_in_vol_3()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, -1.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_in_vol_4()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 1.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_in_vol_5()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 0.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_in_vol_6()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 0.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_1()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_2()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, 1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_3()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_4()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, 1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_5()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, -1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_6()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, -1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_7()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, -1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

void gqt_point_on_corner_8()
{
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, -1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = gqt_point_in_vol_dir(origin, dir, vol_idx);

  CHECK_EQUAL(expected_result, result);
}

int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST(gqt_setup_test); // setup problem
  result += RUN_TEST(gqt_point_in); // point in centre
  // rays fired along cardinal directions
  result += RUN_TEST(gqt_point_in_vol_1); // point in centre
  result += RUN_TEST(gqt_point_in_vol_2); // point in centre
  result += RUN_TEST(gqt_point_in_vol_3); // point in centre
  result += RUN_TEST(gqt_point_in_vol_4); // point in centre
  result += RUN_TEST(gqt_point_in_vol_5); // point in centre
  result += RUN_TEST(gqt_point_in_vol_6); // point in centre
  // rays fired at nodes
  result += RUN_TEST(gqt_point_on_corner_1);
  result += RUN_TEST(gqt_point_on_corner_2);
  result += RUN_TEST(gqt_point_on_corner_3);
  result += RUN_TEST(gqt_point_on_corner_4);

  //result += RUN_TEST(dagmc_point_in({0.0, 0.0, 5.0}); // point in centre
	//result += RUN_TEST(dagmc_point_in({0.0, 0.0, -5.0}); // point in centre
	//result += RUN_TEST(dagmc_point_in({0.0, 5.0, 0.0}); // point in centre
	//result += RUN_TEST(dagmc_point_in({0.0, -5.0, 0.0}); // point in centre
	//result += RUN_TEST(dagmc_point_in({5.0, 0.0, 0.0}); // point in centre
	//result += RUN_TEST(dagmc_point_in({-5.0, 0.0, 0.0}); // point in centre

  delete GQT;
  delete GTT;
  delete MBI;

  return result;
}
