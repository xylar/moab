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
#include "moab/GeomQueryTool.hpp"

using namespace moab;

Interface* MBI;
GeomTopoTool* GTT;
GeomQueryTool* GQT;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

const std::string input_file = TestDir + "/test_geom.h5m";

double eps = 1.0e-6;

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

void gqt_origin_face_rayfire()
{
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {-1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};
  double next_surf_dist;
  EntityHandle next_surf;
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  double expected_next_surf_dist = 5.0;
  CHECK_REAL_EQUAL(expected_next_surf_dist, next_surf_dist, eps);
}

void gqt_outside_face_rayfire()
{
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0}; // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0}; // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 15.0;
  CHECK_REAL_EQUAL(expected_next_surf_dist, next_surf_dist, eps);
}

void gqt_outside_face_rayfire_orient_exit()
{
  GeomQueryTool::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0}; // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0}; // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 15.0;
  CHECK_REAL_EQUAL(expected_next_surf_dist, next_surf_dist, eps);
}

void gqt_outside_face_rayfire_orient_entrance()
{
  GeomQueryTool::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0}; // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0}; // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0.0, -1);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 5.0;
  CHECK_REAL_EQUAL(expected_next_surf_dist, next_surf_dist, eps);
}

void gqt_outside_face_rayfire_history_fail()
{
  GeomQueryTool::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0}; // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0}; // origin at -10 0 0
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  history.reset();

  // ray fired exactly along boundary shared by 2 facets on a single surface,
  // needs two ray_fires to cross, this is expected and ok

  // first ray fire with history
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  // second ray fire with history
  GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);
  // this fire should hit graveyard, i.e. next_surf = 0
  GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  // using history with this geom, there should be no next surface, i.e. 0
  EntityHandle ZERO = 0;
  CHECK_EQUAL(ZERO, next_surf);
}

void gqt_outside_face_rayfire_history()
{
  GeomQueryTool::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = GQT->gttool()->entity_by_id(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0}; // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0}; // origin at -10 0 0
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  history.reset();
  // first ray fire with history
  GQT->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  std::cout << next_surf << " " << history.size() << std::endl;
  // second ray fire with history

  xyz[0] = origin[0] + (next_surf_dist*dir[0]);
  xyz[1] = origin[1] + (next_surf_dist*dir[1]);
  xyz[2] = origin[2] + (next_surf_dist*dir[2]);

  // ray fired execacyl

  GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  GQT->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  // using history with this geom, there should be no next surface, i.e. 0
  EntityHandle ZERO = 0;
  CHECK_EQUAL(ZERO, next_surf);
}

int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST(gqt_setup_test); // setup problem
  // rays fired along cardinal directions
  result += RUN_TEST(gqt_origin_face_rayfire); // point in centre
  result += RUN_TEST(gqt_outside_face_rayfire);
  result += RUN_TEST(gqt_outside_face_rayfire_orient_exit); // fire ray from point outside volume looking for exit intersections
  result += RUN_TEST(gqt_outside_face_rayfire_orient_entrance); // fire ray from point outside volume looking for entrance intersection
  result += RUN_TEST(gqt_outside_face_rayfire_history_fail); // fire ray from point outside geometry using ray history
  result += RUN_TEST(gqt_outside_face_rayfire_history); // fire ray from point outside geometry using ray history

  delete GQT;
  delete GTT;
  delete MBI;

  return result;
}
