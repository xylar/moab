
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"
#include "ReadCGM.hpp"
#include "moab/FileOptions.hpp"
using namespace moab;

#ifdef HAVE_OCC_STEP
std::string input_cube = TestDir + "/io/cube.stp";
#else
std::string input_cube = TestDir + "/io/cube.sat";
#endif

#ifdef HAVE_OCC_STEP
std::string input_cone = TestDir + "/io/cone.stp";
#else
std::string input_cone = TestDir + "/io/cone.sat";
#endif

// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// List of tests in this file
void test_cube_curve_facet();
void test_cone_curve_facet();


void read_file( Interface* moab, bool curve_fatal, const char* input_file, ErrorCode check_val,
		int &curve_fail, int &surface_fail )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();
  EntityHandle fs = 0;

  // set the options
  std::string options;
  #define OPTION_APPEND(X) { if( options.length() ) options += ";"; options += (X); }

  OPTION_APPEND( "CGM_ATTRIBS=no" );
  if(curve_fatal)
    OPTION_APPEND( "FATAL_ON_CURVES" );
  OPTION_APPEND( "VERBOSE_CGM_WARNINGS" );

  // set the file options
  FileOptions opts(options.c_str());

  // new ReadCGM instance
  ReadCGM *RCGM = new ReadCGM(moab);
  ErrorCode rval = RCGM->load_file( input_file, &fs, opts );
  //  CHKERR(rval);
  curve_fail = RCGM->get_failed_curve_count();
  surface_fail = RCGM->get_failed_surface_count();
  //  std::cout << curve_fail << " " << surface_fail << std::endl;

  CHECK_EQUAL(rval,check_val);
  CHECK_EQUAL(curve_fail,0);
  CHECK_EQUAL(surface_fail,0);
}

// Gets the vertex entities from a simple cube file load and checks that the
// correct number of them exist.
void test_cube_curve_facet()
{
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  int curve_fail = 0;
  int surface_fail = 0;
  // should succeed since we are not fatal error on curves
  read_file( mb,false,input_cube.c_str(), MB_SUCCESS,curve_fail, surface_fail );
}

// Gets the vertex entities from a simple cube file load and checks that the
// correct number of them exist.
void test_cone_curve_facet()
{
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  int curve_fail = 0;
  int surface_fail = 0;
  // should expect MB_FAILURE since we fail on curves
  read_file( mb, true,input_cone.c_str(), MB_FAILURE, curve_fail, surface_fail );
}


//void delete_mesh_test();
int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST(test_cube_curve_facet);
  #ifndef HAVE_OCC_STEP
  result += RUN_TEST(test_cone_curve_facet);
  #endif
  return result;
}
