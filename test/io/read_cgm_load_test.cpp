
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef HAVE_OCC_STEP
std::string input_file = TestDir + "/io/dum.stp";
#else
std::string input_file = TestDir + "/io/dum.sat";
#endif

// Checks that a file can be loaded twice without errors
void read_multiple_test()
{
  Core mb;

  ErrorCode rval = mb.load_file( input_file.c_str() );
  CHECK_ERR(rval);
  // second load
  rval = mb.load_file( input_file.c_str() );
  CHECK_ERR(rval);

}

int main(int /* argc */, char** /* argv */)
{
  int result = RUN_TEST( read_multiple_test );

  return result;
}

