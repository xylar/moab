#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"

using namespace moab;

std::string cgnsfile = TestDir + "/io/2d_naca0012.cgns";
std::string cgnsfilew = TestDir + "/io/test.cgns";

void test_read_write();

int main()
{
  int result = 0;

  result += RUN_TEST(test_read_write);

  return result;
}

void test_read_write()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(cgnsfile.c_str());
  if (MB_SUCCESS != rval) std::cerr << "Trouble reading file " << cgnsfile << std::endl;
  CHECK_ERR(rval);

  rval = mb.write_file(cgnsfilew.c_str());
  if (MB_SUCCESS != rval) std::cerr << "Trouble writing file " << cgnsfilew << std::endl;
  CHECK_ERR(rval);
}

