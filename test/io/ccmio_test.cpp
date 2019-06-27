#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"

using namespace moab;

std::string cubfile = TestDir + "/io/singlecyl.cub";
std::string ccmgfiler = TestDir + "/io/singlecyl.ccmg";
std::string ccmgfilew = "singlecyl_tmp.ccmg";

void test_read();
void test_write();

int main()
{
  int result = 0;

  result += RUN_TEST(test_write);
  result += RUN_TEST(test_read);

  return result;
}

void test_write()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(cubfile.c_str());
  if (MB_SUCCESS != rval) std::cerr << "Trouble reading file " << cubfile << std::endl;
  CHECK_ERR(rval);

  rval = mb.write_file(ccmgfilew.c_str());
  if (MB_SUCCESS != rval) std::cerr << "Trouble writing file " << ccmgfilew << std::endl;
  CHECK_ERR(rval);
}

void test_read()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(ccmgfiler.c_str());
  CHECK_ERR(rval);
}

