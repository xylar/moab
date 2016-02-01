#include <iostream>
#include "moab/Interface.hpp"
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
//#include "ReadOBJ.hpp"
#include "moab/Types.hpp"
#include "moab/GeomTopoTool.hpp"


using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
//static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.obj";
//static const char input_cube = 'cube.obj';
#endif

Tag category_tag;
Tag geom_tag;
Tag name_tag;
Tag obj_name_tag;
Tag dim_tag;
Tag id_tag;

Core core;
Interface *mbi= &core;

/*
Errorcode check_vertices()
{
  mbi->get_number_entities_by_dimension(object_meshset);
}
*/

int main()
{
  ErrorCode rval;
  const char* filename = "cube.obj";
  ReadOBJ *robj = new ReadOBJ(mbi);
  rval = robj->load_file(&filename);
//  rval = robj->load_file(&input_cube);
  CHKERR(rval);
}



