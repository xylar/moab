#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/GeomTopoTool.hpp"


using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.obj";
#endif

Tag category_tag;
Tag geom_tag;
Tag name_tag;
Tag obj_name_tag;
Tag dim_tag;
Tag id_tag;

Core core;
Interface *mbi= &core;

Errorcode check_vertices()
{
  mbi->get_number_entities_by_dimension(object_meshset);
}

int main()
{
  ReadOBJ *robj = new ReadOBJ(mbi, rval);
  rval = robj->load_file(filename);
  CHKERR(rval);
}



