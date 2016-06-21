#include <iostream>
#include "moab/Interface.hpp"
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Types.hpp"
#include "moab/GeomTopoTool.hpp"


using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char test[] = STRINGIFY(MESHDIR) "/io/test.obj";
static const char shuttle[] = STRINGIFY(MESHDIR) "/io/shuttle.obj";
#endif

GeomTopoTool *myGeomTool;

Tag geom_tag;
Tag name_tag;
Tag id_tag;


void read_file( Interface& moab, const char* input_file);
void test_check_num_entities();
void test_check_meshsets();
void test_check_groups();


int main()
{
  int result = 0;
  
  result += RUN_TEST(test_check_num_entities);
  result += RUN_TEST(test_check_meshsets);
  result += RUN_TEST(test_check_groups);

  return result;
}


void read_file( Interface& moab, const char* input_file )
{
  ErrorCode rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}


void test_check_num_entities()
{
  ErrorCode rval;
  Core core;
  Interface *mbi = &core;
  read_file(core, test);

  // check that number of verts created is 7 
  Range verts;
  int vert_dim = 0;
  rval =  mbi->get_entities_by_dimension(0, vert_dim, verts);
  CHECK_ERR(rval);
  CHECK_EQUAL(7, (int)verts.size());

  // check that number of tris created is 3
  Range tris;
  int tri_dim = 2;
  rval =  mbi->get_entities_by_dimension(0, tri_dim, tris);
  CHECK_ERR(rval);
  CHECK_EQUAL(3, (int)tris.size());

}

void test_check_meshsets()
{
  ErrorCode rval;
  Core core;
  Interface *mbi = &core;
  read_file(core, test);
 
  myGeomTool = new GeomTopoTool(mbi);
  
  Range ent_sets, mesh_sets;
  rval =  mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag);
  CHECK_ERR(rval);
  rval =  mbi->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag, NULL, 1, ent_sets);

  Range::iterator it;
  Range parents, children; 
  int sense;
  int dim, num_surfs = 0, num_vols = 0;

  for (it = ent_sets.begin(); it != ent_sets.end(); ++it)
    {
      rval =  mbi->tag_get_data(geom_tag, &(*it), 1, &dim);
      
      if (dim == 2)
        { 
          num_surfs++;
          
          // check that one parent is created for each surface
          parents.clear();
          rval =  mbi->get_parent_meshsets(*it, parents);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, (int)parents.size());

          // check that sense of surface wrt parent is FORWARD = 1
          rval = myGeomTool->get_sense(*it, *parents.begin(), sense);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, sense);

        }
      else if (dim == 3)
        { 
          num_vols++;
     
          // check that one child is created for each volume
          children.clear();
          rval =  mbi->get_child_meshsets(*it, children);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, (int)children.size());
        }
    }
  
  // check that two surfaces and two volumes are created 
  CHECK_EQUAL(2, num_surfs);
  CHECK_EQUAL(2, num_vols);
}

  Range::iterator it;
  Range parents, children; 
  int sense;
  int dim, num_surfs = 0, num_vols = 0;
  int num_verts, num_tris;

  for (it = ent_sets.begin(); it != ent_sets.end(); ++it)
    {
      rval =  mbi->tag_get_data(geom_tag, &(*it), 1, &dim);
      
      if (dim == 2)
        { 
          num_surfs++;
          
          // check that one parent is created for each surface
          parents.clear();
          rval =  mbi->get_parent_meshsets(*it, parents);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, (int)parents.size());

          // check that sense of surface wrt parent is FORWARD = 1
          rval = myGeomTool->get_sense(*it, *parents.begin(), sense);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, sense);

          // check that each surface set has correct number of entities
          rval = mbi->get_number_entities_by_type(*it, MBTRI, num_tris);
          CHECK_ERR(rval);
          if (num_tris == 1)
            {
              rval = mbi->get_number_entities_by_dimension(*it, 0, num_verts);
              CHECK_ERR(rval);
              CHECK_EQUAL(3, num_verts);
            }
          else if (num_tris == 4)
            {
              rval = mbi->get_number_entities_by_dimension(*it, 0, num_verts);
              CHECK_ERR(rval);
              CHECK_EQUAL(5, num_verts);
            }
              
        }
      else if (dim == 3)
        { 
          num_vols++;
     
          // check that one child is created for each volume
          children.clear();
          rval =  mbi->get_child_meshsets(*it, children);
          CHECK_ERR(rval);
          CHECK_EQUAL(1, (int)children.size());
        }
    }
  
  // check that two surfaces and two volumes are created 
  CHECK_EQUAL(2, num_surfs);
  CHECK_EQUAL(2, num_vols);
}


void test_check_groups()
{
  ErrorCode rval;
  Core core;
  Interface *mbi = &core;
  read_file(core, shuttle);

  // check that number of tris created is 616
  //  170 tris + 223 quads split into 2 tris = 616 
  Range tris;
  int tri_dim = 2;
  rval =  mbi->get_entities_by_dimension(0, tri_dim, tris);
  CHECK_ERR(rval);
  CHECK_EQUAL(616, (int)tris.size());
  
  // check that 11 mesh sets are created
  // 1 for global vert set + 1 for each of 10 groups
  Range ent_sets, mesh_sets;
  rval =  mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag);
  CHECK_ERR(rval);
  rval =  mbi->get_entities_by_type_and_tag(0, MBENTITYSET, &id_tag, NULL, 1, ent_sets);
  
  CHECK_EQUAL(11, (int)ent_sets.size());
}
