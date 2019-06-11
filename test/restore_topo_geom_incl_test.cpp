#include "moab/Interface.hpp"
#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include <iostream>
#include <map>
#include <set>

using namespace moab;

Tag category_tag;
Tag geom_tag;
Tag name_tag;
Tag obj_name_tag;
Tag dim_tag, id_tag;

bool check_tree ( Interface *mbi, GeomTopoTool *GTT, std::map< int, std::set<int> > &ref_map );
ErrorCode get_all_handles(Interface *mbi);
Range get_children_by_dimension(Interface *mbi, EntityHandle parent, int desired_dimension);
void heappermute(Interface *mbi, int v[], int n, std::map< int, std::set<int> > ref_map, int len);
void swap(int *x, int *y);
void get_cube_info( int cube_id, std::vector<double> &scale, std::vector<double> &trans );
void test_two_cubes();
void test_three_cubes();
void test_four_cubes();

ErrorCode build_cube( Interface *mbi,
                      std::vector<double> scale_vec,
                      std::vector<double> trans_vec,
                      int    object_id,
			          EntityHandle &volume	)
{
  GeomTopoTool *GTT = new GeomTopoTool(mbi);

  ErrorCode rval;

  // Define a 1x1x1 cube centered at orgin

  // coordinates of each corner
  const double coords[] = {
    0.5, -0.5, -0.5,
    0.5,  0.5, -0.5,
   -0.5,  0.5, -0.5,
   -0.5, -0.5, -0.5,
    0.5, -0.5,  0.5,
    0.5,  0.5,  0.5,
   -0.5,  0.5,  0.5,
   -0.5, -0.5,  0.5 };

  // connectivity of 2 triangles per
  //  each face of the cube
  const int connectivity[] = {
    0, 3, 1,  3, 2, 1, // -Z
    0, 1, 4,  5, 4, 1, // +X
    1, 2, 6,  6, 5, 1, // +Y
    6, 2, 3,  7, 6, 3, // -X
    0, 4, 3,  7, 3, 4, // -Y
    4, 5, 7,  5, 6, 7 // +Z
  };


  // Create the geometry
  const int num_verts = 8;
  const int num_tris = 12;
  EntityHandle verts[num_verts], tris[num_tris], surf;

  rval = mbi->create_meshset( MESHSET_SET, surf ); MB_CHK_ERR(rval);
/*
  // scale coords
  int i;
  double scaled_coords[24];
  for ( i = 0; i < num_verts; i++ )
    {
      scaled_coords[3*i]   = coords[3*i]*scale_vec[0];
      scaled_coords[3*i+1] = coords[3*i+1]*scale_vec[1];
      scaled_coords[3*i+2] = coords[3*i+2]*scale_vec[2];
    }

  // translate coords
  double trans_coords[24];
  for ( i = 0; i < num_verts; i++ )
    {
      trans_coords[3*i]   = scaled_coords[3*i] + trans_vec[0];
      trans_coords[3*i+1] = scaled_coords[3*i+1] + trans_vec[1];
      trans_coords[3*i+2] = scaled_coords[3*i+2] + trans_vec[2];
    }
*/
  // transform coords-- scale and translate
  double trans_coords[24];
  for ( int i = 0; i < num_verts; i++ )
    {
      trans_coords[3*i]   = coords[3*i]*scale_vec[0]   + trans_vec[0];
      trans_coords[3*i+1] = coords[3*i+1]*scale_vec[1] + trans_vec[1];
      trans_coords[3*i+2] = coords[3*i+2]*scale_vec[2] + trans_vec[2];
    }

  // create vertices and add to meshset
  for ( int i = 0; i < num_verts; ++i)
    {
      rval = mbi->create_vertex( trans_coords + 3*i, verts[i] ); MB_CHK_ERR(rval);

      rval = mbi->add_entities( surf, &verts[i], 1 ); MB_CHK_ERR(rval);

    }

  // create triangles and add to meshset
  for ( int i = 0; i < num_tris; ++i)
    {
      const EntityHandle conn[] = { verts[connectivity[3*i  ]],
                                    verts[connectivity[3*i+1]],
                                    verts[connectivity[3*i+2]] };
      rval = mbi->create_element( MBTRI, conn, 3, tris[i] ); MB_CHK_ERR(rval);

      rval = mbi->add_entities( surf, &tris[i], 1 ); MB_CHK_ERR(rval);
    }


  // set name, id, geom, and category tags for SURFACE
  rval = mbi->tag_set_data( name_tag, &surf, 1, "Surface\0" ); MB_CHK_ERR(rval);
  std::string object_name;
  rval = mbi->tag_set_data( obj_name_tag, &surf, 1, object_name.c_str() ); MB_CHK_ERR(rval);
  rval = mbi->tag_set_data( id_tag, &surf, 1, &object_id ); MB_CHK_ERR(rval);
  int two = 2;
  rval = mbi->tag_set_data( geom_tag, &surf, 1, &(two) ); MB_CHK_ERR(rval);
  rval = mbi->tag_set_data( category_tag, &surf, 1, "Surface\0" ); MB_CHK_ERR(rval);

  // create volume meshset associated with surface meshset
  // EntityHandle volume;
  rval = mbi->create_meshset( MESHSET_SET, volume ); MB_CHK_ERR(rval);

  // set name, id, geom, and category tags for VOLUME
  rval = mbi->tag_set_data( name_tag, &volume, 1, "Volume\0" ); MB_CHK_ERR(rval);
  rval = mbi->tag_set_data( obj_name_tag, &surf, 1, object_name.c_str() ); MB_CHK_ERR(rval);
  rval = mbi->tag_set_data( id_tag, &volume, 1, &(object_id) ); MB_CHK_ERR(rval);
  int three = 3;
  rval = mbi->tag_set_data( geom_tag, &volume, 1, &(three) ); MB_CHK_ERR(rval);
  rval = mbi->tag_set_data( category_tag, &volume, 1, "Volume\0" ); MB_CHK_ERR(rval);


  // set surface as child of volume
  rval = mbi->add_parent_child( volume, surf ); MB_CHK_ERR(rval);

  // set sense tag
  rval = GTT->set_sense(surf, volume, SENSE_FORWARD); MB_CHK_ERR(rval);

  delete GTT;

  return MB_SUCCESS;
}

int main()
{
  int result = 0;

  result += RUN_TEST(test_two_cubes);
  result += RUN_TEST(test_three_cubes);
  result += RUN_TEST(test_four_cubes);

  return result;

}

ErrorCode get_all_handles(Interface *mbi)
{
  ErrorCode rval;

  rval = mbi->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
                                name_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  MB_CHK_ERR(rval);

  rval = mbi->tag_get_handle( "OBJECT_NAME", 32, MB_TYPE_OPAQUE,
                               obj_name_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  MB_CHK_ERR(rval);

  int negone = -1;
  rval = mbi->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
                                geom_tag, MB_TAG_SPARSE|MB_TAG_CREAT,&negone);
  MB_CHK_ERR(rval);

  id_tag = mbi->globalId_tag();

  rval = mbi->tag_get_handle( CATEGORY_TAG_NAME,
                              CATEGORY_TAG_SIZE,
                              MB_TYPE_OPAQUE,
                              category_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );

  MB_CHK_ERR(rval);
  return MB_SUCCESS;

}

/* This function tests that the tree built by generate_hierarchy is the same
   as the reference tree
*/
bool check_tree ( Interface *mbi, GeomTopoTool *GTT, std::map< int, std::set<int> > &ref_map )
{
  ErrorCode rval;
  int vol_id;
  std::set<int> test_set;

  Range vols;
  rval = GTT->get_gsets_by_dimension(3, vols);
  if (ref_map.size() != vols.size() )
   {
    return false;
   }

  //go through volumes, create sets of children
  for ( Range::iterator it = vols.begin(); it != vols.end() ; ++it)
    {
      //get vol id
      rval = mbi->tag_get_data(id_tag, &(*it), 1, &vol_id ); MB_CHK_ERR(rval);

      //check if test vol in ref map
      if (ref_map.find(vol_id) == ref_map.end())
        {
          return false;
        }

      //put range of child surfaces into set
      Range child_surfs;
      test_set.clear();
      child_surfs = get_children_by_dimension( mbi, *it, 2);

      for (Range::iterator j = child_surfs.begin() ; j != child_surfs.end() ; ++j )
        {
            int child_id;

            rval = mbi->tag_get_data(id_tag, &(*j), 1, &child_id ); MB_CHK_ERR(rval);
            test_set.insert(child_id);
        }

      // compare sets
      if ( test_set != ref_map[vol_id] )
        {
         return false;
        }
    }

  return true;
}

Range get_children_by_dimension(Interface *mbi, EntityHandle parent, int desired_dimension)
{
  ErrorCode rval;
  Range all_children, desired_children;
  Range::iterator it;
  int actual_dimension;

  all_children.clear();
  rval = mbi->get_child_meshsets(parent, all_children);
  MB_CHK_SET_ERR_RET_VAL(rval, "Failed to get child meshsets", all_children);

  for ( it = all_children.begin() ; it != all_children.end() ; ++it)
    {
      rval = mbi->tag_get_data(geom_tag, &(*it), 1, &actual_dimension);
	  MB_CHK_SET_ERR_RET_VAL(rval, "Failed to get geom tag from child meshset", all_children);
      if ( actual_dimension == desired_dimension )
        {
          desired_children.insert(*it);
        }
    }

  return desired_children;

}

/* This function contains info for the scale and translation vectors of
   four different cubes that will be used in the hierarchy testing
*/
void get_cube_info( int cube_id, std::vector<double> &scale, std::vector<double> &trans )
{
  scale.clear();
  trans.clear();

  if ( cube_id == 1 )
    {
      scale.push_back(1);
      scale.push_back(1);
      scale.push_back(1);
      trans.push_back(0);
      trans.push_back(0);
      trans.push_back(0);
    }
  if ( cube_id == 2 )
    {
      scale.push_back(4);
      scale.push_back(4);
      scale.push_back(4);
      trans.push_back(0);
      trans.push_back(0);
      trans.push_back(0);
    }
  if ( cube_id == 3 )
    {
      scale.push_back(8);
      scale.push_back(8);
      scale.push_back(8);
      trans.push_back(0);
      trans.push_back(0);
      trans.push_back(0);
    }
  if ( cube_id == 4 )
    {
      scale.push_back(40);
      scale.push_back(40);
      scale.push_back(40);
      trans.push_back(0);
      trans.push_back(0);
      trans.push_back(0);
    }

}

/*
 * One large cube that contains a smaller one.
 */
void test_two_cubes()
{
  ErrorCode rval;

  Interface *mbi = new Core();

  // get all handles (dimension, id, sense)
  rval = get_all_handles(mbi);
  MB_CHK_ERR_RET(rval);

  int len = 2;
  int num[2] = {1, 2};

  //build reference map
  std::map< int, std::set<int> > ref_map;
  ref_map[1].insert(1);
  ref_map[2].insert(2);
  ref_map[2].insert(1);

  heappermute(mbi, num, len, ref_map, len);

  delete mbi;

}

/*
 * One large cube that contains two others.
 * The two inner cubes are siblings.
 */
void test_three_cubes()
{
  ErrorCode rval;

  Interface *mbi = new Core();

  // get all handles (dimension, id, sense)
  rval = get_all_handles(mbi);
  MB_CHK_ERR_RET(rval);

  int len = 3;
  int num[3] = {1, 2, 3};

  //build reference map
  std::map< int, std::set<int> > ref_map;
  ref_map[1].insert(1);
  ref_map[2].insert(2);
  ref_map[2].insert(1);
  ref_map[3].insert(3);
  ref_map[3].insert(2);

  heappermute(mbi, num, len, ref_map, len);

  delete mbi;
}

/*
 * Four nested cubes of decreasing size; one placed inside the next
 */
void test_four_cubes()
{
  ErrorCode rval;

  Interface *mbi = new Core();

  // get all handles (dimension, id, sense)
  rval = get_all_handles(mbi);
  MB_CHK_ERR_RET(rval);

  int len = 4;
  int num[4] = {1, 2, 3, 4};

  //build reference map
  std::map< int, std::set<int> > ref_map;
  ref_map[1].insert(1);
  ref_map[2].insert(2);
  ref_map[2].insert(1);
  ref_map[3].insert(3);
  ref_map[3].insert(2);
  ref_map[4].insert(4);
  ref_map[4].insert(3);

  heappermute(mbi, num, len, ref_map, len);

  delete mbi;
}


/* Heap's algorithm generates all possible permutations of n objects
   This function is a modification of code found here:
   http://www.sanfoundry.com/c-program-implement-heap-algorithm-permutation-n-numbers
*/
void heappermute(Interface *mbi, int v[], int n, std::map< int, std::set<int> > ref_map, int len)
{

  ErrorCode rval;
  std::vector<double> scale, trans;
  EntityHandle vol;
  Range flat_vols;

  if (n == 1)
    {
      //build cubes
      flat_vols.clear();
      for (int i = 0; i < len; i++)
        {
          get_cube_info(v[i], scale, trans);
          build_cube(mbi, scale, trans, v[i], vol);
		  flat_vols.insert(vol);

        }

	  // construct the topology
	  GeomTopoTool *GTT = new GeomTopoTool(mbi);
      // first build obbs-- necessary for topology construction
      rval = GTT->construct_obb_trees();
      MB_CHK_ERR_RET(rval);
	  rval = GTT->restore_topology_from_geometric_inclusion(flat_vols);
      MB_CHK_ERR_RET(rval);

      //test the topology
      bool result = check_tree(mbi, GTT, ref_map );
      CHECK_EQUAL(1, result);
      delete GTT;

      // delete the geometry so new one can be built;
      rval = mbi->delete_mesh();
      MB_CHK_ERR_RET(rval);

    }

  else
    {
      for (int i = 0; i < n; i++)
        {
          heappermute(mbi, v, n-1, ref_map, len);
          if (n % 2 == 1)
            {
              swap(&v[0], &v[n-1]);
            }

          else
            {
              swap(&v[i], &v[n-1]);
            }
        }
    }

}

void swap(int *x, int *y)
{
  int temp;

  temp = *x;
  *x = *y;
  *y = temp;
}

