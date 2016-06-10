#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/* Input test file: rtttest.rtt
 */
#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/rtttest.rtt";
#else
static const char example[] = "rtttest.rtt";
#endif
// rtttest2.rtt
#ifdef MESHDIR
static const char example_2[] = STRINGIFY(MESHDIR) "/io/rtttest2.rtt";
#else
static const char example_2[] = "rtttest2.rtt";
#endif

void test_loadfile();
void test_meshset_tags();
void test_tets();
void test_tet_tags();
void test_triangles();
void test_triangles_tags();
void test_vertices();
void test_loadfile_2();
void test_meshset_tags_2();
void test_tets_2();
void test_tet_tags_2();
void test_triangles_2();
void test_triangles_tags_2();
void test_vertices_2();

void read_file( Interface& moab, const char* input_file );

int main() {
  int result = 0;
  result += RUN_TEST(test_loadfile);
  result += RUN_TEST(test_meshset_tags);
  result += RUN_TEST(test_tets);
  result += RUN_TEST(test_tet_tags);
  result += RUN_TEST(test_triangles);
  result += RUN_TEST(test_vertices);
  // second model tests
  result += RUN_TEST(test_loadfile_2);
  result += RUN_TEST(test_meshset_tags_2);
  result += RUN_TEST(test_tets_2);
  result += RUN_TEST(test_tet_tags_2);
  result += RUN_TEST(test_triangles_2);
  result += RUN_TEST(test_vertices_2);
  return result;
}

void read_file( Interface& moab, const char* input_file ) {
  ErrorCode rval;
  rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}

void test_loadfile() {
  Core moab;
  read_file( moab, example );
}

void test_meshset_tags() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBENTITYSET,entities);
  CHECK_ERR(rval);

  Tag id_tag;
  // get the tag handle
  rval = moab.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
			     id_tag, MB_TAG_DENSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&id_tag,0,1,entities);
  CHECK_ERR(rval);
  
  // each tet should have the material tag
  int num_vols = 10;
  int num_surfaces = 129;
  int num_sets = entities.size();
  CHECK_EQUAL(num_sets, num_vols+num_surfaces+1+1);

  // from the entities, get the volume meshsets, get it from the dim tag
  Tag dim_tag;
  entities.clear();
  // get the tag handle
  rval = moab.tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
			     dim_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged with dim 3
  int dim3 = 3;
  const void* const tag_vals_dim3[] = { &dim3 };
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&dim_tag,tag_vals_dim3,1,entities);
  CHECK_ERR(rval);

  // there should only be 10 meshsets of dimension 3
  num_sets = entities.size();
  CHECK_EQUAL(num_vols,num_sets);

  entities.clear();
  // get the entities that are tagged with dim 2
  int dim2 = 2;
  const void* const tag_vals_dim2[] = { &dim2 };
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&dim_tag,tag_vals_dim2,1,entities);
  CHECK_ERR(rval);

  // there should only be 129 meshsets of dimension 2
  num_sets = entities.size();
  CHECK_EQUAL(num_surfaces,num_sets);
}

void test_tets() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  // cells = 26710 - number of tets
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTET,entities);
  CHECK_ERR(rval);
  int num_tets = 26710;
  int num_tet_in_moab = entities.size();
  CHECK_EQUAL(num_tet_in_moab, num_tets);
}

void test_tet_tags() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTET,entities);
  CHECK_ERR(rval);

  int num_tets = 26710;
  int num_tet_in_moab = entities.size();
  CHECK_EQUAL(num_tet_in_moab, num_tets);
  
  // get the number of tets tagged with
  entities.clear();
  Tag material_number;
  // get the tag handle
  rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER,
			      material_number, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTET,&material_number,0,1,entities);
  // each tet should have the material tag
  int num_tet_tag = entities.size();
  CHECK_EQUAL(num_tet_tag, num_tets);
  CHECK_ERR(rval);
}

void test_triangles() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTRI,entities);
  CHECK_ERR(rval);

  int num_tri = 6383; // num tris annotated in rtttest.rtt
  int num_tri_in_moab = entities.size();
  CHECK_EQUAL(num_tri_in_moab, num_tri);
}

void test_triangles_tags() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTRI,entities);
  CHECK_ERR(rval);

  int num_tri = 6383; // num tris annotated in rtttest.rtt
  int num_tri_in_moab = entities.size();
  CHECK_EQUAL(num_tri_in_moab, num_tri);

  // get the number of tris tagged with SURFACE_NUMBER
  entities.clear();
  Tag surface_number;
  // get the tag handle
  rval = moab.tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER,
			      surface_number, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTRI,&surface_number,0,1,entities);
  // each tri should have the surface number tag
  int num_tri_tag = entities.size();
  CHECK_EQUAL(num_tri_tag, num_tri);
  CHECK_ERR(rval);

  // get the number of tris tagged with SIDEID_TAG
  entities.clear();
  Tag sideid_tag;
  // get the tag handle
  rval = moab.tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER,
			      sideid_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTRI,&sideid_tag,0,1,entities);
  // each tri should have the sideid tag
  num_tri_tag = entities.size();
  CHECK_EQUAL(num_tri_tag, num_tri);
  CHECK_ERR(rval);
}

void test_vertices() {
  Core moab;
  // load the data into moab
  read_file( moab, example );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBVERTEX,entities);
  CHECK_ERR(rval);

  int num_verts = 5397; // num verts annotated in rtttest.rtt
  int num_verts_in_moab = entities.size();
  CHECK_EQUAL(num_verts_in_moab, num_verts);
}

// This is a duplicate of the same tests as above but with a new rtt file
void test_loadfile_2() {
  Core moab;
  read_file( moab, example_2 );
}

void test_meshset_tags_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBENTITYSET,entities);
  CHECK_ERR(rval);

  Tag id_tag;
  // get the tag handle
  rval = moab.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
			     id_tag, MB_TAG_DENSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&id_tag,0,1,entities);
  CHECK_ERR(rval);
  
  // each tet should have the material tag
  int num_vols = 54;
  int num_surfaces = 11880;
  int num_sets = entities.size();
  CHECK_EQUAL(num_sets, num_vols+num_surfaces+1+1);

  // from the entities, get the volume meshsets, get it from the dim tag
  Tag dim_tag;
  entities.clear();
  // get the tag handle
  rval = moab.tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
			     dim_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged with dim 3
  int dim3 = 3;
  const void* const tag_vals_dim3[] = { &dim3 };
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&dim_tag,tag_vals_dim3,1,entities);
  CHECK_ERR(rval);

  // there should only be 10 meshsets of dimension 3
  num_sets = entities.size();
  CHECK_EQUAL(num_vols,num_sets);

  entities.clear();
  // get the entities that are tagged with dim 2
  int dim2 = 2;
  const void* const tag_vals_dim2[] = { &dim2 };
  rval = moab.get_entities_by_type_and_tag(0,moab::MBENTITYSET,&dim_tag,tag_vals_dim2,1,entities);
  CHECK_ERR(rval);

  // there should only be 129 meshsets of dimension 2
  num_sets = entities.size();
  CHECK_EQUAL(num_surfaces,num_sets);
}

void test_tets_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  // cells = 2617989 - number of tets
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTET,entities);
  CHECK_ERR(rval);
  int num_tets = 2617989;
  int num_tet_in_moab = entities.size();
  CHECK_EQUAL(num_tet_in_moab, num_tets);
}

void test_tet_tags_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTET,entities);
  CHECK_ERR(rval);

  int num_tets = 2617989;
  int num_tet_in_moab = entities.size();
  CHECK_EQUAL(num_tet_in_moab, num_tets);
  
  // get the number of tets tagged with
  entities.clear();
  Tag material_number;
  // get the tag handle
  rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER,
			      material_number, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTET,&material_number,0,1,entities);
  // each tet should have the material tag
  int num_tet_tag = entities.size();
  CHECK_EQUAL(num_tet_tag, num_tets);
  CHECK_ERR(rval);
}

void test_triangles_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTRI,entities);
  CHECK_ERR(rval);

  int num_tri = 581427; // num tris annotated in rtttest.rtt
  int num_tri_in_moab = entities.size();
  CHECK_EQUAL(num_tri_in_moab, num_tri);
}

void test_triangles_tags_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBTRI,entities);
  CHECK_ERR(rval);

  int num_tri = 581427; // num tris annotated in rtttest.rtt
  int num_tri_in_moab = entities.size();
  CHECK_EQUAL(num_tri_in_moab, num_tri);

  // get the number of tris tagged with SURFACE_NUMBER
  entities.clear();
  Tag surface_number;
  // get the tag handle
  rval = moab.tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER,
			      surface_number, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTRI,&surface_number,0,1,entities);
  // each tri should have the surface number tag
  int num_tri_tag = entities.size();
  CHECK_EQUAL(num_tri_tag, num_tri);
  CHECK_ERR(rval);

  // get the number of tris tagged with SIDEID_TAG
  entities.clear();
  Tag sideid_tag;
  // get the tag handle
  rval = moab.tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER,
			      sideid_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);

  // get the entities that are tagged
  rval = moab.get_entities_by_type_and_tag(0,moab::MBTRI,&sideid_tag,0,1,entities);
  // each tri should have the sideid tag
  num_tri_tag = entities.size();
  CHECK_EQUAL(num_tri_tag, num_tri);
  CHECK_ERR(rval);
}

void test_vertices_2() {
  Core moab;
  // load the data into moab
  read_file( moab, example_2 );
  // query the dataset to make sure that there are the correct number of cells
  Range entities;
  ErrorCode rval = moab.get_entities_by_type(0,moab::MBVERTEX,entities);
  CHECK_ERR(rval);

  int num_verts = 457318; // num verts annotated in rtttest.rtt
  int num_verts_in_moab = entities.size();
  CHECK_EQUAL(num_verts_in_moab, num_verts);
}


