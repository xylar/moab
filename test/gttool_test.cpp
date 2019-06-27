/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/**
 * \file gttool_test.cpp
 *
 * \brief test geometrize method from GeomTopoTool
 *
 */
#include "moab/Core.hpp"
#include <iostream>

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "TestUtil.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTagConventions.hpp"

using namespace moab;

std::string filename;
std::string filename2;
std::string ofile;
std::string ofile2;
std::string ofile3;
std::string ofile4;
std::string ofile5;

const char OBB_ROOT_TAG_NAME[] = "OBB_ROOT";
Tag obbRootTag;

bool remove_output_file;
ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet);

ErrorCode create_shell_test(Interface * mb);

ErrorCode duplicate_model_test(Interface * mb);

ErrorCode check_model_test(Interface * mb);

ErrorCode test_root_sets_resize(Interface *mb);

ErrorCode test_delete_obb_tree(Interface *mb);

ErrorCode test_restore_obb_trees(Interface *mb, Interface *mb2, Interface *mb3);

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
    std::cout << "Success";
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}

int main(int argc, char *argv[])
{
  filename = TestDir + "/partBed.smf";
  filename2 = TestDir + "/test_geom.h5m";
  ofile = "output.h5m";
  ofile2 = "shell.h5m";
  ofile3 = "shellCopy.h5m";
  ofile4 = "geom_w_obbs.h5m";
  ofile5 = "geom_missing_obb.h5m";

  remove_output_file = true;
  bool only_check = false;
  bool only_geometrize = false;


  if (argc == 1) {
    std::cout << "Using default input file and output files " << filename << " " << ofile <<
        " " << ofile2 <<  " " << ofile3 << std::endl;
  } else if (argc == 5) {
    filename = argv[1];
    ofile = argv[2];
    ofile2 = argv[3];
    ofile3 = argv[4];
    remove_output_file = false;
  } else if (argc==2) {
    ofile3 = argv[1];// check model only from file
    only_check = true;
    remove_output_file = false; // this is input now, do not delete it
  } else if (argc == 3) {
    filename = argv[1];
    ofile = argv[2];
    only_geometrize = true;
    remove_output_file = false;
  }
    else
  {
    std::cerr << "Usage: " << argv[0] << " [surface_mesh] [mbgeo_file] [shellfile] [copyshellfile] " << std::endl;
    return 1;
  }

  int number_tests_successful = 0;
  int number_tests_failed = 0;

  Core mbcore;
  Interface * mb = &mbcore;

  // use this as a tool to check a model
  if (only_check)
  {
    if (MB_SUCCESS==check_model_test( mb))
      std::cout << ofile3 << " passed gtt check\n";
    return 0;
  }

  mb->load_file(filename.c_str());

  //   FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  std::cout << "geometrize test: ";
  ErrorCode rval = geometrize_test( mb, 0); // just pass the root set
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  if (only_geometrize)
  {
    return number_tests_failed;
  }
  std::cout << "create shell test: ";
  rval = create_shell_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "duplicate model test: ";
  rval = duplicate_model_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "check_model test: ";
  rval = check_model_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "test_rootsets_resize: ";
  Interface* mb2 = new Core();
  rval = test_root_sets_resize(mb2);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  delete mb2;
  std::cout << "\n";

  std::cout << "test_delete_obb_tree: ";
  Interface* mb3 = new Core();
  rval = test_delete_obb_tree(mb3);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  delete mb3;
  std::cout << "\n";

  std::cout << "test_restore_obb_trees: ";
  Interface* mb4 = new Core();
  Interface* mb5 = new Core();
  Interface* mb6 = new Core();
  rval = test_restore_obb_trees(mb4, mb5, mb6);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  delete mb4;
  delete mb5;
  delete mb6;
  std::cout << "\n";

  return number_tests_failed;
}
ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet)
{
  GeomTopoTool gtt(mb);
  EntityHandle outSet;
  ErrorCode rval=gtt.geometrize_surface_set(inputSet, outSet);MB_CHK_SET_ERR(rval, "Can't geometrize the set\n");

  std::cout<<"writing output file: " << ofile.c_str() << " ";
  rval=mb->write_file(ofile.c_str(), 0, 0, &outSet, 1);MB_CHK_SET_ERR(rval, "Can't write output file\n");
  if (remove_output_file)
  {
    remove(ofile.c_str());
  }
  return MB_SUCCESS;
}

ErrorCode create_shell_test(Interface * mb)
{
  // we should be able to delete mesh and create a model from scratch

  ErrorCode rval = mb->delete_mesh();MB_CHK_SET_ERR(rval, "Can't delete existing mesh\n");

  // create some vertices
  double coords [] = { 0, 0, 0,
                     1, 0, 0.1,
                     2, 0, 0,
                     3, 0, -0.1,
                     0, 1, 0,
                     1, 1, 0,
                     2, 1, 0,
                     3, 1, -0.1,
                     0, 2, 0,
                     1, 2, -0.1,
                     2, 2, -0.1,
                     3, 2, -0.2,
                     0, 0, 1,
                     1, 0, 0.9,
                     2, 0.1, 0.85,
                     3, 0.2, 0.8,
                     0, 0.1, 2,
                     1, 0.1, 2,
                     2.1, 0.2, 2.1,
                     3.1, 0.2, 2.1 };

  int nvert = 20;
  Range verts;
  rval = mb->create_vertices(coords, nvert, verts);MB_CHK_SET_ERR(rval, "Can't create vertices\n");

  EntityHandle connec [] = { 1, 2, 5,
                    5, 2, 6,
                    2, 3, 6,
                    6, 3, 7,
                    3, 4, 7,
                    7, 4, 8,
                    5, 6, 9,
                    9, 6, 10,
                    6, 7, 10,
                    10, 7, 11,
                    7, 8, 11,
                    11, 8, 12,  // first face, 1-12
                    13, 14, 1,
                    1, 14, 2,
                    14, 15, 2,
                    2, 15, 3,
                    15, 16, 3,
                    3, 16, 4,
                    17, 18, 13,
                    13, 18, 14,
                    18, 19, 14,
                    14, 19, 15,
                    19, 20, 15,
                    15, 20, 16 // second face, 13-24
                   };
  EntityHandle elem;
  int nbTri = sizeof(connec)/3/sizeof(EntityHandle);
  int i = 0;
  std::vector<EntityHandle> tris;
  for (i=0; i<nbTri; i++)
  {
    mb->create_element(MBTRI, &connec[3*i], 3, elem);
    tris.push_back(elem);
  }


  // create some edges too
  EntityHandle edges [] = { 1, 2,
                          2, 3,
                          3, 4,  // geo edge 1  1:3
                          4, 8,
                          8, 12,  // geo 2         4:5
                          12, 11,
                          11, 10,
                          10, 9, // geo 3 6:8
                          9, 5,
                          5, 1,  // geo 4  9:10
                          1, 13,
                          13, 17, // geo 5  11:12
                          17, 18,
                          18, 19,
                          19, 20, // geo 6  13:15
                          20, 16,
                          16, 4   // geo 7  16:17
                           };
  int nbEdges = sizeof(edges)/2/sizeof(EntityHandle);
  std::vector<EntityHandle> edgs;
  for (i=0; i<nbEdges; i++)
  {
    mb->create_element(MBEDGE, &edges[2*i], 2, elem);
    edgs.push_back(elem);
  }
  // create some sets, and create some ordered sets for edges
  EntityHandle face1, face2;
  rval = mb->create_meshset(MESHSET_SET, face1);MB_CHK_ERR(rval);
  rval = mb->add_entities(face1, &tris[0], 12);MB_CHK_ERR(rval);

  rval = mb->create_meshset(MESHSET_SET, face2);MB_CHK_ERR(rval);
  // next 12 triangles
  rval = mb->add_entities(face2, &tris[12], 12);MB_CHK_ERR(rval);

  // the orientation and senses need to be set for face edges
  moab::GeomTopoTool gTopoTool(mb, false);

  rval = gTopoTool.add_geo_set(face1, 2);MB_CHK_ERR(rval);

  rval = gTopoTool.add_geo_set(face2, 2);MB_CHK_ERR(rval);

  // create some edges
  EntityHandle edge[7]; //edge[0] has EH 1...
 ;
  for (i=0; i<7; i++)
  {
    rval = mb->create_meshset(MESHSET_ORDERED, edge[i]);MB_CHK_ERR(rval);
    rval = gTopoTool.add_geo_set(edge[i], 1);MB_CHK_ERR(rval);
  }

  // first 3 mesh edges...
  rval = mb->add_entities(edge[0], &edgs[0], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[1], &edgs[3], 2);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[2], &edgs[5], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[3], &edgs[8], 2);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[4], &edgs[10], 2);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[5], &edgs[12], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[6], &edgs[15], 2);MB_CHK_ERR(rval);

  // create some sets for vertices; also need to create some for parent/child relationships
  EntityHandle vertSets[6];// start from 0

  for (i=0; i<6; i++)
  {
    rval = mb->create_meshset(MESHSET_SET, vertSets[i]);MB_CHK_ERR(rval);
    rval = gTopoTool.add_geo_set(vertSets[i], 0);MB_CHK_ERR(rval);
  }

  EntityHandle v(1); // first vertex;
  rval = mb->add_entities(vertSets[0], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (4);
  rval = mb->add_entities(vertSets[1], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (9);
  rval = mb->add_entities(vertSets[2], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (12);
  rval = mb->add_entities(vertSets[3], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (17);
  rval = mb->add_entities(vertSets[4], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (20);
  rval = mb->add_entities(vertSets[5], &v, 1);MB_CHK_ERR(rval);

  // need to add parent-child relations between sets
  // edge 1 : 1-2
  rval = mb ->add_parent_child( edge[0], vertSets[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[0], vertSets[1]);MB_CHK_ERR(rval);
  // edge 2 : 2-4
  rval = mb ->add_parent_child( edge[1], vertSets[1]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[1], vertSets[3]);MB_CHK_ERR(rval);

  // edge 3 : 4-3
  rval = mb ->add_parent_child( edge[2], vertSets[3]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[2], vertSets[2]);MB_CHK_ERR(rval);

  // edge 4 : 4-1
  rval = mb ->add_parent_child( edge[3], vertSets[2]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[3], vertSets[0]);MB_CHK_ERR(rval);

  // edge 5 : 1-5
  rval = mb ->add_parent_child( edge[4], vertSets[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[4], vertSets[4]);MB_CHK_ERR(rval);

  // edge 6 : 5-6
  rval = mb ->add_parent_child( edge[5], vertSets[4]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[5], vertSets[5]);MB_CHK_ERR(rval);

  // edge 7 : 6-2
  rval = mb ->add_parent_child( edge[6], vertSets[5]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[6], vertSets[1]);MB_CHK_ERR(rval);

  // face 1: edges 1, 2, 3, 4
  rval = mb ->add_parent_child( face1, edge[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[1]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[2]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[3]);MB_CHK_ERR(rval);

  // face 2: edges 1, 5, 6, 7
  rval = mb ->add_parent_child( face2, edge[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[4]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[5]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[6]);MB_CHK_ERR(rval);

  // set senses !!
  std::vector<EntityHandle> faces;
  faces.push_back(face1); // the face1 has all edges oriented positively
  std::vector<int> senses;
  senses.push_back(moab::SENSE_FORWARD); //

  //faces.push_back(face1);
  //faces.push_back(face2);
  gTopoTool.set_senses(edge[1], faces, senses);
  gTopoTool.set_senses(edge[2], faces, senses);
  gTopoTool.set_senses(edge[3], faces, senses);

  faces[0]=face2; // set senses for edges for face2
  gTopoTool.set_senses(edge[4], faces, senses);
  gTopoTool.set_senses(edge[5], faces, senses);
  gTopoTool.set_senses(edge[6], faces, senses);

  // the only complication is edge1 (edge[0]), that has face1 forward and face 2 reverse
  faces[0]=face1;
  faces.push_back(face2);
  senses.push_back(moab::SENSE_REVERSE); // -1 is reverse; face 2 is reverse for edge1 (0)
  // forward == 0, reverse ==1
  gTopoTool.set_senses(edge[0], faces, senses);

  rval = mb->write_mesh(ofile2.c_str());MB_CHK_ERR(rval);

  rval = mb->delete_mesh();MB_CHK_ERR(rval);

  // now test loading it up
  rval = mb->load_file(ofile2.c_str());MB_CHK_ERR(rval);

  if (remove_output_file)
  {
    remove(ofile2.c_str());
  }
  // do some tests on geometry

  // it would be good to have a method on updating the geom topo tool
  // so we do not have to create another one
  moab::GeomTopoTool gTopoTool2(mb, true);// to find the geomsets
  Range ranges[5];
  rval =  gTopoTool2.find_geomsets(ranges);

  assert(MB_SUCCESS==rval);
  assert(ranges[0].size()==6);
  assert(ranges[1].size()==7);
  assert(ranges[2].size()==2);
  assert(ranges[3].size()==0);
  assert(ranges[4].size()==0);

  return MB_SUCCESS;
}

ErrorCode duplicate_model_test(Interface * mb)
{
  moab::GeomTopoTool gTopoTool2(mb, true);// to find the geomsets

  GeomTopoTool * newModel = NULL;
  ErrorCode rval = gTopoTool2.duplicate_model(newModel);
  if (NULL == newModel || rval!=MB_SUCCESS)
    return MB_FAILURE;

  Range ranges[5];
  rval = newModel->find_geomsets(ranges);MB_CHK_ERR(rval);

  assert(ranges[0].size()==6);
  assert(ranges[1].size()==7);
  assert(ranges[2].size()==2);
  assert(ranges[3].size()==0);

  // write the model to a test file
  EntityHandle rootModelSet = newModel->get_root_model_set();
  std::cout<<"writing duplicated model file: " << ofile3.c_str() << " ";
  rval=mb->write_file(ofile3.c_str(), 0, 0, &rootModelSet, 1);MB_CHK_SET_ERR(rval, "Can't write output files\n");

  delete newModel; // we are done with the new geom topo tool
  // do not delete yet the output file, delay after the next test
  /*if (remove_output_file)
  {
    remove(ofile3.c_str());
  }*/

  return MB_SUCCESS;
}

ErrorCode check_model_test(Interface * mb)
{
  ErrorCode rval = mb->delete_mesh();MB_CHK_SET_ERR(rval, "Can't delete existing mesh\n");

  rval = mb->load_file(ofile3.c_str());MB_CHK_ERR(rval);

  // do some tests on geometry
  // it would be good to have a method on updating the geom topo tool
  // so we do not have to create another one
  if (remove_output_file)
  {
    remove(ofile3.c_str());
  }
  moab::GeomTopoTool gTopoTool(mb, true);

  if (!gTopoTool.check_model())
    return MB_FAILURE;

  return MB_SUCCESS;
}

ErrorCode test_root_sets_resize(Interface *mb) {

  // load the test file
  ErrorCode rval = mb->load_file(filename2.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load input file");

  // create a GTT with all default settings
  moab::GeomTopoTool* gTopoTool = new GeomTopoTool(mb);


  Tag geomTag;

  rval = mb->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1,
					   MB_TYPE_INTEGER, geomTag, MB_TAG_CREAT|MB_TAG_SPARSE);
  MB_CHK_SET_ERR(rval, "Error: Failed to create geometry dimension tag");

  Range surfs;

  const int dim = 2;
  const void* const dim_val[] = { &dim };
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &geomTag,
					       dim_val, 1, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get entity sets by type and tag");

  // in reverse order, add surfaces and construct their trees
  for (Range::reverse_iterator rit = surfs.rbegin(); rit != surfs.rend(); rit++ ) {

    rval = gTopoTool->add_geo_set(*rit, 2);
    MB_CHK_SET_ERR(rval, "Failed to add geometry set to GTT");

    rval = gTopoTool->construct_obb_tree(*rit);
    MB_CHK_SET_ERR(rval, "Failed to construct obb tree for surf " << *rit);

  }

  for(Range::iterator it = surfs.begin(); it != surfs.end(); it++ ) {
    EntityHandle obb_root_set;
    rval = gTopoTool->get_root(*it, obb_root_set);
    MB_CHK_SET_ERR(rval, "Failed to get obb tree root from GTT");

    // make sure the returned root is valid
    CHECK(obb_root_set);
  }

  // clean up GTT
  delete gTopoTool;

  // create a GTT with all default settings
  gTopoTool = new moab::GeomTopoTool(mb, false, 0, false);


  // in reverse order, add surfaces and construct their trees
  for (Range::reverse_iterator rit = surfs.rbegin(); rit != surfs.rend(); rit++ ) {

    rval = gTopoTool->add_geo_set(*rit, 2);
    MB_CHK_SET_ERR(rval, "Failed to add geometry set to GTT");

    rval = gTopoTool->construct_obb_tree(*rit);
    MB_CHK_SET_ERR(rval, "Failed to construct obb tree for surf " << *rit);

  }

  for(Range::iterator it = surfs.begin(); it != surfs.end(); it++ ) {
    EntityHandle obb_root_set;
    rval = gTopoTool->get_root(*it, obb_root_set);
    MB_CHK_SET_ERR(rval, "Failed to get obb tree root from GTT");

    // make sure the returned root is valid
    CHECK(obb_root_set);
  }

  delete gTopoTool;

  rval = mb->delete_mesh();
  MB_CHK_SET_ERR(rval, "Failed to delete mesh in MOAB instance.");

  return MB_SUCCESS;

}
				
ErrorCode test_delete_obb_tree(Interface *mb){

  // Load the test file
  ErrorCode rval = mb->load_file(filename2.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load input file");

  // Create a GTT with all default settings
  moab::GeomTopoTool* gTopoTool = new GeomTopoTool(mb);

  // Get all volumes and surfaces
  Range vols, surfs;
  rval = gTopoTool->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  rval = gTopoTool->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get surface gsets");

  // Build obb tree for volume
  EntityHandle test_vol = vols.front();
  rval = gTopoTool->construct_obb_tree(test_vol);
  MB_CHK_SET_ERR(rval, "Error constructing all trees.");


  // Get the obbRootTag for vol
  rval = mb->tag_get_handle(OBB_ROOT_TAG_NAME, 1,
                            MB_TYPE_HANDLE, obbRootTag,
                            MB_TAG_CREAT|MB_TAG_SPARSE);
  MB_CHK_SET_ERR_CONT(rval, "Error: Failed to create obb root tag");
  EntityHandle gbroot;
  rval = mb->tag_get_data(obbRootTag, &test_vol, 1, &gbroot);
  MB_CHK_SET_ERR(rval, "Failed to get the obb root tag");

  // Test if obb tree in ModelSet
  EntityHandle test_vol_root;
  rval = gTopoTool->get_root(test_vol, test_vol_root);
  MB_CHK_SET_ERR(rval, "Obb root not in ModelSet");

  // CASE 1: Delete vol obb tree including all child surface trees
  rval = gTopoTool->delete_obb_tree(test_vol, false);
  MB_CHK_SET_ERR(rval, "Error deleting volume tree.");

  // Make sure vol tree is gone
  EntityHandle newroot;
  rval = mb->tag_get_data(obbRootTag, &test_vol, 1, &newroot);
  if (MB_SUCCESS == rval){
    return MB_FAILURE;
  }

  // Make sure its child surf trees also gone
  Range::iterator surf_it;
  for(surf_it = surfs.begin(); surf_it != surfs.end(); ++surf_it){
    EntityHandle test_surf_root_gone;
    rval = mb->tag_get_data(obbRootTag, &(*surf_it), 1, &test_surf_root_gone);
    if (MB_SUCCESS == rval){
      return MB_FAILURE;
    }
  }

  // Rebuild vol tree
  rval = gTopoTool->construct_obb_tree(test_vol);
  MB_CHK_SET_ERR(rval, "Error constructing all trees.");

  // CASE 2: Delete just vol, not surf trees
  rval = gTopoTool->delete_obb_tree(test_vol, true);
  MB_CHK_SET_ERR(rval, "Error deleting volume tree.");

  // Make sure vol tree is gone
  rval = mb->tag_get_data(obbRootTag, &test_vol, 1, &gbroot);
  if (MB_SUCCESS == rval){
    return MB_FAILURE;
  }

  // Make sure its child surf trees remain
  for(surf_it = surfs.begin(); surf_it != surfs.end(); ++surf_it){
    EntityHandle test_surf_root;
    rval = mb->tag_get_data(obbRootTag, &(*surf_it), 1, &test_surf_root);
    MB_CHK_SET_ERR(rval, "Problem getting obb root of surface.");
  }

  // CASE 3: Delete surf tree
  EntityHandle test_surf = surfs.front();
  rval = gTopoTool->delete_obb_tree(test_surf, false);
  MB_CHK_SET_ERR(rval, "Error deleting surface tree.");

  // Make sure surf tree is gone
  rval = mb->tag_get_data(obbRootTag, &test_surf, 1, &gbroot);
  if (MB_SUCCESS == rval){
    return MB_FAILURE;
  }

  delete gTopoTool;

  return MB_SUCCESS;
}

ErrorCode test_restore_obb_trees(Interface *mb, Interface *mb2, Interface *mb3){

  // Load the test file
  ErrorCode rval = mb->load_file(filename2.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load input file");

  // Create a GTT with all default settings
  moab::GeomTopoTool* gTopoTool = new GeomTopoTool(mb);

  // Build all obb trees
  rval = gTopoTool->construct_obb_trees();
  MB_CHK_SET_ERR(rval, "Error constructing all trees");
  // Write the file with all obbs
  rval=mb->write_file(ofile4.c_str());
  MB_CHK_SET_ERR(rval, "Can't write output file");

  // Delete a vol obb tree
  Range vols;
  rval = gTopoTool->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  EntityHandle test_vol = vols.front();
  rval = gTopoTool->delete_obb_tree(test_vol, false);
  MB_CHK_SET_ERR(rval, "Error deleting volume tree");
  // Write the file missing an obb
  rval=mb->write_file(ofile5.c_str());
  MB_CHK_SET_ERR(rval, "Can't write output file");

  // Load file containing obbs
  rval = mb2->load_file(ofile4.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load file containing obbs");

  // 1) Check that roots are NOT restored by default GTT settings
  // GeomTopoTool(Interface *impl, bool find_geoments = false, EntityHandle modelRootSet = 0,
  //             bool p_rootSets_vector = true, bool restore_rootSets = true);
  moab::GeomTopoTool* gTopoTool2 = new GeomTopoTool(mb2, false, 0, true, true);

  vols.clear();
  Range surfs;
  rval = gTopoTool2->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  rval = gTopoTool2->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");

  Range gsets;
  gsets.insert_list(surfs.begin(), surfs.end());
  gsets.insert_list(vols.begin(), vols.end());
  EntityHandle test_root2;
  for (Range::iterator rit = gsets.begin(); rit != gsets.end(); ++rit) {
    rval = gTopoTool2->get_root(*rit, test_root2);
    if (MB_SUCCESS == rval){
      return MB_FAILURE;
    }
  }

  // 2) Check that roots ARE restored by setting find_geoments and restore_rootSets to true
  moab::GeomTopoTool* gTopoTool3 = new GeomTopoTool(mb2, true, 0, true, true);

  vols.clear();
  rval = gTopoTool3->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  surfs.clear();
  rval = gTopoTool3->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");

  gsets.clear();
  gsets.insert_list(surfs.begin(), surfs.end());
  gsets.insert_list(vols.begin(), vols.end());
  EntityHandle test_root3;
  for (Range::iterator rit = gsets.begin(); rit != gsets.end(); ++rit) {
    rval = gTopoTool3->get_root(*rit, test_root3);
    MB_CHK_SET_ERR(rval, "Failed to get obb tree root from GTT");
    CHECK(test_root3);
  }

  // 3) Check that roots are deleted and then rebuilt if an obb tree is missing

  // Load file missing obb
  rval = mb3->load_file(ofile5.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load file containing obbs");

  // Create GTT and try to restore OBBs
  moab::GeomTopoTool* gTopoTool4 = new GeomTopoTool(mb3, true, 0, true, true);

  // Check that roots still exist
  vols.clear();
  rval = gTopoTool4->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  surfs.clear();
  rval = gTopoTool4->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");

  gsets.clear();
  gsets.insert_list(surfs.begin(), surfs.end());
  gsets.insert_list(vols.begin(), vols.end());
  EntityHandle test_root4;
  for (Range::iterator rit = gsets.begin(); rit != gsets.end(); ++rit) {
    rval = gTopoTool4->get_root(*rit, test_root4);
    MB_CHK_SET_ERR(rval, "Failed to get obb tree root from GTT");
    CHECK(test_root4);
  }

  // 4) Check that roots exist but rootSets is NOT repopulated when find_geoments = true and restore_rootSets = false
  moab::GeomTopoTool* gTopoTool5 = new GeomTopoTool(mb2, true, 0, true, false);

  // Get the obbRootTag for vol
  rval = mb2->tag_get_handle(OBB_ROOT_TAG_NAME, 1,
                            MB_TYPE_HANDLE, obbRootTag,
                            MB_TAG_CREAT|MB_TAG_SPARSE);
  MB_CHK_SET_ERR_CONT(rval, "Error: Failed to create obb root tag");

  vols.clear();
  rval = gTopoTool5->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");
  surfs.clear();
  rval = gTopoTool5->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Failed to get volume gsets");

  gsets.clear();
  gsets.insert_list(surfs.begin(), surfs.end());
  gsets.insert_list(vols.begin(), vols.end());
  EntityHandle test_root5, tagged_root;
  for (Range::iterator rit = gsets.begin(); rit != gsets.end(); ++rit) {
    // Check that root still exits, but not in rootSet
    rval = mb2->tag_get_data(obbRootTag, &(*rit), 1, &tagged_root);
    MB_CHK_SET_ERR(rval, "Failed to get root from tag");
    CHECK(tagged_root);
    rval = gTopoTool5->get_root(*rit, test_root5);
    if (MB_SUCCESS == rval) return MB_FAILURE;
  }

  delete gTopoTool;
  delete gTopoTool2;
  delete gTopoTool3;
  delete gTopoTool4;
  delete gTopoTool5;

  return MB_SUCCESS;
}
