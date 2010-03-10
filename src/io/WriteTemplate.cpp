/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */


#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "WriteTemplate.hpp"

#include <utility>
#include <algorithm>
#include <time.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "assert.h"
#include "MBInternals.hpp"
#include "ExoIIUtil.hpp"
#include "MBTagConventions.hpp"
#include "MBWriteUtilIface.hpp"

#define INS_ID(stringvar, prefix, id) \
          sprintf(stringvar, prefix, id)

MBWriterIface* WriteTEMPLATE::factory( MBInterface* iface )
  { return new WriteTEMPLATE( iface ); }

WriteTEMPLATE::WriteTEMPLATE(MBInterface *impl) 
    : mbImpl(impl), mCurrentMeshHandle(0)
{
  assert(impl != NULL);

  void* ptr = 0;
  impl->query_interface( "MBWriteUtilIface", &ptr );
  mWriteIface = reinterpret_cast<MBWriteUtilIface*>(ptr);

  // initialize in case tag_get_handle fails below
  //! get and cache predefined tag handles
  int dum_val = 0;
  MBErrorCode result = impl->tag_get_handle(MATERIAL_SET_TAG_NAME,  mMaterialSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mMaterialSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(DIRICHLET_SET_TAG_NAME, mDirichletSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(DIRICHLET_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mDirichletSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(NEUMANN_SET_TAG_NAME,   mNeumannSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(NEUMANN_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mNeumannSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = impl->tag_create(HAS_MID_NODES_TAG_NAME, 4*sizeof(int), MB_TAG_SPARSE, mHasMidNodesTag,
                              dum_val_array);
  }
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mGlobalIdTag,
                              &dum_val);
  
  dum_val = -1;
  result = impl->tag_get_handle("__matSetIdTag", mMatSetIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create("__matSetIdTag", sizeof(int), MB_TAG_DENSE, mMatSetIdTag,
                              &dum_val);
  

  impl->tag_create("WriteTEMPLATE element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

}

WriteTEMPLATE::~WriteTEMPLATE() 
{
  std::string iface_name = "MBWriteUtilIface";
  mbImpl->release_interface(iface_name, mWriteIface);

  mbImpl->tag_delete(mEntityMark);

}

void WriteTEMPLATE::reset_matset(std::vector<WriteTEMPLATE::MaterialSetData> &matset_info)
{
  std::vector<WriteTEMPLATE::MaterialSetData>::iterator iter;
  
  for (iter = matset_info.begin(); iter != matset_info.end(); iter++)
  {
    delete (*iter).elements;
  }
}

MBErrorCode WriteTEMPLATE::write_file(const char *file_name, 
                                      const bool /* overwrite (commented out to remove warning) */,
                                      const FileOptions& opts,
                                      const MBEntityHandle *ent_handles,
                                      const int num_sets,
                                      const std::vector<std::string>&,
                                      const MBTag* ,
                                      int ,
                                      int )
{
  assert(0 != mMaterialSetTag &&
         0 != mNeumannSetTag &&
         0 != mDirichletSetTag);

    // check the file name
  if (NULL == strstr(file_name, ".template"))
    return MB_FAILURE;

  std::vector<MBEntityHandle> matsets, dirsets, neusets, entities;

  fileName = file_name;
  
    // separate into material sets, dirichlet sets, neumann sets

  if (num_sets == 0) {
      // default to all defined sets
    MBRange this_range;
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mMaterialSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(matsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mDirichletSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(dirsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mNeumannSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(neusets));
  }
  else {
    int dummy;
    for (const MBEntityHandle *iter = ent_handles; iter < ent_handles+num_sets; iter++) 
    {
      if (MB_SUCCESS == mbImpl->tag_get_data(mMaterialSetTag, &(*iter), 1, &dummy))
        matsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mDirichletSetTag, &(*iter), 1, &dummy))
        dirsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mNeumannSetTag, &(*iter), 1, &dummy))
        neusets.push_back(*iter);
    }
  }
  
    // if there is nothing to write just return.
  if (matsets.empty() && dirsets.empty() && neusets.empty())
    return MB_FILE_WRITE_ERROR;

  std::vector<WriteTEMPLATE::MaterialSetData> matset_info;
  std::vector<WriteTEMPLATE::DirichletSetData> dirset_info;
  std::vector<WriteTEMPLATE::NeumannSetData> neuset_info;

  MeshInfo mesh_info;
  
  matset_info.clear();
  if(gather_mesh_information(mesh_info, matset_info, neuset_info, dirset_info,
                             matsets, neusets, dirsets) != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }


  // try to open the file after gather mesh info succeeds
  if (/* test for file open failure */ true) {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( initialize_file(mesh_info) != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_nodes(mesh_info.num_nodes, mesh_info.nodes, mesh_info.num_dim) != MB_SUCCESS )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_matsets(mesh_info, matset_info, neuset_info) )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

MBErrorCode WriteTEMPLATE::gather_mesh_information(MeshInfo &mesh_info,
                                                   std::vector<WriteTEMPLATE::MaterialSetData> &matset_info,
                                                   std::vector<WriteTEMPLATE::NeumannSetData> &neuset_info,
                                                   std::vector<WriteTEMPLATE::DirichletSetData> &dirset_info,
                                                   std::vector<MBEntityHandle> &matsets,
                                                   std::vector<MBEntityHandle> &neusets,
                                                   std::vector<MBEntityHandle> &dirsets)
{

  std::vector<MBEntityHandle>::iterator vector_iter, end_vector_iter;

  mesh_info.num_nodes = 0;
  mesh_info.num_elements = 0;
  mesh_info.num_matsets = 0;
  
  int id = 0;

  vector_iter= matsets.begin();
  end_vector_iter = matsets.end();

  mesh_info.num_matsets = matsets.size();

  std::vector<MBEntityHandle> parent_meshsets;

  // clean out the bits for the element mark
  mbImpl->tag_delete(mEntityMark);
  mbImpl->tag_create("WriteTEMPLATE element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

  int highest_dimension_of_element_matsets = 0;

  for(vector_iter = matsets.begin(); vector_iter != matsets.end(); vector_iter++)
  {
       
    WriteTEMPLATE::MaterialSetData matset_data;
    matset_data.elements = new MBRange;

    //for the purpose of qa records, get the parents of these matsets 
    if( mbImpl->get_parent_meshsets( *vector_iter, parent_meshsets ) != MB_SUCCESS )
      return MB_FAILURE;

    // get all Entity Handles in the mesh set
    MBRange dummy_range;
    mbImpl->get_entities_by_handle(*vector_iter, dummy_range, true );

      // find the dimension of the last entity in this range
    MBRange::iterator entity_iter = dummy_range.end();
    entity_iter = dummy_range.end();
    entity_iter--;
    int this_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter));
    entity_iter = dummy_range.begin();
    while (entity_iter != dummy_range.end() &&
           MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter)) != this_dim)
      entity_iter++;
    
    if (entity_iter != dummy_range.end())
      std::copy(entity_iter, dummy_range.end(), mb_range_inserter(*(matset_data.elements)));

    assert(matset_data.elements->begin() == matset_data.elements->end() ||
           MBCN::Dimension(TYPE_FROM_HANDLE(*(matset_data.elements->begin()))) == this_dim);
    
    // get the matset's id
    if(mbImpl->tag_get_data(mMaterialSetTag, &(*vector_iter), 1, &id) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get matset id from a tag for an element matset.");
      return MB_FAILURE;
    }
    
    matset_data.id = id; 
    matset_data.number_attributes = 0;
 
     // iterate through all the elements in the meshset
    MBRange::iterator elem_range_iter, end_elem_range_iter;
    elem_range_iter = matset_data.elements->begin();
    end_elem_range_iter = matset_data.elements->end();

      // get the entity type for this matset, verifying that it's the same for all elements
      // THIS ASSUMES HANDLES SORT BY TYPE!!!
    MBEntityType entity_type = TYPE_FROM_HANDLE(*elem_range_iter);
    end_elem_range_iter--;
    if (entity_type != TYPE_FROM_HANDLE(*(end_elem_range_iter++))) {
      mWriteIface->report_error("Entities in matset %i not of common type", id);
      return MB_FAILURE;
    }

    int dimension = MBCN::Dimension(entity_type);

    if( dimension > highest_dimension_of_element_matsets )
      highest_dimension_of_element_matsets = dimension;

    matset_data.moab_type = mbImpl->type_from_handle(*(matset_data.elements->begin()));
    if (MBMAXTYPE == matset_data.moab_type) return MB_FAILURE;
    
    std::vector<MBEntityHandle> tmp_conn;
    mbImpl->get_connectivity(&(*(matset_data.elements->begin())), 1, tmp_conn);
    matset_data.element_type = 
      ExoIIUtil::get_element_type_from_num_verts(tmp_conn.size(), entity_type, dimension);
    
    if (matset_data.element_type == EXOII_MAX_ELEM_TYPE) {
      mWriteIface->report_error("Element type in matset %i didn't get set correctly", id);
      return MB_FAILURE;
    }
    
    matset_data.number_nodes_per_element = ExoIIUtil::VerticesPerElement[matset_data.element_type];

    // number of nodes for this matset
    matset_data.number_elements = matset_data.elements->size();

    // total number of elements
    mesh_info.num_elements += matset_data.number_elements;

    // get the nodes for the elements
    mWriteIface->gather_nodes_from_elements(*matset_data.elements, mEntityMark, mesh_info.nodes);

    if(!neusets.empty())
    {
      // if there are neusets, keep track of which elements are being written out
      for(MBRange::iterator iter = matset_data.elements->begin(); 
          iter != matset_data.elements->end(); ++iter)
      {
        unsigned char bit = 0x1;
        mbImpl->tag_set_data(mEntityMark, &(*iter), 1, &bit);
      }
    }

    matset_info.push_back( matset_data );
  
  }
 

  //if user hasn't entered dimension, we figure it out
  if( mesh_info.num_dim == 0 )
  {
    //never want 1 or zero dimensions
    if( highest_dimension_of_element_matsets < 2 )
      mesh_info.num_dim = 3;
    else
      mesh_info.num_dim = highest_dimension_of_element_matsets;
  }

  MBRange::iterator range_iter, end_range_iter;
  range_iter = mesh_info.nodes.begin();
  end_range_iter = mesh_info.nodes.end();

  mesh_info.num_nodes = mesh_info.nodes.size(); 

  //------dirsets--------
  
  vector_iter= dirsets.begin();
  end_vector_iter = dirsets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {
    
    WriteTEMPLATE::DirichletSetData dirset_data;
    dirset_data.id = 0;
    dirset_data.number_nodes = 0;

    // get the dirset's id
    if(mbImpl->tag_get_data(mDirichletSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS) {
      mWriteIface->report_error("Couldn't get id tag for dirset %i", id);
      return MB_FAILURE;
    }
    
    dirset_data.id = id; 

    std::vector<MBEntityHandle> node_vector;
    //get the nodes of the dirset that are in mesh_info.nodes
    if( mbImpl->get_entities_by_handle(*vector_iter, node_vector, true) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get nodes in dirset %i", id);
      return MB_FAILURE;
    }

    std::vector<MBEntityHandle>::iterator iter, end_iter;
    iter = node_vector.begin();
    end_iter= node_vector.end();
 
    int j=0; 
    unsigned char node_marked = 0;
    MBErrorCode result;
    for(; iter != end_iter; iter++)
    {
      if (TYPE_FROM_HANDLE(*iter) != MBVERTEX) continue;
      result = mbImpl->tag_get_data(mEntityMark, &(*iter), 1, &node_marked);
      if (MB_SUCCESS != result) {
        mWriteIface->report_error("Couldn't get mark data.");
        return result;
      }
      
      if(node_marked == 0x1) dirset_data.nodes.push_back( *iter );    
      j++;
    } 
    
    dirset_data.number_nodes = dirset_data.nodes.size(); 
    dirset_info.push_back( dirset_data );
  }

  //------neusets--------
  vector_iter= neusets.begin();
  end_vector_iter = neusets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {
    WriteTEMPLATE::NeumannSetData neuset_data;

    // get the neuset's id
    if(mbImpl->tag_get_data(mNeumannSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS)
      return MB_FAILURE;

    neuset_data.id = id; 
    neuset_data.mesh_set_handle = *vector_iter; 
 
    //get the sides in two lists, one forward the other reverse; starts with forward sense
      // by convention
    MBRange forward_elems, reverse_elems;
    if(get_neuset_elems(*vector_iter, 0, forward_elems, reverse_elems) == MB_FAILURE)
      return MB_FAILURE;

    MBErrorCode result = get_valid_sides(forward_elems, 1, neuset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    result = get_valid_sides(reverse_elems, -1, neuset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    
    neuset_data.number_elements = neuset_data.elements.size(); 
    neuset_info.push_back( neuset_data );
  }

  return MB_SUCCESS;
}

MBErrorCode WriteTEMPLATE::get_valid_sides(MBRange &elems, const int sense,
                                           WriteTEMPLATE::NeumannSetData &neuset_data) 
{
    // this is where we see if underlying element of side set element is included in output 

  unsigned char element_marked = 0;
  MBErrorCode result;
  for(MBRange::iterator iter = elems.begin(); iter != elems.end(); iter++)
  {
      // should insert here if "side" is a quad/tri on a quad/tri mesh
    result = mbImpl->tag_get_data(mEntityMark, &(*iter), 1, &element_marked);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get mark data.");
      return result;
    }
    
    if(element_marked == 0x1)
    {
      neuset_data.elements.push_back( *iter );

        // TJT TODO: the sense should really be # edges + 1or2
      neuset_data.side_numbers.push_back((sense == 1 ? 1 : 2));
    }
    else //then "side" is probably a quad/tri on a hex/tet mesh
    {
      std::vector<MBEntityHandle> parents;
      int dimension = MBCN::Dimension( TYPE_FROM_HANDLE(*iter));

        //get the adjacent parent element of "side"
      if( mbImpl->get_adjacencies( &(*iter), 1, dimension+1, false, parents) != MB_SUCCESS ) {
        mWriteIface->report_error("Couldn't get adjacencies for neuset.");
        return MB_FAILURE;
      }
       
      if(!parents.empty())     
      {
          //make sure the adjacent parent element will be output
        for(unsigned int k=0; k<parents.size(); k++)
        {
          result = mbImpl->tag_get_data(mEntityMark, &(parents[k]), 1, &element_marked);
          if (MB_SUCCESS != result) {
            mWriteIface->report_error("Couldn't get mark data.");
            return result;
          }
        
          int side_no, this_sense, this_offset;
          if(element_marked == 0x1 &&
             mbImpl->side_number(parents[k], *iter, side_no, 
                                  this_sense, this_offset) == MB_SUCCESS &&
             this_sense == sense) {
            neuset_data.elements.push_back(parents[k]);
            neuset_data.side_numbers.push_back(side_no+1);
            break;
          }
        }
      }
      else
      {
        mWriteIface->report_error("No parent element exists for element in neuset %i", neuset_data.id);
        return MB_FAILURE;
      }
    }
  }

  return MB_SUCCESS;
}

MBErrorCode WriteTEMPLATE::write_nodes(const int num_nodes, const MBRange& nodes, const int dimension)
{
  //see if should transform coordinates
  MBErrorCode result;
  MBTag trans_tag;
  result = mbImpl->tag_get_handle( MESH_TRANSFORM_TAG_NAME, trans_tag);
  bool transform_needed = true;
  if( result == MB_TAG_NOT_FOUND )
    transform_needed = false;

  int num_coords_to_fill = transform_needed ? 3 : dimension;

  std::vector<double*> coord_arrays(3);
  coord_arrays[0] = new double[num_nodes];
  coord_arrays[1] = new double[num_nodes];
  coord_arrays[2] = NULL;

  if( num_coords_to_fill == 3 ) 
    coord_arrays[2] = new double[num_nodes];
 
  result = mWriteIface->get_node_arrays(dimension, num_nodes, nodes, 
                                        mGlobalIdTag, 0, coord_arrays);
  if(result != MB_SUCCESS)
  {
    delete [] coord_arrays[0];
    delete [] coord_arrays[1];
    if(coord_arrays[2]) delete [] coord_arrays[2];
    return result;
  }

  if( transform_needed )
  {
    double trans_matrix[16]; 
    result = mbImpl->tag_get_data( trans_tag, NULL, 0, trans_matrix ); 
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get transform data.");
      return result;
    }
      
    for( int i=0; i<num_nodes; i++)
    {

      double vec1[3];
      double vec2[3];

      vec2[0] =  coord_arrays[0][i];
      vec2[1] =  coord_arrays[1][i];
      vec2[2] =  coord_arrays[2][i];

      for( int row=0; row<3; row++ )
      {
        vec1[row] = 0.0;
        for( int col = 0; col<3; col++ )
        {
          vec1[row] += ( trans_matrix[ (row*4)+col ] * vec2[col] );
        }
      }

      coord_arrays[0][i] = vec1[0];
      coord_arrays[1][i] = vec1[1];
      coord_arrays[2][i] = vec1[2];

    }
  }


  // write the nodes 

  /* template - write nodes to file here in some way */

  // clean up
  delete [] coord_arrays[0];
  delete [] coord_arrays[1];
  if(coord_arrays[2]) 
    delete [] coord_arrays[2];

  return MB_SUCCESS;

}

MBErrorCode WriteTEMPLATE::write_matsets(MeshInfo & /* mesh_info (commented out to remove warning) */,
                                         std::vector<WriteTEMPLATE::MaterialSetData> &matset_data,
                                         std::vector<WriteTEMPLATE::NeumannSetData> &/* neuset_data (commented out to remove warning) */)
{

  unsigned int i;
  std::vector<int> connect;
  const MBEntityHandle *connecth;
  int num_connecth;
  MBErrorCode result;
  
    // don't usually have anywhere near 31 nodes per element
  connect.reserve(31);
  MBRange::iterator rit;

  WriteTEMPLATE::MaterialSetData matset;
  for (i = 0; i < matset_data.size(); i++) {
    matset = matset_data[i];
    
    int id = matset.id;
      // bogus line to get rid of warning
    if (0 == id) ;

    for (rit = matset.elements->begin(); rit != matset.elements->end(); rit++) {
      
        // get the connectivity of this element
      result = mbImpl->get_connectivity(*rit, connecth, num_connecth);
      if (MB_SUCCESS != result) return result;
      
        // get the vertex ids
      result = mbImpl->tag_get_data(mGlobalIdTag, connecth, num_connecth, &connect[0]);
      if (MB_SUCCESS != result) return result;
      
        // write the data
        /* template - write element connectivity here */

      if(/* template - check for error condition! */ false)
        return MB_FAILURE;
    }
  }

  return MB_SUCCESS;
}

MBErrorCode WriteTEMPLATE::initialize_file(MeshInfo &mesh_info)
{
    // perform the initializations

  int coord_size, ncoords;
  
  coord_size = mesh_info.num_dim;
    /* template - write coord size */

  ncoords = mesh_info.num_nodes;
    /* template - write num nodes*/
  
    /* template - write information on the element types & numbers (depends
       on material and other sets) */

/* node coordinate arrays: */
    /* template - initialize variable to hold coordinate arrays */

   return MB_SUCCESS;
}


MBErrorCode WriteTEMPLATE::open_file(const char* filename)
{
   // not a valid filname
   if(strlen((const char*)filename) == 0)
   {
     mWriteIface->report_error("Output filename not specified");
      return MB_FAILURE;
   }

     /* template - open file & store somewhere */

   // file couldn't be opened
   if(/* template - check for file open error here! */ false)
   {
     mWriteIface->report_error("Cannot open %s", filename);
     return MB_FAILURE;
   }
   return MB_SUCCESS;
}

MBErrorCode WriteTEMPLATE::get_neuset_elems(MBEntityHandle neuset, int current_sense,
                                        MBRange &forward_elems, MBRange &reverse_elems) 
{
  MBRange neuset_elems, neuset_meshsets;

    // get the sense tag; don't need to check return, might be an error if the tag
    // hasn't been created yet
  MBTag sense_tag = 0;
  mbImpl->tag_get_handle("SENSE", sense_tag);

    // get the entities in this set
  MBErrorCode result = mbImpl->get_entities_by_handle(neuset, neuset_elems, true);
  if (MB_FAILURE == result) return result;
  
    // now remove the meshsets into the neuset_meshsets; first find the first meshset,
  MBRange::iterator range_iter = neuset_elems.begin();
  while (TYPE_FROM_HANDLE(*range_iter) != MBENTITYSET && range_iter != neuset_elems.end())
    range_iter++;
  
    // then, if there are some, copy them into neuset_meshsets and erase from neuset_elems
  if (range_iter != neuset_elems.end()) {
    std::copy(range_iter, neuset_elems.end(), mb_range_inserter(neuset_meshsets));
    neuset_elems.erase(range_iter, neuset_elems.end());
  }
  

    // ok, for the elements, check the sense of this set and copy into the right range
    // (if the sense is 0, copy into both ranges)

    // need to step forward on list until we reach the right dimension
  MBRange::iterator dum_it = neuset_elems.end();
  dum_it--;
  int target_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*dum_it));
  dum_it = neuset_elems.begin();
  while (target_dim != MBCN::Dimension(TYPE_FROM_HANDLE(*dum_it)) &&
         dum_it != neuset_elems.end()) 
    dum_it++;

  if (current_sense == 1 || current_sense == 0)
    std::copy(dum_it, neuset_elems.end(), mb_range_inserter(forward_elems));
  if (current_sense == -1 || current_sense == 0)
    std::copy(dum_it, neuset_elems.end(), mb_range_inserter(reverse_elems));
  
    // now loop over the contained meshsets, getting the sense of those and calling this
    // function recursively
  for (range_iter = neuset_meshsets.begin(); range_iter != neuset_meshsets.end(); range_iter++) {

      // first get the sense; if it's not there, by convention it's forward
    int this_sense;
    if (0 == sense_tag ||
        MB_FAILURE == mbImpl->tag_get_data(sense_tag, &(*range_iter), 1, &this_sense))
      this_sense = 1;
      
      // now get all the entities on this meshset, with the proper (possibly reversed) sense
    get_neuset_elems(*range_iter, this_sense*current_sense,
                      forward_elems, reverse_elems);
  }
  
  return result;
}


  