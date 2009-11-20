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

#ifdef __MFC_VER
#pragma warning(disable:4786)
#endif

#include "MBSkinner.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <iostream>
#include "MBUtil.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "AEntityFactory.hpp"

#ifdef M_PI
#  define SKINNER_PI M_PI
#else
#  define SKINNER_PI 3.1415926535897932384626
#endif



MBSkinner::~MBSkinner()
{
  // delete the adjacency tag
}


void MBSkinner::initialize()
{
  // go through and mark all the target dimension entities
  // that already exist as not deleteable
  // also get the connectivity tags for each type
  // also populate adjacency information
  MBEntityType type;
  MBDimensionPair target_ent_types = MBCN::TypeDimensionMap[mTargetDim];

  void* null_ptr = NULL;
  MBErrorCode result;

  result = thisMB->tag_create("skinner adj", sizeof(void*), MB_TAG_DENSE, mAdjTag, &null_ptr);
  assert(MB_SUCCESS == result);

  if(mDeletableMBTag == 0) {
    result = thisMB->tag_create("skinner deletable", 1, MB_TAG_BIT, mDeletableMBTag, NULL);
    assert(MB_SUCCESS == result);
  }
  
  MBRange entities;

    // go through each type at this dimension 
  for(type = target_ent_types.first; type <= target_ent_types.second; ++type)
  {
      // get the entities of this type in the MB
    thisMB->get_entities_by_type(0, type, entities);

      // go through each entity of this type in the MB
      // and set its deletable tag to NO
    MBRange::iterator iter, end_iter;
    end_iter = entities.end();
    for(iter = entities.begin(); iter != end_iter; ++iter)
    {
      unsigned char bit = 0x1;
      result = thisMB->tag_set_data(mDeletableMBTag, &(*iter), 1, &bit);
      assert(MB_SUCCESS == result);
        // add adjacency information too
      if (TYPE_FROM_HANDLE(*iter) != MBVERTEX)
        add_adjacency(*iter);
    }
  }
}

void MBSkinner::deinitialize()
{
  MBErrorCode result = MB_SUCCESS;
  
  if (0 != mDeletableMBTag) {
    result = thisMB->tag_delete( mDeletableMBTag);
    mDeletableMBTag = 0;
    assert(MB_SUCCESS == result);
  }

  // remove the adjaceny tag
  std::vector< std::vector<MBEntityHandle>* > adj_arr;
  std::vector< std::vector<MBEntityHandle>* >::iterator i;
  if (0 != mAdjTag) {
    for (MBEntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
      MBRange entities;
      result = thisMB->get_entities_by_type_and_tag( 0, t, &mAdjTag, 0, 1, entities );
      assert(MB_SUCCESS == result);
      adj_arr.resize( entities.size() );
      result = thisMB->tag_get_data( mAdjTag, entities, &adj_arr[0] );
      assert(MB_SUCCESS == result);
      for (i = adj_arr.begin(); i != adj_arr.end(); ++i)
        delete *i;
    }
  
    result = thisMB->tag_delete(mAdjTag);
    mAdjTag = 0;
    assert(MB_SUCCESS == result);
  }
}


void MBSkinner::add_adjacency(MBEntityHandle entity)
{
  std::vector<MBEntityHandle> *adj = NULL;
  const MBEntityHandle *nodes;
  int num_nodes;
  MBErrorCode result = thisMB->get_connectivity(entity, nodes, num_nodes);
  assert(MB_SUCCESS == result);
  const MBEntityHandle *iter =
    std::min_element(nodes, nodes+num_nodes);

  if(iter == nodes+num_nodes)
    return;

  // add this entity to the node
  if(thisMB->tag_get_data(mAdjTag, iter, 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    adj->push_back(entity);    
  }
  // create a new vector and add it
  else
  {
    adj = new std::vector<MBEntityHandle>;
    adj->push_back(entity);
    result = thisMB->tag_set_data(mAdjTag, iter, 1, &adj);
    assert(MB_SUCCESS == result);
  }
}

void MBSkinner::add_adjacency(MBEntityHandle entity, 
                               const MBEntityHandle *nodes,
                               const int num_nodes)
{
  std::vector<MBEntityHandle> *adj = NULL;
  const MBEntityHandle *iter = 
    std::min_element(nodes, nodes+num_nodes);

  if(iter == nodes+num_nodes)
    return;

  // add this entity to the node
  if(thisMB->tag_get_data(mAdjTag, iter, 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    adj->push_back(entity);    
  }
  // create a new vector and add it
  else
  {
    adj = new std::vector<MBEntityHandle>;
    adj->push_back(entity);
    thisMB->tag_set_data(mAdjTag, iter, 1, &adj);
  }
}

MBErrorCode MBSkinner::find_geometric_skin(MBRange &forward_target_entities) 
{
    // attempts to find whole model skin, using geom topo sets first then
    // normal find_skin function
  bool debug = true;

    // look for geom topo sets
  MBTag geom_tag;
  MBErrorCode result = thisMB->tag_create(GEOM_DIMENSION_TAG_NAME, 4, 
                                            MB_TAG_SPARSE, geom_tag, NULL);

  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
    return result;
  
    // get face sets (dimension = 2)
  MBRange face_sets;
  int two = 2;
  const void *two_ptr = &two;
  result = thisMB->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag, &two_ptr, 1,
                                                 face_sets);

  MBRange::iterator it;
  if (MB_SUCCESS != result)
    return result;
  else if (face_sets.empty())
    return MB_ENTITY_NOT_FOUND;

    // ok, we have face sets; use those to determine skin
  MBRange skin_sets;
  if (debug) std::cout << "Found " << face_sets.size() << " face sets total..." << std::endl;
  
  for (it = face_sets.begin(); it != face_sets.end(); it++) {
    int num_parents;
    result = thisMB->num_parent_meshsets(*it, &num_parents);
    if (MB_SUCCESS != result)
      return result;
    else if (num_parents == 1)
      skin_sets.insert(*it);
  }

  if (debug) std::cout << "Found " << skin_sets.size() << " 1-parent face sets..." << std::endl;

  if (skin_sets.empty())
    return MB_FAILURE;
      
    // ok, we have the shell; gather up the elements, putting them all in forward for now
  for (it = skin_sets.begin(); it != skin_sets.end(); it++) {
    result = thisMB->get_entities_by_handle(*it, forward_target_entities, true);
    if (MB_SUCCESS != result) 
      return result;
  }
        
  return result;
}

MBErrorCode MBSkinner::find_skin( const MBRange& source_entities,
                                  bool get_vertices,
                                  MBRange& output_handles,
                                  MBRange* output_reverse_handles,
                                  bool create_vert_elem_adjs,
                                  bool create_skin_elements )
{
  if (source_entities.empty())
    return MB_SUCCESS;

  MBCore* this_core = dynamic_cast<MBCore*>(thisMB);
  if (this_core && create_vert_elem_adjs && 
      !this_core->a_entity_factory()->vert_elem_adjacencies())
    this_core->a_entity_factory()->create_vert_elem_adjacencies();
    
  if (this_core && this_core->a_entity_factory()->vert_elem_adjacencies())
    return find_skin_vertices( source_entities, 
                               get_vertices ? &output_handles : 0,
                               get_vertices ? 0 : &output_handles,
                               output_reverse_handles,
                               create_skin_elements );
  
  MBRange forward, reverse;
  MBRange prev;
  const int d = MBCN::Dimension(TYPE_FROM_HANDLE(source_entities.front()));
  if (!source_entities.all_of_dimension(d))
    return MB_TYPE_OUT_OF_RANGE;
  
  MBErrorCode rval = thisMB->get_entities_by_dimension( 0, d-1, prev );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = find_skin_noadj( source_entities, forward, reverse );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (get_vertices && !output_reverse_handles) {
    forward.merge( reverse );
    reverse.clear();
  }
  
  if (get_vertices) {
    rval = thisMB->get_connectivity( forward, output_handles );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  if (!create_skin_elements) {
    MBRange new_skin;
    rval = thisMB->get_entities_by_dimension( 0, d-1, new_skin);
    if (MB_SUCCESS != rval)
      return rval;
    new_skin = subtract( new_skin, prev );
    forward = subtract( forward, new_skin );
    reverse = subtract( reverse, new_skin );
    rval = thisMB->delete_entities( new_skin );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  if (!get_vertices) {
    if (output_handles.empty())
      output_handles.swap( forward );
    else
      output_handles.merge( forward );
    if (output_reverse_handles) {
      if (output_reverse_handles->empty())
        output_reverse_handles->swap( reverse );
      else
        output_reverse_handles->merge( reverse );
    }
  }
  
  return MB_SUCCESS;  
}

MBErrorCode MBSkinner::find_skin_noadj(const MBRange &source_entities,
                                 MBRange &forward_target_entities,
                                 MBRange &reverse_target_entities/*,
                                 bool create_vert_elem_adjs*/)
{
  if(source_entities.empty())
    return MB_FAILURE;
  
  // get our working dimensions
  MBEntityType type = thisMB->type_from_handle(*(source_entities.begin()));
  const int source_dim = MBCN::Dimension(type);
  mTargetDim = source_dim - 1;

  // make sure we can handle the working dimensions
  if(mTargetDim < 0 || source_dim > 3)
    return MB_FAILURE;

  //MBCore *this_core = dynamic_cast<MBCore*>(thisMB);
  //bool use_adjs = false;
  //if (!this_core->a_entity_factory()->vert_elem_adjacencies() &&
  //  create_vert_elem_adjs)
  //  this_core->a_entity_factory()->create_vert_elem_adjacencies();
  //
  //if (this_core->a_entity_factory()->vert_elem_adjacencies())
  //  use_adjs = true;
  //
  //else 
  initialize();

  MBRange::const_iterator iter, end_iter;
  end_iter = source_entities.end();
  const MBEntityHandle *conn;
  MBEntityHandle match;

  direction direct;
  MBErrorCode result;
    // assume we'll never have more than 32 vertices on a facet (checked
    // with assert later)
  MBEntityHandle sub_conn[32];
  std::vector<MBEntityHandle> tmp_conn_vec;
  int num_nodes, num_sub_nodes, num_sides;
  const short *sub_indices;
  MBEntityType sub_type;

  // for each source entity
  for(iter = source_entities.begin(); iter != end_iter; ++iter)
  {
    // get the connectivity of this entity
    result = thisMB->get_connectivity(*iter, conn, num_nodes, false, &tmp_conn_vec);
    if (MB_SUCCESS != result)
      return result;
    
    type = thisMB->type_from_handle(*iter);
    MBRange::iterator seek_iter;
    MBRange dum_elems, dum_sub_elems;
    
    // get connectivity of each n-1 dimension entity
    num_sides = MBCN::NumSubEntities( type, mTargetDim );
    for(int i=0; i<num_sides; i++)
    {
      sub_indices = MBCN::SubEntityVertexIndices( type, mTargetDim, i, sub_type, num_sub_nodes );
      assert(num_sub_nodes <= 32);
      for(int j=0; j<num_sub_nodes; j++)
        sub_conn[j] = conn[sub_indices[j]];
      
//      if (use_adjs) {
//        dum_elems.clear();
//        result = thisMB->get_adjacencies(sub_conn, num_sub_nodes, source_dim, false,
//                                         dum_elems);
//        if (MB_SUCCESS != result) return result;
//        dum_elems = intersect( dum_elems, source_entities );
//        if (dum_elems.empty()) {
//          assert(false);  // should never happen
//          return MB_FAILURE;
//        }
//        // if (dum_elems.size() > 2) { 
//          // source entities do not form a valid source-dim patch (t-junction).
//          // do we care?
//        // }
//        
//        if (1 == dum_elems.size()) {
//            // this sub_element is on the skin
//
//            // check for existing entity
//          dum_sub_elems.clear();
//          result = thisMB->get_adjacencies(sub_conn, num_sub_nodes, mTargetDim, false,
//                                           dum_sub_elems);
//          if (MB_SUCCESS != result) return result;
//          if (dum_sub_elems.empty()) {
//              // need to create one
//            MBEntityHandle tmphndl=0;
//            int indices[MB_MAX_SUB_ENTITY_VERTICES];
//            MBEntityType new_type;
//            int num_new_nodes;
//            MBCN::SubEntityNodeIndices( type, num_nodes, mTargetDim, i, new_type, num_new_nodes, indices );
//            for(int j=0; j<num_new_nodes; j++)
//              sub_conn[j] = conn[indices[j]];
//        
//            result = thisMB->create_element(new_type, sub_conn,  
//                                            num_new_nodes, tmphndl);
//            forward_target_entities.insert(tmphndl);
//          }
//          else {
//              // else find the relative sense of this entity to the source_entity in this set
//            int side_no, sense = 0, offset;
//            if (source_entities.find(*dum_elems.begin()) == source_entities.end()) {
//              result = thisMB->side_number(*dum_elems.rbegin(), *dum_sub_elems.begin(),
//                                           side_no, sense, offset);
//            }
//            else {
//              result = thisMB->side_number(*dum_elems.begin(), *dum_sub_elems.begin(),
//                                           side_no, sense, offset);
//            }
//            if (-1 == sense) reverse_target_entities.insert(*dum_sub_elems.begin());
//            else if (1 == sense) forward_target_entities.insert(*dum_sub_elems.begin());
//            else return MB_FAILURE;
//          }
//        }
//      }
//      else {
        
          // see if we can match this connectivity with
          // an existing entity
        find_match( sub_type, sub_conn, num_sub_nodes, match, direct );
        
          // if there is no match, create a new entity
        if(match == 0)
        {
          MBEntityHandle tmphndl=0;
          int indices[MB_MAX_SUB_ENTITY_VERTICES];
          MBEntityType new_type;
          int num_new_nodes;
          MBCN::SubEntityNodeIndices( type, num_nodes, mTargetDim, i, new_type, num_new_nodes, indices );
          for(int j=0; j<num_new_nodes; j++)
            sub_conn[j] = conn[indices[j]];
          result = thisMB->create_element(new_type, sub_conn, num_new_nodes,
                                          tmphndl);
          assert(MB_SUCCESS == result);
          add_adjacency(tmphndl, sub_conn, num_sub_nodes);
          forward_target_entities.insert(tmphndl);
        }
          // if there is a match, delete the matching entity
          // if we can. 
        else
        {
          if ( (seek_iter = forward_target_entities.find(match)) != forward_target_entities.end())
          {
            forward_target_entities.erase(seek_iter);
            remove_adjacency(match);
            if(/*!use_adjs &&*/ entity_deletable(match))
            {
              result = thisMB->delete_entities(&match, 1);
              assert(MB_SUCCESS == result);
            }
          }
          else if ( (seek_iter = reverse_target_entities.find(match)) != reverse_target_entities.end())
          {
            reverse_target_entities.erase(seek_iter);
            remove_adjacency(match);
            if(/*!use_adjs &&*/ entity_deletable(match))
            {
              result = thisMB->delete_entities(&match, 1);
              assert(MB_SUCCESS == result);
            }
          }
          else
          {
            if(direct == FORWARD)
            {
              forward_target_entities.insert(match);
            }
            else
            {
              reverse_target_entities.insert(match);
            }
          }
        }
      //}
    }
  }

  deinitialize();

  return MB_SUCCESS;
}


void MBSkinner::find_match( MBEntityType type, 
                             const MBEntityHandle *conn,
                             const int num_nodes,
                             MBEntityHandle& match,
                             MBSkinner::direction &direct)
{
  match = 0;

  if (type == MBVERTEX) {
    match = *conn;
    direct = FORWARD;
    return;
  }

  const MBEntityHandle *iter = std::min_element(conn, conn+num_nodes);

  std::vector<MBEntityHandle> *adj = NULL;

  MBErrorCode result = thisMB->tag_get_data(mAdjTag, iter, 1, &adj);
  if(result == MB_FAILURE || adj == NULL)
  {
    return;
  }

  std::vector<MBEntityHandle>::iterator jter, end_jter;
  end_jter = adj->end();

  const MBEntityHandle *tmp;
  int num_verts;

  for(jter = adj->begin(); jter != end_jter; ++jter)
  {
    MBEntityType tmp_type;
    tmp_type = thisMB->type_from_handle(*jter);

    if( type != tmp_type )
      continue;

    result = thisMB->get_connectivity(*jter, tmp, num_verts, true);
    assert(MB_SUCCESS == result && num_verts >= num_nodes);
    if(connectivity_match(conn, tmp, num_verts, direct))
    {
      match = *jter;
      break;
    }        
  }
}

bool MBSkinner::connectivity_match( const MBEntityHandle *conn1,
                                     const MBEntityHandle *conn2,
                                     const int num_verts,
                                     MBSkinner::direction &direct)
{
  const MBEntityHandle *iter =
    std::find(conn2, conn2+num_verts, conn1[0]);
  if(iter == conn2+num_verts)
    return false;

  bool they_match = true;

  int i;
  unsigned int j = iter - conn2;
    
  // first compare forward
  for(i = 1; i<num_verts; ++i)
  {
    if(conn1[i] != conn2[(j+i)%num_verts])
    {
      they_match = false;
      break;
    }
  }
  
  if(they_match == true)
  {
    // need to check for reversed edges here
    direct = (num_verts == 2 && j) ? REVERSE : FORWARD;
    return true;
  }
  
  they_match = true;
  
  // then compare reverse
  j += num_verts;
  for(i = 1; i < num_verts; )
  {
    if(conn1[i] != conn2[(j-i)%num_verts])
    {
      they_match = false;
      break;
    }
    ++i;
  }
  if (they_match)
  {
    direct = REVERSE;
  }
  return they_match;
}

  
MBErrorCode MBSkinner::remove_adjacency(MBEntityHandle entity)
{
  std::vector<MBEntityHandle> nodes, *adj = NULL;
  MBErrorCode result = thisMB->get_connectivity(&entity, 1, nodes);
  if (MB_SUCCESS != result) return result;
  std::vector<MBEntityHandle>::iterator iter = 
    std::min_element(nodes.begin(), nodes.end());

  if(iter == nodes.end())
    return MB_FAILURE;

  // remove this entity from the node
  if(thisMB->tag_get_data(mAdjTag, &(*iter), 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    iter = std::find(adj->begin(), adj->end(), entity);
    if(iter != adj->end())
      adj->erase(iter);
  }

  return result;
}

bool MBSkinner::entity_deletable(MBEntityHandle entity)
{
  unsigned char deletable=0;
  MBErrorCode result = thisMB->tag_get_data(mDeletableMBTag, &entity, 1, &deletable);
  assert(MB_SUCCESS == result);
  if(MB_SUCCESS == result && deletable == 1)
    return false;
  return true;
}

MBErrorCode MBSkinner::classify_2d_boundary( const MBRange &boundary,
                                               const MBRange &bar_elements,
                                               MBEntityHandle boundary_edges,
                                               MBEntityHandle inferred_edges,
                                               MBEntityHandle non_manifold_edges,
                                               MBEntityHandle other_edges,
                                               int &number_boundary_nodes)
{
  MBRange bedges, iedges, nmedges, oedges;
  MBErrorCode result = classify_2d_boundary(boundary, bar_elements,
                                             bedges, iedges, nmedges, oedges,
                                             number_boundary_nodes);
  if (MB_SUCCESS != result) return result;
  
    // now set the input meshsets to the output ranges
  result = thisMB->clear_meshset(&boundary_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(boundary_edges, bedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&inferred_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(inferred_edges, iedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&non_manifold_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(non_manifold_edges, nmedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&other_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(other_edges, oedges);
  if (MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}

MBErrorCode MBSkinner::classify_2d_boundary( const MBRange &boundary,
                                               const MBRange &bar_elements,
                                               MBRange &boundary_edges,
                                               MBRange &inferred_edges,
                                               MBRange &non_manifold_edges,
                                               MBRange &other_edges,
                                               int &number_boundary_nodes)
{

  // clear out the edge lists

  boundary_edges.clear();
  inferred_edges.clear();
  non_manifold_edges.clear();
  other_edges.clear();

  number_boundary_nodes = 0;

  // make sure we have something to work with
  if(boundary.empty())
  {
    return MB_FAILURE;
  }
  
  // get our working dimensions
  MBEntityType type = thisMB->type_from_handle(*(boundary.begin()));
  const int source_dim = MBCN::Dimension(type);

  // make sure we can handle the working dimensions
  if(source_dim != 2)
  {
    return MB_FAILURE;
  }
  mTargetDim = source_dim - 1;

  // initialize
  initialize();

  // additional initialization for this routine
  // define a tag for MBEDGE which counts the occurances of the edge below
  // default should be 0 for existing edges, if any

  MBTag count_tag;
  int default_count = 0;
  MBErrorCode result = thisMB->tag_create("mdbskinner count edges", sizeof(int),
                                            MB_TAG_DENSE, count_tag, &default_count);
  assert(MB_SUCCESS == result);

 
  MBRange::const_iterator iter, end_iter;
  end_iter = boundary.end();

  std::vector<MBEntityHandle> conn;
  MBEntityHandle sub_conn[2];
  MBEntityHandle match;

  MBRange edge_list;
  MBRange boundary_nodes;
  MBSkinner::direction direct;
  
  MBEntityType sub_type;
  int num_edge, num_sub_ent_vert;
  const short* edge_verts;
  
  // now, process each entity in the boundary

  for(iter = boundary.begin(); iter != end_iter; ++iter)
  {
    // get the connectivity of this entity
    conn.clear();
    result = thisMB->get_connectivity(&(*iter), 1, conn, false);
    assert(MB_SUCCESS == result);

    // add node handles to boundary_node range
    std::copy(conn.begin(), conn.begin()+MBCN::VerticesPerEntity(type), 
              mb_range_inserter(boundary_nodes));

    type = thisMB->type_from_handle(*iter);
    
    // get connectivity of each n-1 dimension entity (edge in this case)
    const struct MBCN::ConnMap* conn_map = &(MBCN::mConnectivityMap[type][0]);
    num_edge = MBCN::NumSubEntities( type, 1 );
    for(int i=0; i<num_edge; i++)
    {
      edge_verts = MBCN::SubEntityVertexIndices( type, 1, i, sub_type, num_sub_ent_vert );
      assert( sub_type == MBEDGE && num_sub_ent_vert == 2 );
      sub_conn[0] = conn[edge_verts[0]];
      sub_conn[1] = conn[edge_verts[1]];
      int num_sub_nodes = conn_map->num_corners_per_sub_element[i];
      
      // see if we can match this connectivity with
      // an existing entity
      find_match( MBEDGE, sub_conn, num_sub_nodes, match, direct );
  
      // if there is no match, create a new entity
      if(match == 0)
      {
        MBEntityHandle tmphndl=0;
        int indices[MB_MAX_SUB_ENTITY_VERTICES];
        MBEntityType new_type;
        int num_new_nodes;
        MBCN::SubEntityNodeIndices( type, conn.size(), 1, i, new_type, num_new_nodes, indices );
        for(int j=0; j<num_new_nodes; j++)
          sub_conn[j] = conn[indices[j]];
        
        result = thisMB->create_element(new_type, sub_conn,  
                                        num_new_nodes, tmphndl);
        assert(MB_SUCCESS == result);
        add_adjacency(tmphndl, sub_conn, num_sub_nodes);
        //target_entities.insert(tmphndl);
        edge_list.insert(tmphndl);
        int count;
        result = thisMB->tag_get_data(count_tag, &tmphndl, 1, &count);
        assert(MB_SUCCESS == result);
        count++;
        result = thisMB->tag_set_data(count_tag, &tmphndl, 1, &count);
        assert(MB_SUCCESS == result);

      }
      else
      {
        // We found a match, we must increment the count on the match
        int count;
        result = thisMB->tag_get_data(count_tag, &match, 1, &count);
        assert(MB_SUCCESS == result);
        count++;
        result = thisMB->tag_set_data(count_tag, &match, 1, &count);
        assert(MB_SUCCESS == result);

        // if the entity is not deletable, it was pre-existing in
        // the database.  We therefore may need to add it to the
        // edge_list.  Since it will not hurt the range, we add
        // whether it was added before or not
        if(!entity_deletable(match))
        {
          edge_list.insert(match);
        }
      }
    }
  }

  // Any bar elements in the model should be classified separately
  // If the element is in the skin edge_list, then it should be put in
  // the non-manifold edge list.  Edges not in the edge_list are stand-alone
  // bars, and we make them simply boundary elements

  if (!bar_elements.empty())
  {
    MBRange::iterator bar_iter;
    for(iter = bar_elements.begin(); iter != bar_elements.end(); ++iter)
    {
      MBEntityHandle handle = *iter;
      bar_iter = edge_list.find(handle);
      if (bar_iter != edge_list.end())
      {
        // it is in the list, erase it and put in non-manifold list
        edge_list.erase(bar_iter);
        non_manifold_edges.insert(handle);
      }
      else
      {
        // not in the edge list, make it a boundary edge
        boundary_edges.insert(handle);
      }
    }
  }

  // now all edges should be classified.  Go through the edge_list,
  // and put all in the appropriate lists

  MBRange::iterator edge_iter, edge_end_iter;
  edge_end_iter = edge_list.end();
  int count;
  for(edge_iter = edge_list.begin(); edge_iter != edge_end_iter; edge_iter++)
  {
    // check the count_tag
    result = thisMB->tag_get_data(count_tag, &(*edge_iter), 1, &count);
    assert(MB_SUCCESS == result);
    if (count == 1)
    {
      boundary_edges.insert(*edge_iter);
   }
    else if (count == 2)
    {
      other_edges.insert(*edge_iter);
    }
    else
    {
      non_manifold_edges.insert(*edge_iter);
    }
  }

  // find the inferred edges from the other_edge_list

  double min_angle_degrees = 20.0;
  find_inferred_edges(const_cast<MBRange&> (boundary), other_edges, inferred_edges, min_angle_degrees);

  // we now want to remove the inferred_edges from the other_edges

  MBRange temp_range;
 
  std::set_difference(other_edges.begin(), other_edges.end(),
                      inferred_edges.begin(), inferred_edges.end(),
                      mb_range_inserter(temp_range),
                      std::less<MBEntityHandle>() );

  other_edges = temp_range;

  // get rid of count tag and deinitialize

  result = thisMB->tag_delete(count_tag);
  assert(MB_SUCCESS == result);
  deinitialize();

  // set the node count
  number_boundary_nodes = boundary_nodes.size();

  return MB_SUCCESS;
} 

void MBSkinner::find_inferred_edges(MBRange &skin_boundary,
                                     MBRange &candidate_edges,
                                     MBRange &inferred_edges,
                                     double reference_angle_degrees)
{

  // mark all the entities in the skin boundary
  MBTag mark_tag;
  MBErrorCode result = thisMB->tag_create("find inferred edges mark", 1, MB_TAG_BIT, mark_tag, NULL);
  assert(MB_SUCCESS == result);
  for(MBRange::iterator mark_iter = skin_boundary.begin();
      mark_iter != skin_boundary.end(); ++mark_iter)
  {
    unsigned char bit = true;
    result = thisMB->tag_set_data(mark_tag, &(*mark_iter), 1, &bit);
    assert(MB_SUCCESS == result);
  }

  // find the cosine of the reference angle

  double reference_cosine = cos(reference_angle_degrees*SKINNER_PI/180.0);
  
  // check all candidate edges for an angle greater than the minimum

  MBRange::iterator iter, end_iter = candidate_edges.end();
  std::vector<MBEntityHandle> adjacencies;
  std::vector<MBEntityHandle>::iterator adj_iter;
  MBEntityHandle face[2];

  for(iter = candidate_edges.begin(); iter != end_iter; ++iter)
  {

    // get the 2D elements connected to this edge
    adjacencies.clear();
    result = thisMB->get_adjacencies(&(*iter), 1, 2, false, adjacencies);
    if (MB_SUCCESS != result) 
      continue;

    // there should be exactly two, that is why the edge is classified as nonBoundary
    // and manifold

    int faces_found = 0;
    for(adj_iter = adjacencies.begin(); adj_iter != adjacencies.end() && faces_found < 2; ++adj_iter)
    {
      // we need to find two of these which are in the skin
      unsigned char is_marked = 0;
      result = thisMB->tag_get_data(mark_tag, &(*adj_iter), 1, &is_marked);
      assert(MB_SUCCESS == result);
      if(is_marked)
      {
        face[faces_found] = *adj_iter;
        faces_found++;
      } 
    }

//    assert(faces_found == 2 || faces_found == 0);
    if (2 != faces_found) 
      continue;

    // see if the two entities have a sufficient angle

    if ( has_larger_angle(face[0], face[1], reference_cosine) )
    {
       inferred_edges.insert(*iter);
    }
  }
  
  result = thisMB->tag_delete(mark_tag);
  assert(MB_SUCCESS == result);
}

bool MBSkinner::has_larger_angle(MBEntityHandle &entity1,
                                 MBEntityHandle &entity2,
                                 double reference_angle_cosine)
{
  // compare normals to get angle.  We assume that the surface quads
  // which we test here will be approximately planar

  double norm[2][3];
  MBUtil::normal(thisMB, entity1, norm[0][0], norm[0][1], norm[0][2]);
  MBUtil::normal(thisMB, entity2, norm[1][0], norm[1][1], norm[1][2]);

  double cosine = norm[0][0] * norm[1][0] + norm[0][1] * norm[1][1] + norm[0][2] * norm[1][2];

  if (cosine < reference_angle_cosine)
  {
    return true;
  }


  return false;
}

  // get skin entities of prescribed dimension
MBErrorCode MBSkinner::find_skin(const MBRange &entities,
                                 int dim,
                                 MBRange &skin_entities,
                                 bool create_vert_elem_adjs) 
{
  MBRange tmp_skin;
  MBErrorCode result = find_skin(entities, (dim==0), tmp_skin, 0, 
                                 create_vert_elem_adjs, true);
  if (MB_SUCCESS != result || tmp_skin.empty()) return result;
  
  if (tmp_skin.all_of_dimension(dim)) {
    if (skin_entities.empty())
      skin_entities.swap(tmp_skin);
    else
      skin_entities.merge(tmp_skin);
  }
  else {
    result = thisMB->get_adjacencies( tmp_skin, dim, true, skin_entities, 
                                      MBInterface::UNION );
  }
  
  return result;
}

MBErrorCode MBSkinner::find_skin_vertices( const MBRange& entities,
                                           MBRange* skin_verts,
                                           MBRange* skin_elems,
                                           MBRange* skin_rev_elems,
                                           bool create_skin_elems,
                                           bool corners_only )
{
  MBErrorCode rval;
  if (entities.empty())
    return MB_SUCCESS;
  
  const int dim = MBCN::Dimension(TYPE_FROM_HANDLE(entities.front()));
  if (dim < 1 || dim > 3 || !entities.all_of_dimension(dim))
    return MB_TYPE_OUT_OF_RANGE;
  
    // are we skinning all entities
  size_t count = entities.size();
  int num_total;
  rval = thisMB->get_number_entities_by_dimension( 0, dim, num_total );
  if (MB_SUCCESS != rval)
    return rval;
  bool all = (count == (size_t)num_total);
  
    // Create a bit tag for fast intersection with input entities range. 
    // If we're skinning all the entities in the mesh, we really don't
    // need the tag.  To save memory, just create it with a default value
    // of one and don't set it.  That way MOAB will return 1 for all 
    // entities.
  MBTag tag;
  char bit = all ? 1 : 0;
  rval = thisMB->tag_create( NULL, 1, MB_TAG_BIT, tag, &bit );
  if (MB_SUCCESS != rval)
    return rval;
  
    // tag all entities in input range
  if (!all) {
    std::vector<unsigned char> vect(count, 1);
    rval = thisMB->tag_set_data( tag, entities, &vect[0] );
    if (MB_SUCCESS != rval) {
      thisMB->tag_delete(tag);
      return rval;
    }
  }
  
  switch (dim) {
    case 1:
      if (skin_verts)
        rval = find_skin_vertices_1D( tag, entities, *skin_verts );
      else if (skin_elems)
        rval = find_skin_vertices_1D( tag, entities, *skin_elems );
      else
        rval = MB_SUCCESS;
      break;
    case 2:
      rval = find_skin_vertices_2D( tag, entities, skin_verts, 
                                    skin_elems, skin_rev_elems, 
                                    create_skin_elems, corners_only );
      break;
    case 3:
      rval = find_skin_vertices_3D( tag, entities, skin_verts, 
                                    skin_elems, skin_rev_elems, 
                                    create_skin_elems, corners_only );
      break;
    default:
      rval = MB_TYPE_OUT_OF_RANGE;
      break;
  }
  
  thisMB->tag_delete(tag);
  return rval;
}

MBErrorCode MBSkinner::find_skin_vertices_1D( MBTag tag,
                                              const MBRange& edges,
                                              MBRange& skin_verts )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle>::iterator i;
  MBRange::iterator hint = skin_verts.begin();
  
  if (!edges.all_of_dimension(1))
    return MB_TYPE_OUT_OF_RANGE;
  
    // get all the vertices
  MBRange verts;
  rval = thisMB->get_connectivity( edges, verts, true );
  if (MB_SUCCESS != rval)
    return rval;
  
  std::vector<char> tag_vals;
  std::vector<MBEntityHandle> adj;
  int count;
  for (MBRange::const_iterator it = verts.begin(); it != verts.end(); ++it) {
    adj.clear();
    rval = thisMB->get_adjacencies( &*it, 1, 1, false, adj );
    if (MB_SUCCESS != rval) return rval;
    if (adj.empty())
      continue;

        // remove those not in the input list
    tag_vals.resize( adj.size() );
    rval = thisMB->tag_get_data( tag, &adj[0], adj.size(), &tag_vals[0] );
    if (MB_SUCCESS != rval) return rval;
    count = std::count( tag_vals.begin(), tag_vals.end(), '\001' );
    
    if (count == 1) {
      hint = skin_verts.insert( hint, *it );
    }
  }
  
  return MB_SUCCESS;
}


// A Container for storing a list of element sides adjacent
// to a vertex.  The template parameter is the number of 
// corners for the side.
template <unsigned CORNERS> 
class AdjSides 
{
public:
  struct Side {
    MBEntityHandle handles[CORNERS-1];
    MBEntityHandle adj_elem;
    bool skin() const { return 0 != adj_elem; }
    
    Side( const MBEntityHandle* array, int idx,
          MBEntityHandle adj, unsigned short  ) 
      : adj_elem(adj) /*, elem_side(side)*/ 
    {
      switch (CORNERS) {
        default:
          assert(false);
          break;
        case 4: handles[2] = array[(idx+3)%CORNERS];
        case 3: handles[1] = array[(idx+2)%CORNERS];
        case 2: handles[0] = array[(idx+1)%CORNERS];
      }
      if (CORNERS == 3 && handles[1] > handles[0])
        std::swap( handles[0], handles[1] );
      if (CORNERS == 4 && handles[2] > handles[0])
        std::swap( handles[0], handles[2] );
    }
    
    Side( const MBEntityHandle* array,  int idx,
          MBEntityHandle adj, unsigned short ,
          const short* indices ) 
      : adj_elem(adj) /*, elem_side(side)*/ 
    {
      switch (CORNERS) {
        default:
          assert(false);
          break;
        case 4: handles[2] = array[indices[(idx+3)%CORNERS]];
        case 3: handles[1] = array[indices[(idx+2)%CORNERS]];
        case 2: handles[0] = array[indices[(idx+1)%CORNERS]];
      }
      if (CORNERS == 3 && handles[1] > handles[0])
        std::swap( handles[0], handles[1] );
      if (CORNERS == 4 && handles[2] > handles[0])
        std::swap( handles[0], handles[2] );
    }
   
    bool operator==( const Side& other ) const 
    {
      switch (CORNERS) {
        default:
          assert(false);
          return false;
        case 4:
          return handles[0] == other.handles[0] 
              && handles[1] == other.handles[1]
              && handles[2] == other.handles[2];
        case 3:
          return handles[0] == other.handles[0] 
              && handles[1] == other.handles[1];
        case 2:
          return handles[0] == other.handles[0];
      }
    }
  };

private:

  std::vector<Side> data;
  size_t skin_count;
  
public:

  typedef typename std::vector<Side>::iterator iterator;
  typedef typename std::vector<Side>::const_iterator const_iterator;
  const_iterator begin() const { return data.begin(); }
  const_iterator end() const { return data.end(); }
  
  void clear() { data.clear(); skin_count = 0; }
  bool empty() const { return data.empty(); }
  
  AdjSides() : skin_count(0) {}
  
  size_t num_skin() const { return skin_count; }
  
  void insert( const MBEntityHandle* handles, int skip_idx,
               MBEntityHandle adj_elem, unsigned short elem_side )
  {
    Side side( handles, skip_idx, adj_elem, elem_side );
    iterator p = std::find( data.begin(), data.end(), side );
    if (p == data.end()) {
      data.push_back( side );
      ++skin_count;
    }
    else if (p->adj_elem) {
      p->adj_elem = 0;
      --skin_count;
    }
  }
  
  void insert( const MBEntityHandle* handles,  int skip_idx,
               MBEntityHandle adj_elem, unsigned short elem_side,
               const short* indices )
  {
    Side side( handles, skip_idx, adj_elem, elem_side, indices );
    iterator p = std::find( data.begin(), data.end(), side );
    if (p == data.end()) {
      data.push_back( side );
      ++skin_count;
    }
    else if (p->adj_elem) {
      p->adj_elem = 0;
      --skin_count;
    }
  }
  
  bool find_and_unmark( const MBEntityHandle* other, int skip_index, MBEntityHandle& elem_out ) 
  {
    Side s( other, skip_index, 0, 0 );
    iterator p = std::find( data.begin(), data.end(), s );
    if (p == data.end() || !p->adj_elem)
      return false;
    else {
      elem_out = p->adj_elem;
      p->adj_elem = 0;
      --skin_count;
      return true;
    }
  }
};

MBErrorCode MBSkinner::create_side( MBEntityHandle elem,
                                    MBEntityType side_type,
                                    const MBEntityHandle* side_conn,
                                    MBEntityHandle& side_elem )
{
  const int max_side = 9;
  const MBEntityHandle* conn;
  int len, side_len, side, sense, offset, indices[max_side];
  MBErrorCode rval;
  MBEntityType type = TYPE_FROM_HANDLE(elem), tmp_type;
  const int ncorner = MBCN::VerticesPerEntity( side_type );
  const int d = MBCN::Dimension(side_type);
  
  rval = thisMB->get_connectivity( elem, conn, len, false );
  if (MB_SUCCESS != rval) return rval;
 
  MBCN::SideNumber( type, conn, side_conn, ncorner, d, side, sense, offset );
  MBCN::SubEntityNodeIndices( type, len, d, side, tmp_type, side_len, indices );
  assert(side_len <= max_side);
  assert(side_type == tmp_type);
  
  //NOTE: re-create conn array even when no higher-order nodes
  //      because we want it to always be forward with respect
  //      to the side ordering.
  MBEntityHandle side_conn_full[max_side];
  for (int i = 0; i < side_len; ++i)
    side_conn_full[i] = conn[indices[i]];
  
  return thisMB->create_element( side_type, side_conn_full, side_len, side_elem );
}

bool MBSkinner::edge_reversed( MBEntityHandle face,
                               const MBEntityHandle* edge_ends )
{
  const MBEntityHandle* conn;
  int len, idx;
  MBErrorCode rval = thisMB->get_connectivity( face, conn, len, true );
  if (MB_SUCCESS != rval) {
    assert(false);
    return false;
  }
  idx = std::find( conn, conn+len, edge_ends[0] ) - conn;
  if (idx == len) {
    assert(false);
    return false;
  }
  return (edge_ends[1] == conn[(idx+len-1)%len]);
}

bool MBSkinner::face_reversed( MBEntityHandle region,
                               const MBEntityHandle* face_corners,
                               MBEntityType face_type )
{
  const MBEntityHandle* conn;
  int len, side, sense, offset;
  MBErrorCode rval = thisMB->get_connectivity( region, conn, len, true );
  if (MB_SUCCESS != rval) {
    assert(false);
    return false;
  }
  short r = MBCN::SideNumber( TYPE_FROM_HANDLE(region), conn, face_corners, 
                              MBCN::VerticesPerEntity(face_type),
                              MBCN::Dimension(face_type),
                              side, sense, offset );
  assert(0 == r);
  return (!r && sense == -1);
}

MBErrorCode MBSkinner::find_skin_vertices_2D( MBTag tag,
                                              const MBRange& faces,
                                              MBRange* skin_verts,
                                              MBRange* skin_edges,
                                              MBRange* reversed_edges,
                                              bool create_edges,
                                              bool corners_only )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle>::iterator i, j;
  MBRange::iterator hint;
  if (skin_verts)
    hint = skin_verts->begin();
  std::vector<MBEntityHandle> storage;
  const MBEntityHandle *conn;
  int len;
  bool find_edges = skin_edges || create_edges;
  MBEntityHandle face;
  
  if (!faces.all_of_dimension(2))
    return MB_TYPE_OUT_OF_RANGE;
  
    // get all the vertices
  MBRange verts;
  rval = thisMB->get_connectivity( faces, verts, true );
  if (MB_SUCCESS != rval)
    return rval;
  
  std::vector<char> tag_vals;
  std::vector<MBEntityHandle> adj;
  AdjSides<2> adj_edges;
  for (MBRange::const_iterator it = verts.begin(); it != verts.end(); ++it) {
    bool higher_order = false;
  
      // get all adjacent faces
    adj.clear();
    rval = thisMB->get_adjacencies( &*it, 1, 2, false, adj );
    if (MB_SUCCESS != rval) return rval;
    if (adj.empty())
      continue;

      // remove those not in the input list
    i = j = adj.begin();
    tag_vals.resize( adj.size() );
    rval = thisMB->tag_get_data( tag, &adj[0], adj.size(), &tag_vals[0] );
    if (MB_SUCCESS != rval) return rval;

    i = j = adj.begin();
    for (; i != adj.end(); ++i) 
      if (tag_vals[i - adj.begin()])
        *(j++) = *i;
    adj.erase( j, adj.end() );
    
      // For each adjacent face, check the edges adjacent to the current vertex
    adj_edges.clear();    // other vertex for adjacent edges
    for (i = adj.begin(); i != adj.end(); ++i) {
      rval = thisMB->get_connectivity( *i, conn, len, false, &storage );
      if (MB_SUCCESS != rval) return rval;
      
      MBEntityHandle prev, next; // vertices of two adjacent edge-sides
      const int idx = std::find(conn, conn+len, *it) - conn;
      assert(idx != len);

      if (TYPE_FROM_HANDLE(*i) == MBTRI && len > 3) {
        len = 3;
        higher_order = true;
        if (idx > 2) // skip higher-order nodes for now
          continue;
      }
      else if (TYPE_FROM_HANDLE(*i) == MBQUAD && len > 4) {
        len = 4;
        higher_order = true;
        if (idx > 3) // skip higher-order nodes for now
          continue;
      }

      const int prev_idx = (idx + len - 1)%len;
      prev = conn[prev_idx];
      next = conn[(idx+1)%len];
      adj_edges.insert( &prev, 1, *i, prev_idx );
      adj_edges.insert( &next, 1, *i, idx );
    }
    
      // If vertex is not on skin, advance to next vertex
    if (0 == adj_edges.num_skin())
      continue;
    
      // Put skin vertex in output list
    if (skin_verts) {
      hint = skin_verts->insert( hint, *it );
 
        // Add mid edge nodes to vertex list
      if (!corners_only && higher_order) {
        for (AdjSides<2>::const_iterator p = adj_edges.begin(); p != adj_edges.end(); ++p) {
          if (p->skin()) {
            face = p->adj_elem;
            MBEntityType type = TYPE_FROM_HANDLE(face);

            rval = thisMB->get_connectivity( face, conn, len, false );
            if (MB_SUCCESS != rval) return rval;
            if (!MBCN::HasMidEdgeNodes( type, len ))
              continue;

            MBEntityHandle ec[2] = { *it, p->handles[0] };
            int side, sense, offset;
            MBCN::SideNumber( type, conn, ec, 2, 1, side, sense, offset );
            offset = MBCN::HONodeIndex( type, len, 1, side );
            assert(offset >= 0 && offset < len);
            skin_verts->insert( conn[offset] );
          }
        }
      }
    }
 
    if (find_edges) {
        // Search list of adjacent edges for any that are on the skin
      adj.clear();
      rval = thisMB->get_adjacencies( &*it, 1, 1, false, adj );
      if (MB_SUCCESS != rval) return rval;
      for (i = adj.begin(); i != adj.end(); ++i) {
        rval = thisMB->get_connectivity( *i, conn, len, true );
        if (MB_SUCCESS != rval) return rval;

          // bool equalality expression will be evaluate to the 
          // index of *it in the conn array.  Note that the order
          // of the terms in the if statement is important.  We want
          // to unmark any existing skin edges even if we aren't returning
          // them.  Otherwise we'll end up creating duplicates if create_edges
          // is true.
        if (adj_edges.find_and_unmark( conn, (conn[1] == *it), face ) && skin_edges) {
          if (reversed_edges && edge_reversed( face, conn ))
            reversed_edges->insert( *i );
          else
            skin_edges->insert( *i );
        }
      }
    }
    
    if (create_edges && adj_edges.num_skin()) {
        // Create any skin edges that don't exist
      for (AdjSides<2>::const_iterator p = adj_edges.begin(); p != adj_edges.end(); ++p) {
        if (p->skin()) {
          MBEntityHandle edge, ec[] = { *it, p->handles[0] };
          rval = create_side( p->adj_elem, MBEDGE, ec, edge );
          if (MB_SUCCESS != rval) return rval;
          if (skin_edges)
            skin_edges->insert( edge );
        }
      }
    }

  } // end for each vertex
  
  return MB_SUCCESS;
}
  

MBErrorCode MBSkinner::find_skin_vertices_3D( MBTag tag,
                                              const MBRange& entities,
                                              MBRange* skin_verts,
                                              MBRange* skin_faces,
                                              MBRange* reversed_faces,
                                              bool create_faces,
                                              bool corners_only )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle>::iterator i, j;
  MBRange::iterator hint;
  if (skin_verts)
    hint = skin_verts->begin();
  std::vector<MBEntityHandle> storage, storage2;
  const MBEntityHandle *conn, *conn2;
  int len, len2;
  bool find_faces = skin_faces || create_faces;
  int clen, side, sense, offset, indices[9];
  MBEntityType face_type;
  MBEntityHandle elem;
  
  if (!entities.all_of_dimension(3))
    return MB_TYPE_OUT_OF_RANGE;
  
  MBRange verts;
  rval = thisMB->get_adjacencies( entities, 0, false, verts, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  
  AdjSides<4> adj_quads;
  AdjSides<3> adj_tris;  
  AdjSides<2> adj_poly;  
  std::vector<char> tag_vals;
  std::vector<MBEntityHandle> adj;
  for (MBRange::const_iterator it = verts.begin(); it != verts.end(); ++it) {
    bool higher_order = false;
  
      // get all adjacent elements
    adj.clear();
    rval = thisMB->get_adjacencies( &*it, 1, 3, false, adj );
    if (MB_SUCCESS != rval) return rval;
    if (adj.empty())
      continue;
      
      // remove those not in the input list
    i = j = adj.begin();
    tag_vals.resize( adj.size() );
    rval = thisMB->tag_get_data( tag, &adj[0], adj.size(), &tag_vals[0] );
    if (MB_SUCCESS != rval) return rval;
    for (; i != adj.end(); ++i) 
      if (tag_vals[i - adj.begin()])
        *(j++) = *i;
    adj.erase( j, adj.end() );
      
      // Build lists of sides of 3D element adjacent to the current vertex
    adj_quads.clear(); // store three other vertices for each adjacent quad face
    adj_tris.clear();  // store two other vertices for each adjacent tri face
    adj_poly.clear();  // store handle of each adjacent polygonal face
    int idx;
    for (i = adj.begin(); i != adj.end(); ++i) {
      const MBEntityType type = TYPE_FROM_HANDLE(*i);
      
        // Special case for POLYHEDRA
      if (type == MBPOLYHEDRON) {
        rval = thisMB->get_connectivity( *i, conn, len );
        if (MB_SUCCESS != rval) return rval;
        for (int k = 0; k < len; ++k) {
          rval = thisMB->get_connectivity( conn[k], conn2, len2, true, &storage2 );
          if (MB_SUCCESS != rval) return rval;
          idx = std::find( conn2, conn2+len2, *it) - conn2;
          if (idx == len2) // vertex not in this face
            continue;
          
            // Treat 3- and 4-vertex faces specially, so that
            // if the mesh contains both elements and polyhedra,
            // we don't miss one type adjacent to the other.
          switch (len2) {
            case 3:
              adj_tris.insert( conn2, idx, *i, k );
              break;
            case 4:
              adj_quads.insert( conn2, idx, *i, k );
              break;
            default:
              adj_poly.insert( conn+k, 1, *i, k );
              break;
            }
        }
      }
      else {
        rval = thisMB->get_connectivity( *i, conn, len, false, &storage );
        if (MB_SUCCESS != rval) return rval;

        idx = std::find(conn, conn+len, *it) - conn;
        assert(idx != len);
        
        if (len > MBCN::VerticesPerEntity( type )) {
          higher_order =true;
            // skip higher-order nodes for now
          if (idx >= MBCN::VerticesPerEntity( type )) 
            continue;
        }

        const int num_faces = MBCN::NumSubEntities( type, 2 );
        for (int f = 0; f < num_faces; ++f) {
          int num_vtx;
          const short* face_indices = MBCN::SubEntityVertexIndices(type, 2, f, face_type, num_vtx );
          const short face_idx = std::find(face_indices, face_indices+num_vtx, (short)idx) - face_indices;
          if (face_idx == num_vtx)
            continue; // current vertex not in this face

          assert(num_vtx <= 4);
          switch (face_type) {
            case MBTRI:
              adj_tris.insert( conn, face_idx, *i, f, face_indices );
              break;
            case MBQUAD:
              adj_quads.insert( conn, face_idx, *i, f, face_indices );
              break;
            default:
              return MB_TYPE_OUT_OF_RANGE;
          }
        }
      }
    } // end for (adj[3])
    
      // If vertex is not on skin, advance to next vertex
    if (0 == (adj_tris.num_skin() + adj_quads.num_skin() + adj_poly.num_skin()))
      continue;
    
      // Put skin vertex in output list
    if (skin_verts) {
      hint = skin_verts->insert( hint, *it );
 
        // Add mid-edge and mid-face nodes to vertex list
      if (!corners_only && higher_order) {
        for (AdjSides<3>::const_iterator t = adj_tris.begin(); t != adj_tris.end(); ++t) {
          if (t->skin()) {
            elem = t->adj_elem;
            MBEntityType type = TYPE_FROM_HANDLE(elem);

            rval = thisMB->get_connectivity( elem, conn, len, false );
            if (MB_SUCCESS != rval) return rval;
            if (!MBCN::HasMidNodes( type, len ))
              continue;

            MBEntityHandle ec[3] = { *it, t->handles[0], t->handles[1] };
            MBCN::SideNumber( type, conn, ec, 3, 2, side, sense, offset );
            MBCN::SubEntityNodeIndices( type, len, 2, side, face_type, clen, indices );
            assert(MBTRI == face_type);
            for (int k = 3; k < clen; ++k)
              skin_verts->insert( conn[indices[k]] );
          }
        }
        for (AdjSides<4>::const_iterator q = adj_quads.begin(); q != adj_quads.end(); ++q) {
          if (q->skin()) {
            elem = q->adj_elem;
            MBEntityType type = TYPE_FROM_HANDLE(elem);

            rval = thisMB->get_connectivity( elem, conn, len, false );
            if (MB_SUCCESS != rval) return rval;
            if (!MBCN::HasMidNodes( type, len ))
              continue;

            MBEntityHandle ec[4] = { *it, q->handles[0], q->handles[1], q->handles[2] };
            MBCN::SideNumber( type, conn, ec, 4, 2, side, sense, offset );
            MBCN::SubEntityNodeIndices( type, len, 2, side, face_type, clen, indices );
            assert(MBQUAD == face_type);
            for (int k = 4; k < clen; ++k)
              skin_verts->insert( conn[indices[k]] );
          }
        }
      }
    }

    if (find_faces) {
        // Search list of adjacent faces for any that are on the skin
      adj.clear();
      rval = thisMB->get_adjacencies( &*it, 1, 2, false, adj );
      if (MB_SUCCESS != rval) return rval;

      for (i = adj.begin(); i != adj.end(); ++i) {
        rval = thisMB->get_connectivity( *i, conn, len, true );
        if (MB_SUCCESS != rval) return rval;
        const int idx = std::find( conn, conn+len, *it ) - conn;
        assert(idx != len);
          // Note that the order of the terms in the if statements below
          // is important.  We want to unmark any existing skin faces even 
          // if we aren't returning them.  Otherwise we'll end up creating 
          // duplicates if create_faces is true.
        if (3 == len) {
          if (adj_tris.find_and_unmark( conn, idx, elem ) && skin_faces) {
            if (reversed_faces && face_reversed( elem, conn, MBTRI ))
              reversed_faces->insert( *i );
            else
              skin_faces->insert( *i );
          }
        }
        else if (4 == len) {
          if (adj_quads.find_and_unmark( conn, idx, elem ) && skin_faces) {
            if (reversed_faces && face_reversed( elem, conn, MBQUAD ))
              reversed_faces->insert( *i );
            else
              skin_faces->insert( *i );
          }
        }
        else {
          if (adj_poly.find_and_unmark( &*i, 1, elem ) && skin_faces)
            skin_faces->insert( *i );
        }
      }
    }

    if (!create_faces)
      continue;

      // Polyhedra always have explictly defined faces, so
      // there is no way we could need to create such a face.
    assert(0 == adj_poly.num_skin());
    
      // Create any skin tris that don't exist
    if (adj_tris.num_skin()) {
      for (AdjSides<3>::const_iterator t = adj_tris.begin(); t != adj_tris.end(); ++t) {
        if (t->skin()) {
          MBEntityHandle tri, c[3] = { *it, t->handles[0], t->handles[1] };
          rval = create_side( t->adj_elem, MBTRI, c, tri );
          if (MB_SUCCESS != rval) return rval;
          if (skin_faces)
            skin_faces->insert( tri );
        }
      }
    }
    
      // Create any skin quads that don't exist
    if (adj_quads.num_skin()) {
      for (AdjSides<4>::const_iterator q = adj_quads.begin(); q != adj_quads.end(); ++q) {
        if (q->skin()) {
          MBEntityHandle quad, c[4] = { *it, q->handles[0], q->handles[1], q->handles[2] };
          rval = create_side( q->adj_elem, MBQUAD, c, quad );
          if (MB_SUCCESS != rval) return rval;
          if (skin_faces)
            skin_faces->insert( quad );
        }
      }
    }
  } // end for each vertex
  
  return MB_SUCCESS;
}
