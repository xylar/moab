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


/**********************************************
 * Filename   :     SparseTagCollection.cpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/


#include <memory.h>
#include <algorithm>

#include "SparseTagCollection.hpp"
#include "MBRange.hpp"
#include "TagCompare.hpp"
#include "VarLenTag.hpp"

SparseTagCollection::SparseTagCollection(int data_size)
{
  mDataSize = data_size;
}

SparseTagCollection::~SparseTagCollection()
{
  void* tag_data = NULL;

  myMapType::iterator tag_iterator;
  for(tag_iterator = mData.begin(); tag_iterator != mData.end(); ++tag_iterator)
  {
    tag_data = tag_iterator->second;
    if (mDataSize == MB_VARIABLE_LENGTH)
      reinterpret_cast<VarLenTag*>(tag_data)->clear();
    if(tag_data != NULL)
      mAllocator.destroy(tag_data);
  }
  mData.clear();
}

MBErrorCode SparseTagCollection::set_data(const MBEntityHandle entity_handle, const void* data)
{
  if (mDataSize == MB_VARIABLE_LENGTH)
    return MB_VARIABLE_DATA_LENGTH;

  MBErrorCode ret_val = MB_TAG_NOT_FOUND;

#ifdef HAVE_UNORDERED_MAP
  myMapType::iterator iterator = mData.find(entity_handle);
#else
  myMapType::iterator iterator = mData.lower_bound(entity_handle);
#endif
  
  // data space already exists
  if (iterator!= mData.end() && iterator->first == entity_handle)
  {
    memcpy( iterator->second, data, mDataSize);
    ret_val = MB_SUCCESS;
  }
  // we need to make some data space
  else 
  {
    void* new_data = mAllocator.allocate(mDataSize);
    memcpy(new_data, data, mDataSize);
    mData.insert(iterator, std::pair<const MBEntityHandle,void*>(entity_handle, new_data));
    ret_val = MB_SUCCESS;
  }

  return ret_val;
}

MBErrorCode SparseTagCollection::get_data(const MBEntityHandle entity_handle, void* data)
{
  if (mDataSize == MB_VARIABLE_LENGTH)
    return MB_VARIABLE_DATA_LENGTH;

  myMapType::iterator iter = mData.find(entity_handle);

  if(iter == mData.end())
    return MB_TAG_NOT_FOUND;
  
  memcpy(data, iter->second, mDataSize);
  return MB_SUCCESS;
}

MBErrorCode SparseTagCollection::set_data(const MBEntityHandle entity_handle, const void* data, int size)
{
#ifdef HAVE_UNORDERED_MAP
  myMapType::iterator iterator = mData.find(entity_handle);
#else
  myMapType::iterator iterator = mData.lower_bound(entity_handle);
#endif
  
  if (mDataSize == MB_VARIABLE_LENGTH) {
    if (size == 0) {
        // ignore return value: still success if entity not found
      remove_data( entity_handle );
      return MB_SUCCESS;
    }
  
    if (iterator == mData.end() || iterator->first != entity_handle) {
      void* new_data = mAllocator.allocate(sizeof(VarLenTag));
      new (new_data) VarLenTag;
      iterator = mData.insert( iterator, std::pair<const MBEntityHandle,void*>(entity_handle, new_data) );
    }
    reinterpret_cast<VarLenTag*>(iterator->second)->set( data, size );
  }
  else {
    if (size != 0 && size != mDataSize)
      return MB_INVALID_SIZE;
  
    if (iterator == mData.end() || iterator->first != entity_handle) {
      void* new_data = mAllocator.allocate(mDataSize);
      iterator = mData.insert( iterator, std::pair<const MBEntityHandle,void*>(entity_handle, new_data) );
    }
    memcpy( iterator->second, data, mDataSize);
  }

  return MB_SUCCESS;
}

MBErrorCode SparseTagCollection::get_data(const MBEntityHandle entity_handle, const void*& data, int& size)
{
  myMapType::iterator iter = mData.find(entity_handle);

  if(iter == mData.end())
    return MB_TAG_NOT_FOUND;
  
  if (mDataSize == MB_VARIABLE_LENGTH) {
    const VarLenTag* vtag = reinterpret_cast<const VarLenTag*>(iter->second);
    size = vtag->size();
    data = vtag->data();
  }
  else {
    size = mDataSize;
    data = iter->second;
  }
  return MB_SUCCESS;
}



MBErrorCode SparseTagCollection::remove_data( const MBEntityHandle entity_handle )
{
  myMapType::iterator iterator = mData.find(entity_handle);

  if(iterator != mData.end())
  {
    if (mDataSize == MB_VARIABLE_LENGTH)
      reinterpret_cast<VarLenTag*>(iterator->second)->clear();
    mAllocator.destroy(iterator->second);
    mData.erase(iterator);
    return MB_SUCCESS;
  }
  return MB_ENTITY_NOT_FOUND;
}


//! get number of entities of type
MBErrorCode SparseTagCollection::get_number_entities(MBEntityType type, int& num_entities)
{
  num_entities = 0;
  myMapType::iterator iter;
  for(iter = mData.begin(); iter != mData.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(iter->first) == type)
      num_entities++;
  }
  return MB_SUCCESS;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagCollection::get_entities(MBEntityType type, MBRange &entities)
{
  myMapType::iterator iter;
  for(iter = mData.begin(); iter != mData.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(iter->first) == type)
      entities.insert(iter->first);    
  }
  return MB_SUCCESS;
}


//! gets all entity handles that match a type, tag and tag value
MBErrorCode SparseTagCollection::get_entities_with_tag_value(
                                                    const TagInfo& tag_info,
                                                    MBEntityType type, 
                                                    MBRange &entities, 
                                                    const void* tag_value,
                                                    int value_size)
{
  myMapType::iterator iter, end;
#ifdef HAVE_UNORDERED_MAP
  iter = mData.begin();
  end = mData.end();
  while(iter != end) {
    while (iter != end && TYPE_FROM_HANDLE(iter->first) != type)
      ++iter;
    myMapType::iterator iter2 = iter;
    while (iter2 != end && TYPE_FROM_HANDLE(iter2->first) == type)
      ++iter2;
    find_tag_values_equal( tag_info, tag_value, value_size, iter, iter2, entities );
    iter = iter2;
  }
    
#else
  int junk;
  iter = mData.lower_bound( CREATE_HANDLE( type, MB_START_ID, junk ) );
  end = mData.upper_bound( CREATE_HANDLE( type, MB_END_ID, junk ) );
  find_tag_values_equal( tag_info, tag_value, value_size, iter, end, entities );
#endif
  return MB_SUCCESS;
}




