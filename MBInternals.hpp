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


#ifndef MB_INTERNALS_HPP
#define MB_INTERNALS_HPP

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#ifndef IS_BUILDING_MB
#error "MBInternals.hpp isn't supposed to be included into an application"
#endif

#include <vector>

#include "MBInterface.hpp"

/*! Define MBEntityHandle for both 32 bit and 64 bit systems.
 *  The decision to use 64 bit handles must be made at compile time.
 *  \bug we should probably have an Int64 typedef
 *
 *  MBEntityHandle format:
 *  0xXYYYYYYY  (assuming a 32-bit handle.  Top 4 bits reserved on a 64 bit system)
 *  X - reserved for entity type.  This system can only handle 15 different types
 *  Y - Entity id space.  Max id is over 200M
 *
 *  Note that for specialized databases (such as all hex) 16 bits are not
 *  required for the entity type and the id space can be increased to over 2B.
 */
#define MB_TYPE_WIDTH 4
#define MB_PROC_AND_ID_WIDTH (8*sizeof(MBEntityHandle)-MB_TYPE_WIDTH)
#define MB_ID_WIDTH (MB_PROC_AND_ID_WIDTH-MB_PROC_WIDTH)
#define MB_TYPE_MASK ((MBEntityHandle)0xF << MB_PROC_AND_ID_WIDTH)
//             2^MB_TYPE_WIDTH-1 ------^

#define MB_START_ID 1              //!< All entity id's currently start at 1
#define MB_END_ID MB_ID_MASK //!< Last id is the complement of the MASK

extern unsigned MB_PROC_WIDTH;
extern unsigned MB_PROC_RANK;
extern unsigned MB_PROC_COUNT;
extern MBEntityHandle MB_PROC_MASK;
extern MBEntityHandle MB_ID_MASK;

/**\brief Initialize data for parallel mode
 *
 * Set data used for constructing entity handles in parallel mode.
 */
void MB_SET_PROC( int num_cpu, int rank );

//! Given a type and an id create a handle.  
inline MBEntityHandle CREATE_HANDLE(const int type, const MBEntityHandle id, const int proc, int& err) 
{
  err = 0; //< Assume that there is a real error value defined somewhere

  if (id > MB_END_ID || type > MBMAXTYPE)
  {
    err = 1;   //< Assume that there is a real error value defined somewhere
    return 1;  //<You've got to return something.  What do you return?
  }
  
  MBEntityHandle type_bits = (MBEntityHandle)type << MB_PROC_AND_ID_WIDTH;
  MBEntityHandle proc_bits = ((MBEntityHandle)proc << MB_ID_WIDTH)&MB_PROC_MASK;
  return (MBEntityHandle)(type_bits|proc_bits|id);
}

//! Given a type and an id create a handle.  
inline MBEntityHandle CREATE_HANDLE(int type, MBEntityHandle id, int& err) 
{
  return CREATE_HANDLE( type, id, MB_PROC_RANK, err );
}

//! Get the entity id out of the handle.
inline unsigned long ID_FROM_HANDLE (MBEntityHandle handle)
{
  return (handle & MB_ID_MASK);
}

//! Get the type out of the handle.  Can do a simple shift because
//! handles are unsigned (therefore shifting fills with zero's)
inline MBEntityType TYPE_FROM_HANDLE(MBEntityHandle handle) 
{
  return static_cast<MBEntityType> (handle >> MB_PROC_AND_ID_WIDTH);
}

//! Get the type out of the handle.  Can do a simple shift because
//! handles are unsigned (therefore shifting fills with zero's)
inline int PROC_FROM_HANDLE(MBEntityHandle handle) 
{
  return (handle & MB_PROC_MASK) >> MB_ID_WIDTH;
}

typedef unsigned long MBTagId;

/* MBTag format
 * 0xXXZZZZZZ  ( 32 bits total )
 * Z - reserved for internal sub-tag id
 * X - reserved for internal properties & lookup speed
 */
#define MB_TAG_PROP_WIDTH 2
#define MB_TAG_SHIFT_WIDTH (8*sizeof(MBTagId)-MB_TAG_PROP_WIDTH)
#define MB_TAG_PROP_MASK (~((MBTagId)0x3 << MB_TAG_SHIFT_WIDTH))
//        2^MB_TAG_PROP_WIDTH-1 ------^

inline unsigned long ID_FROM_TAG_HANDLE(MBTag tag_handle) 
{
  return static_cast<MBTagId>( (MBTagId)tag_handle & MB_TAG_PROP_MASK );
}

inline MBTagType PROP_FROM_TAG_HANDLE(MBTag tag_handle)
{
  return static_cast<MBTagType>( (MBTagId)tag_handle >> MB_TAG_SHIFT_WIDTH );
}

inline MBTag TAG_HANDLE_FROM_ID(MBTagId tag_id, MBTagType prop) 
{
//  assert( tag_id <= MB_TAG_PROP_MASK );
  return reinterpret_cast<MBTag>( ((MBTagId)prop)<<MB_TAG_SHIFT_WIDTH | tag_id );
}

#endif
