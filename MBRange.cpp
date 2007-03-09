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

/****************************************************
 * File     :      MBRange.cpp
 *
 * Purpose  :      Stores contiguous or partially
 *                 contiguous values in an optimized
 *                 fashion.  Partially contiguous
 *                 accessing patterns is also optimized.
 *
 * Creator  :      Clinton Stimpson
 *
 * Date     :      15 April 2002
 *
 *******************************************************/


#include <assert.h>
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "MBCN.hpp"
#include <iostream>

#ifdef HAVE_BOOST_POOL_SINGLETON_POOL_HPP
#  include <boost/pool/singleton_pool.hpp>
   typedef boost::singleton_pool< MBRange::PairNode , sizeof(MBRange::PairNode) > 
    PairAlloc;
   static inline MBRange::PairNode* alloc_pair()
    { return new (PairAlloc::malloc()) MBRange::PairNode; }
   static inline MBRange::PairNode* alloc_pair(MBRange::PairNode* n, MBRange::PairNode* p, MBEntityHandle f, MBEntityHandle s)
    { return new (PairAlloc::malloc()) MBRange::PairNode(n,p,f,s); }
   static inline void free_pair( MBRange::PairNode* node )
    { node->~PairNode(); PairAlloc::free( node ); }
#else
   static inline MBRange::PairNode* alloc_pair()
    { return new MBRange::PairNode; }
   static inline MBRange::PairNode* alloc_pair(MBRange::PairNode* n, MBRange::PairNode* p, MBEntityHandle f, MBEntityHandle s)
    { return new MBRange::PairNode(n,p,f,s); }
   static inline void free_pair( MBRange::PairNode* node )
    { delete node; }
#endif

/*! 
  returns the number of values this list represents
 */
MBEntityID MBRange::size() const
{
  // go through each pair and add up the number of values
  // we have.
  MBEntityID size=0;
  for(PairNode* iter = mHead.mNext; iter != &mHead; iter = iter->mNext)
  {
    size += ((iter->second - iter->first) + 1);
  }
  return size;
}

/*!
  advance iterator
*/
MBRange::const_iterator& MBRange::const_iterator::operator+=( MBEntityID sstep )
{
    // Check negative now to avoid infinite loop below.
  if (sstep < 0)
  {
    return operator-=( -sstep );
  }
  MBEntityHandle step = sstep;
  
    // Handle current PairNode.  Either step is within the current
    // node or need to remove the remainder of the current node
    // from step.
  MBEntityHandle this_node_rem = mNode->second - mValue;
  if (this_node_rem >= step)
  {
    mValue += step;
    return *this;
  }
  step -= this_node_rem + 1;

    // For each node we are stepping past, decrement step
    // by the size of the node.
  PairNode* node = mNode->mNext;
  MBEntityHandle node_size = node->second - node->first + 1;
  while (step >= node_size)
  {
    step -= node_size;
    node = node->mNext;
    node_size = node->second - node->first + 1;
  }
  
    // Advance into the resulting node by whatever is
    // left in step.
  mNode = node;
  mValue = mNode->first + step;
  return *this;
}
  
    
 
/*!
  regress iterator
*/
MBRange::const_iterator& MBRange::const_iterator::operator-=( MBEntityID sstep )
{
    // Check negative now to avoid infinite loop below.
  if (sstep < 0)
  {
    return operator+=( -sstep );
  }
  MBEntityHandle step = sstep;
  
    // Handle current PairNode.  Either step is within the current
    // node or need to remove the remainder of the current node
    // from step.
  MBEntityHandle this_node_rem = mValue - mNode->first;
  if (this_node_rem >= step)
  {
    mValue -= step;
    return *this;
  }
  step -= this_node_rem + 1;

    // For each node we are stepping past, decrement step
    // by the size of the node.
  PairNode* node = mNode->mPrev;
  MBEntityHandle node_size = node->second - node->first + 1;
  while (step >= node_size)
  {
    step -= node_size;
    node = node->mPrev;
    node_size = node->second - node->first + 1;
  }
  
    // Advance into the resulting node by whatever is
    // left in step.
  mNode = node;
  mValue = mNode->second - step;
  return *this;
}
  
 
  

  //! another constructor that takes an initial range
MBRange::MBRange( MBEntityHandle val1, MBEntityHandle val2 )
{
  mHead.mNext = mHead.mPrev = alloc_pair(&mHead, &mHead, val1, val2);
  mHead.first = mHead.second = 0;
}

  //! copy constructor
MBRange::MBRange(const MBRange& copy)
{
    // set the head node to point to itself
  mHead.mNext = mHead.mPrev = &mHead;
  mHead.first = mHead.second = 0;

  const PairNode* copy_node = copy.mHead.mNext;
  PairNode* new_node = &mHead;
  for(; copy_node != &(copy.mHead); copy_node = copy_node->mNext)
  {
    PairNode* tmp_node = alloc_pair(new_node->mNext, new_node, copy_node->first,
                                      copy_node->second);
    new_node->mNext->mPrev = tmp_node;
    new_node->mNext = tmp_node;
    new_node = tmp_node;
  }
}

  //! clears the contents of the list 
void MBRange::clear()
{
  PairNode* tmp_node = mHead.mNext;
  while(tmp_node != &mHead)
  {
    PairNode* to_delete = tmp_node;
    tmp_node = tmp_node->mNext;
    free_pair( to_delete );
  }
  mHead.mNext = &mHead;
  mHead.mPrev = &mHead;
}

MBRange& MBRange::operator=(const MBRange& copy)
{
  clear();
  const PairNode* copy_node = copy.mHead.mNext;
  PairNode* new_node = &mHead;
  for(; copy_node != &(copy.mHead); copy_node = copy_node->mNext)
  {
    PairNode* tmp_node = alloc_pair(new_node->mNext, new_node, copy_node->first,
                                      copy_node->second);
    new_node->mNext->mPrev = tmp_node;
    new_node->mNext = tmp_node;
    new_node = tmp_node;
  }
  return *this;
}


/*!
  inserts a single value into this range
*/

MBRange::iterator MBRange::insert(MBEntityHandle val)
{

  // if this is empty, just add it and return an iterator to it
  if(val == 0)
    return end();

  if(&mHead == mHead.mNext)
  {
    mHead.mNext = mHead.mPrev = alloc_pair(&mHead, &mHead, val, val);
    return iterator(mHead.mNext, val);
  }
  
  // find the location in the list where we can safely insert
  // new items and keep it ordered
  PairNode* jter = mHead.mNext;
  for( ; (jter != &mHead) && (jter->second < val); jter=jter->mNext);
  PairNode* iter = jter;
  jter = jter->mPrev;

  // if this val is already in the list
  if( (iter->first <= val && iter->second >= val) && (iter != &mHead) )
  {
    // return an iterator pointing to this location
    return iterator( iter, val );
  }

  // one of a few things can happen at this point
  // 1. this range needs to be backwardly extended
  // 2. the previous range needs to be forwardly extended
  // 3. a new range needs to be added

  
  // extend this range back a bit
  else if( (iter->first == (val+1)) && (iter != &mHead) )
  {
    iter->first = val;
    // see if we need to merge two ranges
    if( (iter != mHead.mNext) && (jter->second == (val-1)))
    {
      jter->second = iter->second;
      iter->mPrev->mNext = iter->mNext;
      iter->mNext->mPrev = iter->mPrev;
      free_pair( iter );
      return iterator( jter, val );
    }
    else
    {
      return iterator( iter, val );
    }

  }
  // extend the previous range forward a bit
  else if( (jter->second == (val-1)) && (iter != mHead.mNext) )
  {
    jter->second = val;
    return iterator(jter, val);
  }
  // make a new range
  else
  {
    PairNode* new_node = alloc_pair(iter, iter->mPrev, val, val);
    iter->mPrev = new_node->mPrev->mNext = new_node;
    return iterator(new_node, val);
  }

}


/*!
  inserts a range of values
*/
MBRange::iterator MBRange::insert(MBEntityHandle val1, MBEntityHandle val2)
{

  if(val1 == 0 || val1 > val2)
    return end();

  return insert( begin(), val1, val2 );
}

MBRange::iterator MBRange::insert( MBRange::iterator prev,
                                   MBEntityHandle val1, 
                                   MBEntityHandle val2 )
{
  assert( val1 <= val2 && val1 );

  // Empty 
  if (mHead.mNext == &mHead)
  {
    assert( prev == end() );
    PairNode* new_node = alloc_pair( &mHead, &mHead, val1, val2 );
    mHead.mNext = mHead.mPrev = new_node;
    return iterator( mHead.mNext, val1 );
  }
  
  PairNode* iter = prev.mNode;
  assert( iter != &mHead );
  assert( iter->mPrev == &mHead || iter->mPrev->second+1 < val1 );
  
  // Input range is before beginning?
  if (iter->mPrev == &mHead && val2 < iter->first - 1)
  {
    PairNode* new_node = alloc_pair( iter, &mHead,  val1, val2 );
    mHead.mNext = iter->mPrev = new_node;
    return iterator( mHead.mNext, val1 );
  }
  
  // Find first intersecting list entry, or the next entry
  // if none intersects.
  while (iter != &mHead && iter->second+1 < val1)
    iter = iter->mNext;
  
  // Need to insert new pair (don't intersect any existing pair)?
  if (iter == &mHead || iter->first-1 > val2)
  {
    PairNode* new_node = alloc_pair( iter, iter->mPrev, val1, val2 );
    iter->mPrev = iter->mPrev->mNext = new_node;
    return iterator( iter->mPrev, val1 );
  }
  
  // Make the first intersecting pair the union of itself with [val1,val2]
  if (iter->first > val1)
    iter->first = val1;
  if (iter->second >= val2)  
    return iterator( iter, val1 );
  iter->second = val2;
  
  // Merge any remaining pairs that intersect [val1,val2]
  while (iter->mNext != &mHead && iter->mNext->first <= val2 + 1)
  {
    PairNode* dead = iter->mNext;
    iter->mNext = dead->mNext;
    dead->mNext->mPrev = iter;
    
    if (dead->second > val2)
      iter->second = dead->second;
    free_pair( dead );
  }
  
  return iterator( iter, val1 );
}
    

/*!
  erases an item from this list and returns an iterator to the next item
*/

MBRange::iterator MBRange::erase(iterator iter)
{
  // one of a few things could happen
  // 1. shrink a range
  // 2. split a range
  // 3. remove a range

  if(iter == end())
    return end();

  // the iterator most likely to be returned
  iterator new_iter = iter;
  ++new_iter;

  PairNode* kter = iter.mNode;
  
  // just remove the range
  if(kter->first == kter->second)
  {
    kter->mNext->mPrev = kter->mPrev;
    kter->mPrev->mNext = kter->mNext;
    free_pair( kter );
    return new_iter;
  }
  // shrink it
  else if(kter->first == iter.mValue)
  {
    kter->first++;
    return new_iter;
  }
  // shrink it the other way
  else if(kter->second == iter.mValue)
  {
    kter->second--;
    return new_iter;
  }
  // split the range
  else
  {
    PairNode* new_node = alloc_pair(iter.mNode, iter.mNode->mPrev, kter->first, iter.mValue-1);
    new_node->mPrev->mNext = new_node->mNext->mPrev = new_node;
    iter.mNode->first = iter.mValue+1;
    return new_iter;
  }

}

void MBRange::delete_pair_node( PairNode* node )
{
  if (node != &mHead) { // pop_front() and pop_back() rely on this check
    node->mPrev->mNext = node->mNext;
    node->mNext->mPrev = node->mPrev;
    free_pair( node );
  }
}

  //! remove first entity from range
MBEntityHandle MBRange::pop_front()
{
  MBEntityHandle retval = front();
  if (mHead.mNext->first == mHead.mNext->second) // need to remove pair from range
    delete_pair_node( mHead.mNext );
  else 
    ++(mHead.mNext->first); // otherwise just adjust start value of pair

  return retval;
}

  //! remove last entity from range
MBEntityHandle MBRange::pop_back()
{
  MBEntityHandle retval = back();
  if (mHead.mPrev->first == mHead.mPrev->second) // need to remove pair from range
    delete_pair_node( mHead.mPrev );
  else
    --(mHead.mPrev->second); // otherwise just adjust end value of pair

  return retval;
}

/*!
  finds a value in the list.
  this method is preferred over other algorithms because
  it can be found faster this way.
*/
MBRange::const_iterator MBRange::find(MBEntityHandle val) const
{
  // iterator through the list
  PairNode* iter = mHead.mNext;
  for( ; iter != &mHead && (val > iter->second); iter=iter->mNext );
  return ((iter->second >= val) && (iter->first <= val)) ? const_iterator(iter,val) : end();
}

/*!
  merges another MBRange with this one
*/


void MBRange::merge( const MBRange& range )
{
  merge( range.begin(), range.end() );
}

void MBRange::merge( MBRange::const_iterator begin,
                     MBRange::const_iterator end )
{
  if (begin == end)
    return;
  
  PairNode* node = begin.mNode;
  if (end.mNode == node)
  {
    insert( *begin, (*end)-1 );
    return;
  }
  
  MBRange::iterator hint = insert( *begin, node->second );
  node = node->mNext;
  while (node != end.mNode)
  {
    hint = insert( hint, node->first, node->second );
    node = node->mNext;
  }
  
  if (*end > node->first)
  {
    if (*end <= node->second)
      insert( hint, node->first, *(end) - 1 );
    else
      insert( hint, node->first, node->second );
  }
}

  
  

#include <algorithm>


// checks the range to make sure everything is A-Ok.
void MBRange::sanity_check() const
{
  if(empty())
    return;

  const PairNode* node = mHead.mNext;
  std::vector<const PairNode*> seen_before;
  bool stop_it = false;
  
  for(; stop_it == false; node = node->mNext)
  {
    // have we seen this node before?
    assert(std::find(seen_before.begin(), seen_before.end(), node) == seen_before.end());
    seen_before.push_back(node);

    // is the connection correct?
    assert(node->mNext->mPrev == node);

    // are the values right?
    assert(node->first <= node->second);
    if(node != &mHead && node->mPrev != &mHead)
      assert(node->mPrev->second < node->first);

    if(node == &mHead)
      stop_it = true;

  }

}

// for debugging
void MBRange::print() const
{
  print( std::cout );
}

void MBRange::print( std::ostream& stream ) const
{
  if (empty()) {
    stream << "\tempty" << std::endl;
    return;
  }
  
  for (const_pair_iterator i = const_pair_begin(); i != const_pair_end(); ++i) {
    MBEntityType t1 = TYPE_FROM_HANDLE( i->first );
    MBEntityType t2 = TYPE_FROM_HANDLE( i->second );
  
    stream << "\t" << MBCN::EntityTypeName( t1 ) << " " << ID_FROM_HANDLE( i->first );
    if(i->first != i->second) {
      stream << " - ";
      if (t1 != t2) 
        stream << MBCN::EntityTypeName( t2 ) << " ";
      stream << ID_FROM_HANDLE( i->second );
    }
    stream << std::endl;
  }
}

  // intersect two ranges, placing the results in the return range
#define MAX(a,b) (a < b ? b : a)
#define MIN(a,b) (a > b ? b : a)
MBRange MBRange::intersect(const MBRange &range2) const 
{
  pair_iterator r_it[2] = {pair_iterator(begin()), pair_iterator(range2.begin())};
  MBEntityHandle low_it, high_it;
  
  MBRange lhs;
  
    // terminate the while loop when at least one "start" iterator is at the
    // end of the list
  while (r_it[0] != end() && r_it[1] != range2.end()) {
    
    if (r_it[0]->second < r_it[1]->first)
        // 1st subrange completely below 2nd subrange
      r_it[0]++;
    else if (r_it[1]->second < r_it[0]->first) 
        // 2nd subrange completely below 1st subrange
      r_it[1]++;
    
    else {
        // else ranges overlap; first find greater start and lesser end
      low_it = MAX(r_it[0]->first, r_it[1]->first);
      high_it = MIN(r_it[0]->second, r_it[1]->second);
      
        // insert into result
      lhs.insert(low_it, high_it);
      
        // now find bounds of this insertion and increment corresponding iterator
      if (high_it == r_it[0]->second) r_it[0]++;
      if (high_it == r_it[1]->second) r_it[1]++;
    }
  }
  
  return lhs;
}

MBRange MBRange::subtract(const MBRange &range2) const 
{
    // brain-dead implementation right now
  MBRange res = *this;
  for (MBRange::const_iterator rit = range2.begin(); rit != range2.end(); rit++)
    res.erase(*rit);

  return res;
}

MBRange::const_iterator MBRange::lower_bound(MBRange::const_iterator first,
                                             MBRange::const_iterator last,
                                             MBEntityHandle val)
{
    // Find the first pair whose end is >= val
  PairNode* iter;
  for (iter = first.mNode; iter != last.mNode; iter = iter->mNext)
  {
    if (iter->second >= val)
    {
        // This is the correct pair.  Either 'val' is in the range, or
        // the range starts before 'val' and iter->first IS the lower_bound.
      if (iter->first > val)
        return const_iterator(iter, iter->first);
      return const_iterator(iter, val);
    }
  }
  
  if (iter->first >= val)
    return const_iterator( iter, iter->first );
  else if(*last > val)
    return const_iterator( iter, val );
  else
    return last;
}

MBRange::const_iterator MBRange::upper_bound(MBRange::const_iterator first,
                                             MBRange::const_iterator last,
                                             MBEntityHandle val)
{
  MBRange::const_iterator result = lower_bound( first, last, val );
  if (result != last && *result == val)
    ++result;
  return result;
}

MBRange::const_iterator MBRange::lower_bound( MBEntityType type ) const
{
  int err;
  MBEntityHandle handle = CREATE_HANDLE( type, 0, err );
  return err ? end() : lower_bound( begin(), end(), handle );
}
MBRange::const_iterator MBRange::upper_bound( MBEntityType type ) const
{
    // if (type+1) overflows, err will be true and we return end().
  int err; 
  MBEntityHandle handle = CREATE_HANDLE( type + 1, 0, err );
  return err ? end() : lower_bound( begin(), end(), handle );
}
std::pair<MBRange::const_iterator, MBRange::const_iterator>
MBRange::equal_range( MBEntityType type ) const
{
  std::pair<MBRange::const_iterator, MBRange::const_iterator> result;
  int err;
  MBEntityHandle handle = CREATE_HANDLE( type, 0, err );
  result.first = err ? end() : lower_bound( begin(), end(), handle );
    // if (type+1) overflows, err will be true and we return end().
  handle = CREATE_HANDLE( type+1, 0, err );
  result.second = err ? end() : lower_bound( result.first, end(), handle );
  return result;
}
  
bool MBRange::all_of_type( MBEntityType type ) const
{
  return empty() 
      || (TYPE_FROM_HANDLE(front()) == type
       && TYPE_FROM_HANDLE(back()) == type);
}

bool MBRange::all_of_dimension( int dimension ) const
{
  return empty() 
      || (MBCN::Dimension(TYPE_FROM_HANDLE(front())) == dimension
       && MBCN::Dimension(TYPE_FROM_HANDLE(back())) == dimension);
}

unsigned MBRange::num_of_type( MBEntityType type ) const
{
  const_pair_iterator iter = const_pair_begin();
  while(iter != const_pair_end() && TYPE_FROM_HANDLE((*iter).second) < type)
    ++iter;
  
  unsigned count = 0;
  for ( ; iter != const_pair_end(); ++iter)
  {
    MBEntityType start_type = TYPE_FROM_HANDLE((*iter).first);
    MBEntityType end_type = TYPE_FROM_HANDLE((*iter).second);
    if (start_type > type)
      break;
   
    int sid = start_type < type ? 1 : ID_FROM_HANDLE((*iter).first);
    int eid = end_type > type ? MB_END_ID : ID_FROM_HANDLE((*iter).second);
    count += eid - sid + 1;
  }

  return count;
}
  
unsigned MBRange::num_of_dimension( int dim ) const
{
  const_pair_iterator iter = const_pair_begin();
  while(iter != const_pair_end() && MBCN::Dimension(TYPE_FROM_HANDLE((*iter).second)) < dim)
    ++iter;
  
  int junk;
  unsigned count = 0;
  for ( ; iter != const_pair_end(); ++iter)
  {
    int start_dim = MBCN::Dimension(TYPE_FROM_HANDLE((*iter).first));
    int end_dim = MBCN::Dimension(TYPE_FROM_HANDLE((*iter).second));
    if (start_dim > dim)
      break;
      
    MBEntityHandle sh = start_dim < dim ? 
                        CREATE_HANDLE( MBCN::TypeDimensionMap[dim].first, 1, junk ) :
                        (*iter).first;
    MBEntityHandle eh = end_dim > dim ?
                        CREATE_HANDLE( MBCN::TypeDimensionMap[dim].second, MB_END_ID, junk ) :
                        (*iter).second;
    count += eh - sh + 1;
  }

  return count;
}
  
    


//! swap the contents of this range with another one
//! THIS FUNCTION MUST NOT BE INLINED, THAT WILL ELIMINATE RANGE_EMPTY AND THIS_EMPTY
//! BY SUBSTITUTION AND THE FUNCTION WON'T WORK RIGHT!
void MBRange::swap( MBRange &range )
{
    // update next/prev nodes of head of both ranges
  bool range_empty = (range.mHead.mNext == &(range.mHead));
  bool this_empty = (mHead.mNext == &mHead);

  range.mHead.mNext->mPrev = (range_empty ? &(range.mHead) : &mHead);
  range.mHead.mPrev->mNext = (range_empty ? &(range.mHead) : &mHead);
  mHead.mNext->mPrev = (this_empty ? &mHead : &(range.mHead));
  mHead.mPrev->mNext = (this_empty ? &mHead : &(range.mHead));

    // switch data in head nodes of both ranges
  PairNode *range_next = range.mHead.mNext, *range_prev = range.mHead.mPrev;
  range.mHead.mNext = (this_empty ? &(range.mHead) : mHead.mNext);
  range.mHead.mPrev = (this_empty ? &(range.mHead) : mHead.mPrev);
  mHead.mNext = (range_empty ? &mHead : range_next);
  mHead.mPrev = (range_empty ? &mHead : range_prev);

}

MBEntityHandle MBRange::operator[](MBEntityID index) const
{
  MBRange::const_iterator i = begin();
  i += index;
  return *i;
}

    //! return a subset of this range, by type
MBRange MBRange::subset(const MBEntityType t) 
{
  MBRange result;
  result.merge( lower_bound(t), upper_bound(t) );
  return result;
}


bool operator==( const MBRange& r1, const MBRange& r2 )
{
  MBRange::const_pair_iterator i1, i2;
  i1 = r1.const_pair_begin();
  i2 = r2.const_pair_begin();
  for ( ; i1 != r1.const_pair_end(); ++i1, ++i2) 
    if (i2 == r2.const_pair_end() ||
        i1->first != i2->first ||
        i1->second != i2->second)
      return false;
  return i2 == r2.const_pair_end();
}

