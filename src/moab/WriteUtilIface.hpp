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


#ifndef MOAB_WRITE_UTIL_IFACE_HPP
#define MOAB_WRITE_UTIL_IFACE_HPP


#include "moab/Range.hpp"
#include <vector>

namespace moab {

//! Interface implemented in MOAB which provides memory for mesh reading utilities
class MB_DLL_EXPORT WriteUtilIface
{
public:

    //! constructor
  WriteUtilIface(){}

    //! destructor
  virtual ~WriteUtilIface(){}
  
    //! Check if the specified file already exists.
    //! Returns MB_SUCCESS if file does not exist, MB_ALREADY_ALLOCATED
    //! if file does exist, or MB_FAILURE for some other error condition.
  virtual ErrorCode check_doesnt_exist( const char* file_name ) = 0;

    //! Given information about the nodes to be written, and pointers to memory
    //! to which coordinates will be written, writes coordinate data there, and
    //! also assigns global ids to nodes & writes to a tag, if a tag is specified
    //! \param num_arrays Number of coordinate arrays requested
    //! \param num_nodes Number of nodes to be written
    //! \param entities Range of nodes to be written
    //! \param node_id_tag Tag used to write ids to nodes
    //! \param start_node_id Starting value for node ids
    //! \param arrays Pointers to memory where coordinate data will be written
    //! \return status Return status
  virtual ErrorCode get_node_arrays(
    const int num_arrays,
    const int num_nodes, 
    const Range& entities, 
    Tag node_id_tag,
    const int start_node_id,
    std::vector<double*>& arrays
    ) = 0;
      
  /** Get an array of coordinate values for nodes
   *
   * Given a range of node handles, retreive a single coordinate
   * value for each. 
   *
   * Failure conditions:
   *  - invalid entity handles (not vertices, non-existant entity, etc.)
   *  - range is empty (<code>iter == end</code>)
   *  - <code>output_array</code> is null
   *  - insufficient space in <code>output_array</code>
   *
   *\param which_array  The coordinate to retreive (0-&gt;X, 1-&gt;Y, 2-&gt;Z)
   *\param begin        The first node handle.
   *\param end          One past the last node handle.
   *\param output_size  The size of <code>output_array</code>.
   *\param output_array The memory in which to write the node coordinates.
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_node_array(
      const int which_array, 
      Range::const_iterator begin,
      const Range::const_iterator end,
      const size_t output_size,
      double* const output_array
      ) = 0;

    //! Given information about elements to be written and a pointer to memory
    //! where connectivity for those elements should be written, writes connectivity
    //! to that memory; uses node ids stored in a tag during call to <em>get_node_arrays</em>
    //! function
    //! \param num_elements Number of elements to be written
    //! \param verts_per_element Number of vertices per element
    //! \param node_id_tag Tag used to store node ids
    //! \param entities Range of elements to be written
    //! \param element_id_tag Tag which should be used to store element ids
    //! \param start_element_id Starting value for element ids
    //! \param array Pointer to memory where connectivity data will be written
    //! \return status Return status
  virtual ErrorCode get_element_array(
    const int num_elements, 
    const int verts_per_element,
    Tag node_id_tag,
    const Range& entities, 
    Tag element_id_tag,
    int start_element_id,
    int* array
    ) = 0;

  /** Get connectivity for elements 
   *
   * Get the connectivity list for a range of elements.
   *
   * Failure cases:
   *  - Passed range is empty (<code>begin == end</code>).
   *  - <code>vertices_per_elem</code> is less than one
   *  - <code>element_array</code> is null.
   *  - The range contains invalid handles (non-existant entities,
   *      not an element, etc.)
   *  - Retreiving ID tag for an entity failed.
   *  - Insufficient space in passed array.
   *
   *\param begin        The first element handle
   *\param end          One past the last element handle
   *\param vertices_per_elem Number of vertices to retreive for each
   *                    element.  If the element has more vertices, the
   *                    element connectivity will be truncated.  If 
   *                    <code>vertices_per_elem</code> is greater than the
   *                    number of nodes for an eleement, the data will be
   *                    padded with zeros.
   *\param node_id_tag  A tag with integer values.  
   *\param array_size   The length of <code>element_array</code>
   *\param element_array The memory location at which to store the 
   *                    connectivity list.
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_element_array(
      Range::const_iterator begin,
      const Range::const_iterator end,
      const int vertices_per_elem,
      Tag node_id_tag,
      const size_t array_size, 
      int *const element_array
      ) = 0;

  /** Get connectivity for elements 
   *
   * Get the connectivity list for a range of elements.
   *
   * Failure cases:
   *  - Passed range is empty (<code>begin == end</code>).
   *  - <code>vertices_per_elem</code> is less than one
   *  - <code>element_array</code> is null.
   *  - The range contains invalid handles (non-existant entities,
   *      not an element, etc.)
   *  - Insufficient space in passed array.
   *
   *\param begin        The first element handle
   *\param end          One past the last element handle
   *\param vertices_per_elem Number of vertices to retreive for each
   *                    element.  If the element has more vertices, the
   *                    element connectivity will be truncated.  If 
   *                    <code>vertices_per_elem</code> is greater than the
   *                    number of nodes for an eleement, the data will be
   *                    padded with zeros.
   *\param array_size   The length of <code>element_array</code>
   *\param element_array The memory location at which to store the 
   *                    connectivity list.
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_element_array(
      Range::const_iterator begin,
      const Range::const_iterator end,
      const int vertices_per_elem,
      const size_t array_size, 
      EntityHandle *const element_array
      ) = 0;


  /** Get poly (polygon or polyhedron) connectivity size
   *\param begin  First iterator in range of poly
   *\param end    One past last in range of poly.
   *\param connectivity_size  The lenght of the connectivity list
   *              For the specified range of polyhedra.
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_poly_array_size(
      Range::const_iterator begin,
      const Range::const_iterator end,
      int& connectivity_size 
      ) = 0;
   

  /** Get poly (polygon or polyhedron) connectivity.
   *
   * Connectivity is returned in two arrays.  The first is
   * an array of global IDs that is the concatenation of the
   * connectivity for the entire range of polys.  The second
   * is the last index of the connectivity data for each poly
   * in the global ID array.
   *
   * This function will add as many polys as possible to the
   * passed arrays given the sizes of those arrays.  It will
   * then pass back position at which it stoped and the sizes
   * of the data written to the arrays.
   *
   * Failure cases:
   *  - Passed range is empty (<code>begin == end</code>).
   *  - <code>element_array</code> or <code>index_array</code> is null.
   *  - The range contains invalid handles (non-existant entities,
   *      not an poly, etc.)
   *  - Retreiving ID tag for an entity failed.
   *
   *\param iter               As input, the first element handle.
   *                          As output, one past the last element handle
   *                          for which data was written to the arrays.
   *\param end                The iterator at which to stop.
   *\param node_id_tag        A tag with integer values.  
   *\param element_array_len  As input, length of <code>element_array</code>.
   *                          As output, the number of entries written in that
   *                          array.
   *\param element_array      The memory location at which to store the 
   *                          connectivity list.
   *\param index_array_len    As input, the length of <code>index_array</code>.
   *                          As output, the number of entries written in that
   *                          array.
   *\param index_array        The memory location at which to store offsets.
   *\param index_offset       Value to offset (add to) index values.  As output
   *                          the input value plus the amount of data 
   *                          written to the element array.  (The value you
   *                          presumably want to pass to the next call.)
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_poly_arrays(
      Range::const_iterator& iter,
      const Range::const_iterator end,
      const Tag node_id_tag,
      size_t& element_array_len,
      int *const element_array,
      size_t& index_array_len,
      int *const index_array,
      int& index_offset
      ) = 0;

    //! given elements to be written, gather all the nodes which define those elements
    //! \param elements Range of elements to be written
    //! \param node_bit_mark_tag Bit tag to use to identify nodes
    //! \param nodes Range of nodes gathered from elements (returned)
    //! \return status Return status
  virtual ErrorCode gather_nodes_from_elements(
    const Range& elements,
    const Tag node_bit_mark_tag,
    Range& nodes
    ) = 0;

    //! assign ids to input entities starting with start_id, written to id_tag
    //! if id_tag is zero, assigns to GLOBAL_ID_TAG_NAME
    //! \param elements Entities to be written
    //! \param id_tag Tag used to store entity id
    //! \param start_id Starting value for entity ids
    //! \return status Return status
  virtual ErrorCode assign_ids(Range &elements,
                                 Tag id_tag,
                                 const int start_id) = 0;
  

  /** Get explict adjacencies 
   *
   * Get explicit adjacences stored in database.
   * Does not create any explicit adjacencies or search for
   * implicit ones.
   *
   *\param entity  The entity to retreive adjacencies for.
   *\param id_tag  The global ID tag
   *\param adj     The output list of global IDs of adjacent entities.
   *\author Jason Kraftcheck
   */
  virtual ErrorCode get_adjacencies(
      EntityHandle entity,
      Tag id_tag,
      std::vector<int>& adj 
  ) = 0;
  
  
  virtual ErrorCode get_adjacencies( EntityHandle entity,
                               const EntityHandle*& adj_array,
                               int& num_adj ) = 0;

  
  /**\brief Re-order outgoing element connectivity
   *
   * Permute the connectivity of each element such that the node
   * order is that of the target file format rather than that of MBCN.
   *\param order The permutation to use.  Must be an array of 'node_per_elem'
   *             integers and be a permutation of the values [0..node_per_elem-1].
   *             Such that for a single element:
   *             target_conn[i] == mbcn_conn[order[i]]
   *\param conn  The connectivty array to re-order
   *\param num_elem  The number of elements in the connectivity array
   *\param node_per_elem The number of nodes in each element's connectivity list.
   */
  template <typename T> static inline 
  void reorder( const int* order, T* conn, int num_elem, int node_per_elem );


  /**\brief Get list of tags to write.
   *
   * Get the list of tags to write to the file, possibly using
   * an optional user-specifed tag list.
   *
   *\author Jason Kraftcheck
   */
  virtual ErrorCode 
  get_tag_list( std::vector<Tag>& result_list,
                const Tag* user_tag_list = 0, 
                int user_tag_list_length = 0,
                bool include_variable_length_tags = true ) = 0;

    //! if an error occured when reading the mesh, report it to MB
    //! it makes sense to have this as long as Interface has a write_mesh function
    //! \return status Return status
  virtual ErrorCode report_error( const std::string& error ) = 0;
  
    //! overloaded report_error behaves like the above
    //! \return status Return status
  virtual ErrorCode report_error( const char* error, ... ) 
#ifdef __GNUC__
__attribute__((format(printf,2,3)))
#endif
  = 0;

};

  
template <typename T> inline 
void WriteUtilIface::reorder( const int* order, T* conn, 
                                int num_elem, int node_per_elem )
{
  std::vector<T> elem(node_per_elem);
  T* const end = conn + num_elem*node_per_elem;
  while (conn != end) {
    std::copy( conn, conn+node_per_elem, elem.begin() );
    for (int j = 0; j < node_per_elem; ++j, ++conn) 
      *conn = elem[order[j]];
  }
}

} // namespace moab 

#endif 

