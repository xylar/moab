/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */



#ifndef MOAB_GEOM_TOPO_TOOL_HPP
#define MOAB_GEOM_TOPO_TOOL_HPP

#include "moab/Forward.hpp"
#include "moab/Range.hpp"


#include <map>

namespace moab {

// forward declare this class to avoid the header leaking in here
class OrientedBoxTreeTool;

/** \class GeomTopoTool
 * \brief Tool for interpreting geometric topology sets in MOAB database
 * Tool for interpreting geometric topology sets in MOAB database; see MOAB metadata_info
 * document for information on how geometric topology sets are read and represented.
 */
class GeomTopoTool
{
public:
  GeomTopoTool(Interface *impl, bool find_geoments = false, EntityHandle modelRootSet = 0);
  ~GeomTopoTool();
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
  ErrorCode restore_topology();

  //! Store sense of entity relative to wrt_entity.
     //!\return MB_MULTIPLE_ENTITIES_FOUND if surface already has a forward volume.
     //!        MB_SUCCESS if successful
     //!        otherwise whatever internal error code occured.
  
   ErrorCode set_sense( EntityHandle entity,
                        EntityHandle wrt_entity,
                        int sense);

     //! Get the sense of entity with respect to wrt_entity
     //! Returns MB_ENTITY_NOT_FOUND if no relationship found
   ErrorCode get_sense( EntityHandle entity,
                        EntityHandle wrt_entity,
                        int & sense );
     //! Get the senses of a surface with respect to its volumes
  ErrorCode get_surface_senses(EntityHandle surface_ent,
			       EntityHandle &forward_vol,
			       EntityHandle &reverse_vol);
  
     //! Set the senses of a surface with respect to its volumes  
  ErrorCode set_surface_senses(EntityHandle surface_ent,
			       EntityHandle forward_vol,
			       EntityHandle reverse_vol);

  ErrorCode get_senses (EntityHandle entity,
    std::vector<EntityHandle> &wrt_entities,
    std::vector<int> &senses);

  ErrorCode set_senses (EntityHandle entity,
                        std::vector<EntityHandle> &wrt_entities,
                        std::vector<int> &senses);
  
  /** Get the volume on the other side of a surface
   *
   * @param A surface to query
   * @param old_volume A volume on one side of surface
   * @param new_volume Output parameter for volume on the other side of surface
   * @return MB_SUCCESS if new_volume was set successfully, error if not.
   */
  ErrorCode next_vol( EntityHandle surface, EntityHandle old_volume,
                      EntityHandle& new_volume );

  
  // Retrieve geometry sets of desired dimension from model set
  //  0 = verts, 1 = curves, 2 = surfs, 3 = vols
  ErrorCode get_gsets_by_dimension( int dim, Range &gset);
    
  // Build obb tree for the entity set given; entity can be surface or volume
  ErrorCode construct_obb_tree(EntityHandle eh);

      // get the corners of the OBB for a given volume
  ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]);

    // get the center point and three vectors for the OBB of a given volume
  ErrorCode getobb(EntityHandle volume, double center[3],
                     double axis1[3], double axis2[3], double axis3[3]);

    /** \brief get the other (d-1)-dimensional entity bounding a set across a (d-2)-dimensional entity
     *
     * Given a d-dimensional entity and one (d-1)-dimensional entity, return the (d-1) dimensional
     * entity across a specified (d-2)-dimensional entity.  For example, given a surface, edge, and vertex,
     * returns the other edge bounding the surface sharing the vertex.  In the case of degenerate results,
     * e.g. two loops bounding a surface and sharing a vertex, tries to step in positively-oriented
     * direction.  This won't always work; in those cases, will return MB_MULTIPLE_ENTITIES_FOUND.
     *
     * In the special case where bounded is a curve, then not_this can be a vertex and across zero.
     * This function returns the other vertex on the curve.
     */
  ErrorCode other_entity(EntityHandle bounded, EntityHandle not_this, EntityHandle across,
                         EntityHandle &other);

    /** \brief return the dimension of the set, or -1 if it's not a geom_dimension set
     */
  int dimension(EntityHandle this_set);
  
  // used mostly for debugging purposes
  int global_id(EntityHandle this_set);

  // map from dimension & global ID to EntityHandle
  EntityHandle entity_by_id(int dimension, int id);

  ErrorCode find_geomsets(Range *ranges = NULL);

  // Build obb trees for all surfaces and volumes in model set
  // if make_one_vol true, joins trees from all surfaces in model into single
  // volume obb tree
  ErrorCode construct_obb_trees(bool make_one_vol = false);
  
  /* Relies on future work in OBBTreeTool before final implementation */
  //ErrorCode delete_obb_tree(EntityHandle eh);
  
  ErrorCode remove_root(EntityHandle vol_or_surf);
  
  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle &root);

  EntityHandle get_one_vol_root();

  OrientedBoxTreeTool *obb_tree() {return obbTree;}

  // this could make the obb tree out of date
  ErrorCode add_geo_set(EntityHandle set, int dimension, int global_id  = 0);

  // will assume no geo sets are defined for this surface
  // will output a mesh_set that contains everything (all sets of interest), for proper output
  ErrorCode geometrize_surface_set(EntityHandle surface, EntityHandle & output);

  // get the implicit complement handle
  ErrorCode get_implicit_complement(EntityHandle &implicit_complement, bool create_if_missing = false);

  ErrorCode is_owned_set(EntityHandle eh);
  
  // this would be a deep copy, into a new geom topo tool
  // sets will be duplicated, but entities not
  // modelSet will be a new one;
  // will take as input a pointer to a std::vector of gents (surfaces and volumes, usually),
  // which will serve to filter the gents from modelSet (only dependents will be part of the new gtt)
  // if the pointer is null, all gsets in the original modelSet are duplicated

  ErrorCode duplicate_model(GeomTopoTool *& duplicate, std::vector<EntityHandle> * pvGEnts = NULL);

  EntityHandle get_root_model_set() { return modelSet; }

  bool check_model();
  // should be used instead of keeping multiple ranges, for example in FBEngine
  const Range * geoRanges() { return geomRanges ; }

  Interface* get_moab_instance() { return mdbImpl; }

  Tag get_sense_tag();

  Tag get_gid_tag();
  
  Tag get_geom_tag();

  bool have_obb_tree();
  
private:
  Interface *mdbImpl;
  Tag sense2Tag;
  Tag senseNEntsTag, senseNSensesTag;
  Tag geomTag;
  Tag gidTag;
  Tag nameTag;
  // the model set encompasses a full topological model
  EntityHandle modelSet;
  Range geomRanges[5];// add one more dimension, for set of gentities; by default, they will
                      // have geom_dimension 4
  int maxGlobalId[5]; // one max global id for each dimension
  bool updated;

  OrientedBoxTreeTool* obbTree;
  EntityHandle setOffset;
  std::vector<EntityHandle> rootSets;

  bool contiguous;
  std::map<EntityHandle, EntityHandle>  mapRootSets;
  EntityHandle oneVolRootSet;

    //! Creates a volume for undefined space in the model
    // The implicit complement is composed of all surfaces that only
    // have one parent volume, i.e. surfaces that are in contact with the outside
    // world
  ErrorCode create_implicit_complement(EntityHandle &implicit_complement);
  
    //! compute vertices inclusive and put on tag on sets in geom_sets
  ErrorCode construct_vertex_ranges(const Range &geom_sets,
				      const Tag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  ErrorCode separate_by_dimension(const Range &geom_sets);

    //! verify global id tag
  ErrorCode check_gid_tag(bool create = false);

    //! verify geometry tag
  ErrorCode check_geom_tag(bool create = false);

    //! verify sense face tag
  ErrorCode check_face_sense_tag(bool create = false);

    //! verify sense edge tags
  ErrorCode check_edge_sense_tags(bool create = false);

    //! Set the contigous variable
    //  If it has changed, update the storage of the rootsets 
  void set_contiguous(bool new_value);  

    //! Test if the entity sets are contiguous or not
  ErrorCode update_contiguous();

};

// get the root of the obbtree for a given entity
inline ErrorCode GeomTopoTool::get_root(EntityHandle vol_or_surf, EntityHandle &root) 
{
   if(contiguous)
   {
     unsigned int index = vol_or_surf - setOffset;
     root = (index < rootSets.size() ? rootSets[index] : 0);
   }
   else
      root = mapRootSets[vol_or_surf];
   return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

inline ErrorCode GeomTopoTool::remove_root(EntityHandle vol_or_surf) {
   if(contiguous)
   {
     unsigned int index = vol_or_surf - setOffset;
     if( index < rootSets.size() ) {
       rootSets[index] = 0;
     }
     else {
       return MB_INDEX_OUT_OF_RANGE;
     }
   }
   else {
      mapRootSets[vol_or_surf] = 0;
   }
   
   return MB_SUCCESS;
}
  
inline EntityHandle GeomTopoTool::get_one_vol_root()
{
  return oneVolRootSet;
}

  inline Tag GeomTopoTool::get_sense_tag() { check_face_sense_tag(true); return sense2Tag; }
  
  inline Tag GeomTopoTool::get_gid_tag() { check_gid_tag(true); return gidTag; }
  
  inline Tag GeomTopoTool::get_geom_tag() { check_geom_tag(true); return geomTag; }

} // namespace moab 

#endif

