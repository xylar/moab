#include "iGeom_MOAB.hpp"
#include "iMesh.h"
#include "MBInterface.hpp"
#include "GeomTopoTool.hpp"
#include "MBOrientedBox.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBCartVect.hpp"

iBase_Error iGeom_LAST_ERROR;
bool i_created = false; // if interface is created
bool t_created = false;
MBRange _my_gsets[4];
GeomTopoTool* _my_geomTopoTool = NULL;

//int treeOffset = 0; // if obb tree is created, set to the first root set
//MBRange _my_treeRootSets;
//MBOrientedBoxTreeTool _my_obbTree;

#define COPY_RANGE(r, vec) {                      \
    MBEntityHandle *tmp_ptr = reinterpret_cast<MBEntityHandle*>(vec);	\
    std::copy(r.begin(), r.end(), tmp_ptr);}
 
static inline void
iGeom_processError(iBase_ErrorType code, const char* desc);

static void
iGeom_get_adjacent_entities(iGeom_Instance instance,
			    const MBEntityHandle from, 
			    const int to_dim,
			    MBRange &adj_ents, int* err);

void iGeom_getDescription( iGeom_Instance instance,
                           char* descr,
                           int* err,
                           int descr_len )
{
  unsigned int len = MIN(strlen(iGeom_LAST_ERROR.description), ((unsigned int) descr_len));
  strncpy(descr, iGeom_LAST_ERROR.description, len);
  descr[len] = '\0';
  RETURN(iBase_SUCCESS);
}

void iGeom_getErrorType(iGeom_Instance instance,
                        /*out*/ int *error_type, 
                        int *err)
{
  *error_type = iGeom_LAST_ERROR.error_type;
  RETURN(iBase_SUCCESS);
}

void iGeom_newGeom( char const* options,
                    iGeom_Instance* instance_out,
                    int* err,
                    int options_len ) 
{
  if (*instance_out && !(reinterpret_cast<MBInterface*>(instance_out))) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Passed in instance must be an MBInterface*.");
  }

  // make a new imesh instance
  iMesh_newMesh(options, reinterpret_cast<iMesh_Instance*>(instance_out),
		err, options_len);
  ERRORR("Failure to create instance.");

  i_created = true;

  RETURN(iBase_SUCCESS);
}
  
void iGeom_dtor( iGeom_Instance instance, int* err ) 
{
  if (i_created) {
    iMesh_dtor(IMESH_INSTANCE(instance), err);
    ERRORR("Failed to destruct instance.");
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_load( iGeom_Instance instance,
                 char const* name,
                 char const* options,
                 int* err,
                 int name_len,
                 int options_len ) 
{
  // load mesh-based geometry
  iMesh_load(IMESH_INSTANCE(instance), NULL, name, options, err, 
             name_len, options_len);
  ERRORR("Failure to load geometry.");
  
  // keep mesh-based geometries in MBRange
  GETGTT(instance);
  MBErrorCode rval = _my_geomTopoTool->find_geomsets(_my_gsets);
  MBERRORR("Failure to keep geometry list.");

  RETURN(iBase_SUCCESS);
}

  
void iGeom_save( iGeom_Instance instance,
                 char const* name,
                 char const* options,
                 int* err,
                 int name_len,
                 int options_len ) 
{
  iMesh_save(IMESH_INSTANCE(instance), NULL, name, options, err, 
             name_len, options_len);
  ERRORR("Failed get save geomtry instance.");

  RETURN(iBase_SUCCESS);
}

void iGeom_getRootSet( iGeom_Instance instance,
                       iBase_EntitySetHandle* root_set,
                       int* err )
{
  iMesh_getRootSet(IMESH_INSTANCE(instance), root_set, err);
  ERRORR("Failed get root set.");

  RETURN(iBase_SUCCESS);
}

void iGeom_getBoundBox( iGeom_Instance,
                        double* min_x,
                        double* min_y,
                        double* min_z,
                        double* max_x,
                        double* max_y,
                        double* max_z,
                        int* err )
{
  RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntities( iGeom_Instance instance,
                        iBase_EntitySetHandle set_handle,
                        int entity_type,
                        iBase_EntityHandle** entity_handles,
                        int* entity_handles_allocated,
                        int* entity_handles_size,
                        int* err ) 
{
  int i;
  if (0 > entity_type || 4 < entity_type) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Bad entity type.");
  }
  else if (entity_type < 4) {
    *entity_handles_size = _my_gsets[entity_type].size();
    CHECK_SIZE(*entity_handles, *entity_handles_allocated,
	       *entity_handles_size, iBase_EntityHandle, NULL);
    COPY_RANGE(_my_gsets[entity_type], *entity_handles);
  }
  else {
    *entity_handles_size = 0;
    MBRange total_range;
    for (i = 0; i < 4; i++) {
      total_range.merge(_my_gsets[i]);
    }
    *entity_handles_size = total_range.size();
    CHECK_SIZE(*entity_handles, *entity_handles_allocated,
	       *entity_handles_size, iBase_EntityHandle, NULL);
    COPY_RANGE(total_range, *entity_handles);
  }
  
  RETURN(iBase_SUCCESS);
}
  
void iGeom_getNumOfType( iGeom_Instance instance,
                         iBase_EntitySetHandle set_handle,
                         int entity_type,
                         int* num_out,
                         int* err ) 
{
  if (0 > entity_type || 3 < entity_type) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Bad entity type.");
  }
  *num_out = _my_gsets[entity_type].size();
  
  RETURN(iBase_SUCCESS);
}

void iGeom_getEntType( iGeom_Instance instance,
                       iBase_EntityHandle entity_handle,
                       int* type,
                       int* err ) 
{
  for (int i = 0; i < 4; i++) {
    if (_my_gsets[i].find(MBH_cast(entity_handle)) != _my_gsets[i].end()) {
      *type = i;
      RETURN(iBase_SUCCESS);
    }
  }
  
  *err = iBase_INVALID_ENTITY_TYPE;
  ERRORR("Entity not a geometry entity.");
}

void iGeom_getArrType( iGeom_Instance instance,
                       iBase_EntityHandle const* entity_handles,
                       int entity_handles_size,
                       int** type,
                       int* type_allocated,
                       int* type_size,
                       int* err ) 
{
  CHECK_SIZE(*type, *type_allocated, *type_size, int, NULL);

  int tmp_err;
  *err = iBase_SUCCESS;
  
  for (int i = 0; i < entity_handles_size; i++) {
    iGeom_getEntType(instance, entity_handles[i], *type + i, &tmp_err);
    if (iBase_SUCCESS != tmp_err) {
      *err = tmp_err;
      ERRORR("Failed to get entity type in iGeom_getArrType.");
    }
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getEntAdj( iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      int to_dimension,
                      iBase_EntityHandle** adj_entities,
                      int* adj_entities_allocated,
                      int* adj_entities_size,
                      int* err ) 
{
  MBRange adjs;
  MBEntityHandle this_ent = MBH_cast(entity_handle);

  // get adjacent
  iGeom_get_adjacent_entities(instance, this_ent, to_dimension,
			      adjs, err);
  ERRORR("Failed to get adjacent entities in iGeom_getEntAdj.");
  
  // copy adjacent entities
  CHECK_SIZE(*adj_entities, *adj_entities_allocated,
	     (int) adjs.size(), iBase_EntityHandle, NULL);
  COPY_RANGE(adjs, *adj_entities);

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrAdj( iGeom_Instance instance,
                      iBase_EntityHandle const* entity_handles,
                      int entity_handles_size,
                      int requested_entity_type,
                      iBase_EntityHandle** adj_entity_handles,
                      int* adj_entity_handles_allocated,
                      int* adj_entity_handles_size,
                      int** offset,
                      int* offset_allocated,
                      int* offset_size,
                      int* err ) 
{
  // check offset array size
  MBRange temp_range, total_range;
  CHECK_SIZE(*offset, *offset_allocated, entity_handles_size + 1, int, NULL);
  *offset_size = entity_handles_size + 1;
  
  // get adjacent entities
  for (int i = 0; i < entity_handles_size; ++i) {
    (*offset)[i] = total_range.size();
    temp_range.clear();
    iGeom_get_adjacent_entities(instance, MBH_cast(entity_handles[i]),
				requested_entity_type, temp_range, err);
    ERRORR("Failed to get adjacent entities in iGeom_getEntAdj.");

    total_range.merge(temp_range);
  }
  int nTot = total_range.size();
  (*offset)[entity_handles_size] = nTot;

  // copy adjacent entities
  CHECK_SIZE(*adj_entity_handles, *adj_entity_handles_allocated,
	     nTot, iBase_EntityHandle, NULL);
  COPY_RANGE(total_range, *adj_entity_handles);
  *adj_entity_handles_size = nTot;

  RETURN(iBase_SUCCESS);
}
  
void iGeom_getEnt2ndAdj( iGeom_Instance instance,
                         iBase_EntityHandle entity_handle,
                         int bridge_dimension,
                         int to_dimension,
                         iBase_EntityHandle** adjacent_entities,
                         int* adjacent_entities_allocated,
                         int* adjacent_entities_size,
                         int* err )
{
  MBRange to_ents, bridge_ents, tmp_ents;
  iGeom_get_adjacent_entities(instance, MBH_cast(entity_handle), bridge_dimension, bridge_ents, err);
  ERRORR("Failed to get adjacent entities in iGeom_getEnt2ndAdj.");
  
  MBRange::iterator iter, jter, kter, end_jter;
  MBRange::iterator end_iter = bridge_ents.end();
  for (iter = bridge_ents.begin(); iter != end_iter; iter++) {
    iGeom_get_adjacent_entities(instance, *iter, to_dimension, tmp_ents, err);
    ERRORR("Failed to get adjacent entities in iGeom_getEnt2ndAdj.");
    
    for (jter = tmp_ents.begin(); jter != end_jter; jter++) {
      if (to_ents.find(*jter) == to_ents.end()) {
	to_ents.insert(*jter);
      }
    }
    tmp_ents.clear();
  }

  *adjacent_entities_size = to_ents.size();
  CHECK_SIZE(*adjacent_entities, *adjacent_entities_allocated,
	     *adjacent_entities_size, iBase_EntityHandle, NULL);
  COPY_RANGE(to_ents, *adjacent_entities);

  RETURN(iBase_SUCCESS);
}

void iGeom_getArr2ndAdj( iGeom_Instance instance,
                         iBase_EntityHandle const* entity_handles,
                         int entity_handles_size,
                         int order_adjacent_key,
                         int requested_entity_type,
                         iBase_EntityHandle** adj_entity_handles,
                         int* adj_entity_handles_allocated,
                         int* adj_entity_handles_size,
                         int** offset,
                         int* offset_allocated,
                         int* offset_size,
                         int* err )
{
  CHECK_SIZE(*offset, *offset_allocated, entity_handles_size + 1, int, NULL);
  MBRange bridge_range, temp_range, entity_range, total_range;
   
  for (int i = 0; i < entity_handles_size; ++i) {
    bridge_range.clear();
    entity_range.clear();
    iGeom_get_adjacent_entities(instance, MBH_cast(entity_handles[i]),
				order_adjacent_key, bridge_range, err);
    ERRORR("Failed to get adjacent entities in iGeom_getArr2ndAdj.");

    MBRange::iterator iter, jter, end_jter;
    MBRange::iterator end_iter = bridge_range.end();
    for (iter = bridge_range.begin(); iter != end_iter; iter++) {
      temp_range.clear();
      iGeom_get_adjacent_entities(instance, *iter,
				  requested_entity_type, temp_range, err);
      ERRORR("Failed to get adjacent entities in iGeom_getArr2ndAdj.");

      for (jter = temp_range.begin(); jter != end_jter; jter++) {
	if (entity_range.find(*jter) == entity_range.end()) {
	  entity_range.insert(*jter);
	}
      }
    }

    (*offset)[i] = total_range.size();
    total_range.merge(entity_range);
  }
  *adj_entity_handles_size = total_range.size();
  (*offset)[entity_handles_size] = *adj_entity_handles_size;

  CHECK_SIZE(*adj_entity_handles, *adj_entity_handles_allocated,
	     *adj_entity_handles_size, iBase_EntityHandle, NULL);
  COPY_RANGE(total_range, *adj_entity_handles);

  RETURN(iBase_SUCCESS);
}

void iGeom_isEntAdj( iGeom_Instance instance,
                     iBase_EntityHandle entity_handle1,
                     iBase_EntityHandle entity_handle2,
                     int* are_adjacent,
                     int* err )
{
  int type1, type2;
  iGeom_getEntType(instance, entity_handle1, &type1, err);
  ERRORR("Failed to get entity type in iGeom_isEntAdj.");
  iGeom_getEntType(instance, entity_handle2, &type2, err);
  ERRORR("Failed to get entity type in iGeom_isEntAdj.");

  MBErrorCode rval;
  MBRange adjs;
  if (type1 < type2) {
    rval = MBI->get_parent_meshsets(MBH_cast(entity_handle1), adjs, type2 - type1);
    MBERRORR("Failed to get parent meshsets in iGeom_isEntAdj.");
  }
  else {
    rval = MBI->get_child_meshsets(MBH_cast(entity_handle1), adjs, type2 - type1);
    MBERRORR("Failed to get child meshsets in iGeom_isEntAdj.");
  }
  
  *are_adjacent = adjs.find(MBH_cast(entity_handle2)) != _my_gsets[type2].end();

  RETURN(iBase_SUCCESS);
}

void iGeom_isArrAdj( iGeom_Instance instance,
                     iBase_EntityHandle const* entity_handles_1,
                     int entity_handles_1_size,
                     iBase_EntityHandle const* entity_handles_2,
                     int entity_handles_2_size,
                     int** is_adjacent_info,
                     int* is_adjacent_info_allocated,
                     int* is_adjacent_info_size,
                     int* err )
{
  int index1 = 0;
  int index2 = 0;
  size_t index1_step, index2_step;
  int count;
    
  // If either list contains only 1 entry, compare that entry with
  // every entry in the other list.
  if (entity_handles_1_size == entity_handles_2_size) {
    index1_step = index2_step = 1;
    count = entity_handles_1_size;
  }
  else if (entity_handles_1_size == 1) {
    index1_step = 0;
    index2_step = 1;
    count = entity_handles_2_size;
  }
  else if (entity_handles_2_size == 1) {
    index1_step = 1;
    index2_step = 0;
    count = entity_handles_1_size;
  }
  else {
    RETURN(iBase_INVALID_ENTITY_COUNT);
  }
  
  CHECK_SIZE(*is_adjacent_info, *is_adjacent_info_allocated,
	     count, int, NULL);

  for (int i = 0; i < count; ++i)
  {
    iGeom_isEntAdj(instance, entity_handles_1[index1],
		   entity_handles_2[index2],
		   &((*is_adjacent_info)[i]), err);
    ERRORR("Failed to check if entities are adjacent in iGeom_isArrAdj.");

    index1 += index1_step;
    index2 += index2_step;
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getEntClosestPt( iGeom_Instance instance,
                            iBase_EntityHandle entity_handle,
                            double near_x, 
                            double near_y,
                            double near_z,
                            double* on_x,
                            double* on_y,
                            double* on_z,
                            int* err )
{
  MBErrorCode rval;
  int type;
  iGeom_getEntType(instance, entity_handle, &type, err);
  ERRORR("Failed to get entity type.");

  if (type == 0) {
    iGeom_getVtxCoord(instance, entity_handle,
		      on_x, on_y, on_z, err);
    ERRORR("Failed to get vertex coordinates.");
  }
  else if (type == 1) {
    // just copy over the coordinates
    // should be modified
    *on_x = near_x;
    *on_y = near_y;
    *on_z = near_z;
  }
  else if (type == 2 || type == 3) {
    if (!t_created) {
      GETGTT(instance);
      rval = _my_geomTopoTool->construct_obb_trees();
      MBERRORR("Failed to construct obb tree.");
      t_created = true;
    }
    
    double point[3] = {near_x, near_y, near_z};
    double point_out[3];
    MBEntityHandle root, facet_out;
    _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
    MBERRORR("Failed to get tree root in iGeom_getEntClosestPt.");
    rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
							     point_out,
							     facet_out);
    MBERRORR("Failed to get closest point in iGeom_getEntClosestPt.");
    
    *on_x = point_out[0];
    *on_y = point_out[1];
    *on_z = point_out[2];
  }
  else RETURN(iBase_INVALID_ENTITY_TYPE);

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrClosestPt( iGeom_Instance instance,
                            iBase_EntityHandle const* entity_handles,
                            int entity_handles_size,
                            int storage_order,
                            double const* near_coordinates,
                            int near_coordinates_size,
                            double** on_coordinates,
                            int* on_coordinates_allocated,
                            int* on_coordinates_size,
                            int* err )
{
  CHECK_SIZE(*on_coordinates, *on_coordinates_allocated,
	     near_coordinates_size, double, NULL);
  for (int i = 0; i < entity_handles_size; i++) {
    if (storage_order == iBase_INTERLEAVED) {
      iGeom_getEntClosestPt(instance, entity_handles[i], near_coordinates[3*i],
			    near_coordinates[3*i + 1], near_coordinates[3*i + 2],
			    on_coordinates[3*i], on_coordinates[3*i + 1],
			    on_coordinates[3*i + 2], err);
    }
    else if (storage_order == iBase_BLOCKED) {
      iGeom_getEntClosestPt(instance, entity_handles[i], near_coordinates[i],
			    near_coordinates[i + entity_handles_size],
			    near_coordinates[i + 2*entity_handles_size],
			    on_coordinates[i], on_coordinates[i + entity_handles_size],
			    on_coordinates[i + 2*entity_handles_size], err);
    }
    ERRORR("Failed to get closest point.");
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getEntNrmlXYZ( iGeom_Instance instance,
                          iBase_EntityHandle entity_handle,
                          double x,
                          double y,
                          double z,
                          double* nrml_i,
                          double* nrml_j,
                          double* nrml_k,
                          int* err )
{
  // just do for surface and volume
  int type;
  iGeom_getEntType(instance, entity_handle, &type, err);
  ERRORR("Failed to get entity type in iGeom_getEntNrmlXYZ.");

  if (type != 2 && type != 3) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Entities passed into gentityNormal must be face or volume.");
  }

  // get closest location and facet
  double point[3] = {x, y, z};
  double point_out[3];
  MBEntityHandle root, facet_out;
  _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
  MBErrorCode rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
							   point_out,
							   facet_out);
  MBERRORR("Failed to get closest location in iGeom_getEntNrmlXYZ.");

  // get facet normal
  const MBEntityHandle* conn;
  int len, sense;
  MBCartVect coords[3], normal;
  rval = MBI->get_connectivity(facet_out, conn, len);
  MBERRORR("Failed to get triangle connectivity in iGeom_getEntNrmlXYZ.");
  if (len != 3) RETURN(iBase_FAILURE);
  
  rval = MBI->get_coords(conn, len, coords[0].array());
  MBERRORR("Failed to get triangle coordinates in iGeom_getEntNrmlXYZ.");
  
  coords[1] -= coords[0];
  coords[2] -= coords[0];
  normal = coords[1] * coords[2];
  normal.normalize();
  *nrml_i = normal[0];
  *nrml_j = normal[1];
  *nrml_k = normal[2];

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrNrmlXYZ( iGeom_Instance instance,
                          iBase_EntityHandle const* entity_handles,
                          int entity_handles_size,
                          int storage_order,
                          double const* coordinates,
                          int coordinates_size,
                          double** normals,
                          int* normals_allocated,
                          int* normals_size,
                          int* err )
{
  // set up iteration according to storage order.
  // allow either gentity_handles or near_coordinates to contain
  // only one value, where that single value is applied for every
  // entry in the other list.
  size_t index = 0;
  size_t coord_step, norm_step = 1, ent_step;
  int count;
  if (3*entity_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = entity_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = entity_handles_size;
  }
  else if (entity_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    *err = iBase_INVALID_ENTITY_COUNT;
    ERRORR("Mismatched array sizes");
  }

  // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*normals, *normals_allocated, 3*count, double, NULL);
  
  const double *coord_x, *coord_y, *coord_z;
  double *norm_x, *norm_y, *norm_z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    norm_x = *normals;
    norm_y = norm_x + count;
    norm_z = norm_y + count;
    norm_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    norm_x = *normals;
    norm_y = norm_x+1;
    norm_z = norm_x+2;
    coord_step *= 3;
    norm_step = 3;
  }
  
  for (int i = 0; i < count; ++i)
  {
    iGeom_getEntNrmlXYZ(instance, entity_handles[index],
			*coord_x, *coord_y, *coord_z,
			norm_x, norm_y, norm_z, err);
    ERRORR("Failed to get entity normal of point.");

    //entities += ent_step;
    index += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    norm_x += norm_step;
    norm_y += norm_step;
    norm_z += norm_step;
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getEntNrmlPlXYZ( iGeom_Instance instance,
                            iBase_EntityHandle entity_handle,
                            double x,
                            double y,
                            double z,
                            double* pt_x,
                            double* pt_y,
                            double* pt_z,
                            double* nrml_i,
                            double* nrml_j,
                            double* nrml_k,
                            int* err )
{
  // just do for surface and volume
  int type;
  iGeom_getEntType(instance, entity_handle, &type, err);
  ERRORR("Failed to get entity type in iGeom_getEntNrmlPlXYZ.");

  if (type != 2 && type != 3) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Entities passed into gentityNormal must be face or volume.");
  }

  // get closest location and facet
  double point[3] = {x, y, z};
  double point_out[3];
  MBEntityHandle root, facet_out;
  _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
  MBErrorCode rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
								       point_out,
								       facet_out);
  MBERRORR("Failed to get closest location in iGeom_getEntNrmlPlXYZ.");

  // get closest point
  *pt_x = point_out[0];
  *pt_y = point_out[1];
  *pt_z = point_out[2];

  // get facet normal
  const MBEntityHandle* conn;
  int len, sense;
  MBCartVect coords[3], normal;
  rval = MBI->get_connectivity(facet_out, conn, len);
  MBERRORR("Failed to get triangle connectivity in iGeom_getEntNrmlPlXYZ.");
  if (len != 3) RETURN(iBase_FAILURE);
  
  rval = MBI->get_coords(conn, len, coords[0].array());
  MBERRORR("Failed to get triangle coordinates in iGeom_getEntNrmlPlXYZ.");
  
  coords[1] -= coords[0];
  coords[2] -= coords[0];
  normal = coords[1] * coords[2];
  normal.normalize();
  *nrml_i = normal[0];
  *nrml_j = normal[1];
  *nrml_k = normal[2];

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrNrmlPlXYZ( iGeom_Instance instance,
                            iBase_EntityHandle const* entity_handles,
                            int entity_handles_size,
                            int storage_order,
                            double const* near_coordinates,
                            int near_coordinates_size,
                            double** on_coordinates,
                            int* on_coordinates_allocated,
                            int* on_coordinates_size,
                            double** normals,
                            int* normals_allocated,
                            int* normals_size,
                            int* err )
{
  // set up iteration according to storage order.
  // allow either gentity_handles or near_coordinates to contain
  // only one value, where that single value is applied for every
  // entry in the other list.
  size_t index = 0;
  size_t near_step, on_step = 1, ent_step;
  int count;
  if (3*entity_handles_size == near_coordinates_size) {
    near_step = ent_step = 1;
    count = entity_handles_size;
  }
  else if (near_coordinates_size == 3) {
    near_step = 0;
    ent_step = 1;
    count = entity_handles_size;
  }
  else if (entity_handles_size == 1) {
    near_step = 1;
    ent_step = 0;
    count = near_coordinates_size / 3;
  }
  else {
    *err = iBase_INVALID_ENTITY_COUNT;
    ERRORR("Mismatched array sizes");
  }

  // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*on_coordinates, *on_coordinates_allocated, 3*count, double, NULL);
  CHECK_SIZE(*normals, *normals_allocated, 3*count, double, NULL);
  
  const double *near_x, *near_y, *near_z;
  double *on_x, *on_y, *on_z;
  double *norm_x, *norm_y, *norm_z;
  if (storage_order == iBase_BLOCKED) {
    near_x = near_coordinates;
    near_y = near_x + near_coordinates_size/3;
    near_z = near_y + near_coordinates_size/3;
    on_x = *on_coordinates;
    on_y = on_x + count;
    on_z = on_y + count;
    norm_x = *normals;
    norm_y = norm_x + count;
    norm_z = norm_y + count;
    on_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    near_x = near_coordinates;
    near_y = near_x+1;
    near_z = near_x+2;
    on_x = *on_coordinates;
    on_y = on_x+1;
    on_z = on_x+2;
    norm_x = *normals;
    norm_y = norm_x+1;
    norm_z = norm_x+2;
    near_step *= 3;
    on_step = 3;
  }
  
  for (int i = 0; i < count; ++i)
  {
    iGeom_getEntNrmlPlXYZ(instance, entity_handles[index],
			  *near_x, *near_y, *near_z,
			  on_x, on_y, on_z,
			  norm_x, norm_y, norm_z, err);
    ERRORR("Failed to get entity normal of point.");

    //entities += ent_step;
    index += ent_step;
    near_x += near_step;
    near_y += near_step;
    near_z += near_step;
    on_x += on_step;
    on_y += on_step;
    on_z += on_step;
    norm_x += on_step;
    norm_y += on_step;
    norm_z += on_step;
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getEntTgntXYZ( iGeom_Instance,
                          iBase_EntityHandle entity_handle,
                          double x,
                          double y,
                          double z,
                          double* tgnt_i,
                          double* tgnt_j,
                          double* tgnt_k,
                          int* err )
{
  RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getArrTgntXYZ( iGeom_Instance,
                          iBase_EntityHandle const* entity_handles,
                          int entity_handles_size,
                          int storage_order,
                          double const* coordinates,
                          int coordinates_size,
                          double** tangents,
                          int* tangents_allocated,
                          int* tangents_size,
                          int* err )
{
  RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntBoundBox( iGeom_Instance instance,
                           iBase_EntityHandle entity_handle,
                           double* min_x,
                           double* min_y,
                           double* min_z,
                           double* max_x,
                           double* max_y,
                           double* max_z,
                           int* err )
{
  MBErrorCode rval;
  int type;
  iGeom_getEntType(instance, entity_handle, &type, err);
  ERRORR("Failed to get entity type.");

  if (type == 0) {
    iGeom_getVtxCoord(instance, entity_handle,
		      min_x, min_y, min_z, err);
    ERRORR("Failed to get vertex coordinates.");
    max_x = min_x;
    max_y = min_y;
    max_z = min_z;
  }
  else if (type == 1) {
    *err = iBase_NOT_SUPPORTED;
    ERRORR("iGeom_getEntBoundBox is not supported for Edge entity type.");
  }
  else if (type == 2 || type == 3) {
    if (!t_created) {
      GETGTT(instance);
      rval = _my_geomTopoTool->construct_obb_trees();
      MBERRORR("Failed to construct obb tree.");
      t_created = true;
    }
    
    MBEntityHandle root;
    MBOrientedBox box;
    rval = _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
    MBERRORR("Failed to get tree root in iGeom_getEntBoundBox.");
    rval = _my_geomTopoTool->obb_tree()->box(root, box);
    MBERRORR("Failed to get closest point in iGeom_getEntBoundBox.");
    
    MBCartVect min, max;
    min = box.center - box.length[0]*box.axis[0] - box.length[1]*box.axis[1]
      - box.length[2]*box.axis[2];
    max = box.center + box.length[0]*box.axis[0] + box.length[1]*box.axis[1]
      + box.length[2]*box.axis[2];
    *min_x = min[0];
    *min_y = min[1];
    *min_z = min[2];
    *max_x = max[0];
    *max_y = max[1];
    *max_z = max[2];
  }
  else RETURN(iBase_INVALID_ENTITY_TYPE);

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrBoundBox( iGeom_Instance instance,
                           iBase_EntityHandle const* entity_handles,
                           int entity_handles_size,
                           int storage_order,
                           double** min_corner,
                           int* min_corner_allocated,
                           int* min_corner_size,
                           double** max_corner,
                           int* max_corner_allocated,
                           int* max_corner_size,
                           int* err )
{
  // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*min_corner, *min_corner_allocated, 3*entity_handles_size, double, NULL);
  CHECK_SIZE(*max_corner, *max_corner_allocated, 3*entity_handles_size, double, NULL);
  
  size_t step, init;
  if (storage_order == iBase_BLOCKED) {
    step = 1;
    init = entity_handles_size;
  }
  else {
    step = 3;
    init = 1;
  }
  double *min_x, *min_y, *min_z, *max_x, *max_y, *max_z;
  min_x = *min_corner;
  max_x = *max_corner;
  min_y = min_x + init;
  max_y = max_x + init;
  min_z = min_y + init;
  max_z = max_y + init;
  
  for (int i = 0; i < entity_handles_size; ++i)
  {
    iGeom_getEntBoundBox(instance, entity_handles[i],
			 min_x, min_y, min_z,
			 max_x, max_y, max_z, err);
    ERRORR("Failed to get entity bounding box in iGeom_getArrBoundBox.");
    
    min_x += step;
    max_x += step;
    min_y += step;
    max_y += step;
    min_z += step;
    max_z += step;
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getVtxCoord( iGeom_Instance instance,
                        iBase_EntityHandle vertex_handle,
                        double* x,
                        double* y,
                        double* z,
                        int* err )
{
  int type;
  iGeom_getEntType(instance, vertex_handle, &type, err);
  ERRORR("Failed to get entity type in iGeom_getVtxCoord.");

  if (type != 0) {
    *err = iBase_INVALID_ENTITY_TYPE;
    ERRORR("Entity is not a vertex type.");
  }

  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0, verts_size;
  double *vert_coords = NULL;
  int vert_coords_alloc = 0, vert_coords_size;
  iMesh_getEntities(IMESH_INSTANCE(instance), reinterpret_cast<iBase_EntitySetHandle> (vertex_handle),
		    iBase_VERTEX, iMesh_POINT, &verts, &verts_alloc, &verts_size, err);
  ERRORR("Failed to get vertices.");

  if (verts_size != 1) {
    *err = iBase_FAILURE;
    ERRORR("Vertex has multiple points.");
  }

  iMesh_getVtxCoord(IMESH_INSTANCE(instance), verts[0],
		    x, y, z, err);
  ERRORR("Failed to get vertex coordinate.");

  RETURN(iBase_SUCCESS);
}

void iGeom_getVtxArrCoords( iGeom_Instance instance,
                            iBase_EntityHandle const* entity_handles,
                            int entity_handles_size,
                            int storage_order,
                            double** coordinates,
                            int* coordinates_allocated,
                            int* coordinates_size,
                            int* err )
{
  // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*coordinates, *coordinates_allocated, 3*entity_handles_size, double, NULL);

  double *x, *y, *z;
  size_t step;
  if (storage_order == iBase_BLOCKED) {
    x = *coordinates;
    y = x + entity_handles_size;
    z = y + entity_handles_size;
    step = 1;
  }
  else {
    x = *coordinates;
    y = x + 1;
    z = x + 2;
    step = 3;
  }
  
  for (int i = 0; i < entity_handles_size; i++) {
    iGeom_getVtxCoord(instance, entity_handles[i], x, y, z, err);
    x += step;
    y += step;
    z += step;
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_getPntRayIntsct( iGeom_Instance,
                            double x,
                            double y,
                            double z,
                            double dir_x,
                            double dir_y,
                            double dir_z,
                            iBase_EntityHandle** intersect_entity_handles,
                            int* intersect_entity_handles_allocated,
                            int* intersect_entity_hangles_size,
                            int storage_order,
                            double** intersect_coords,
                            int* intersect_coords_allocated,
                            int* intersect_coords_size,
                            double** param_coords,
                            int* param_coords_allocated,
                            int* param_coords_size,
                            int* err )
{
 }

void iGeom_getPntArrRayIntsct( iGeom_Instance,
                               int storage_order,
                               const double* coords,
                               int coords_size,
                               const double* directions,
                               int directions_size,
                               iBase_EntityHandle** intersect_entity_handles,
                               int* intersect_entity_handles_allocated,
                               int* intersect_entity_hangles_size,
                               int** offset,
                               int* offset_allocated,
                               int* offset_size,
                               double** intersect_coords,
                               int* intersect_coords_allocated,
                               int* intersect_coords_size,
                               double** param_coords,
                               int* param_coords_allocated,
                               int* param_coords_size,
                               int* err ){ }
void iGeom_getEntNrmlSense( iGeom_Instance,
                            iBase_EntityHandle face,
                            iBase_EntityHandle region,
                            int* sense_out,
                            int* err ){ }
void iGeom_getArrNrmlSense( iGeom_Instance,
                            iBase_EntityHandle const* face_handles,
                            int face_handles_size,
                            iBase_EntityHandle const* region_handles,
                            int region_handles_size,
                            int** sense,
                            int* sense_allocated,
                            int* sense_size,
                            int* err ){ }
void iGeom_getEgFcSense( iGeom_Instance,
                         iBase_EntityHandle edge,
                         iBase_EntityHandle face,
                         int* sense_out,
                         int* err ){ }
void iGeom_getEgFcArrSense( iGeom_Instance,
                            iBase_EntityHandle const* edge_handles,
                            int edge_handles_size,
                            iBase_EntityHandle const* face_handles,
                            int face_handles_size,
                            int** sense,
                            int* sense_allocated,
                            int* sense_size,
                            int* err ){ }
void iGeom_getEgVtxSense( iGeom_Instance,
                          iBase_EntityHandle edge,
                          iBase_EntityHandle vertex1,
                          iBase_EntityHandle vertex2,
                          int* sense_out,
                          int* err ){ }
void iGeom_getEgVtxArrSense( iGeom_Instance,
                             iBase_EntityHandle const* edge_handles,
                             int edge_handles_size,
                             iBase_EntityHandle const* vertex_handles_1,
                             int veretx_handles_1_size,
                             iBase_EntityHandle const* vertex_handles_2,
                             int vertex_handles_2_size,
                             int** sense,
                             int* sense_allocated,
                             int* sense_size,
                             int* err ){ }
void iGeom_measure( iGeom_Instance,
                    iBase_EntityHandle const* entity_handles,
                    int entity_handles_size,
                    double** measures,
                    int* measures_allocated,
                    int* measures_size,
                    int* err ){ }
void iGeom_getFaceType( iGeom_Instance,
                        iBase_EntityHandle face_handle,
                        char* face_type,
                        int* err,
                        int* face_type_length){ }
void iGeom_getParametric( iGeom_Instance,
                          int* is_parametric,
                          int* err ){ }
void iGeom_isEntParametric( iGeom_Instance,
                            iBase_EntityHandle entity_handle,
                            int* parametric,
                            int* err ){ }
void iGeom_isArrParametric( iGeom_Instance,
                            iBase_EntityHandle const* entity_handles,
                            int entity_handles_size,
                            int** is_parametric,
                            int* is_parametric_allocated,
                            int* is_parametric_size,
                            int* err ){ }
void iGeom_getArrTolerance( iGeom_Instance,
                            iBase_EntityHandle const* entity_handles,
                            int entity_handles_size,
                            double** tolerances,
                            int* tolerances_allocated,
                            int* tolerances_size,
                            int* err ){ }

void iGeom_initEntIter( iGeom_Instance,
                        iBase_EntitySetHandle entity_set_handle,
                        int entity_dimension,
                        iGeom_EntityIterator* entity_iterator,
                        int* err ){ }

void iGeom_initEntArrIter( iGeom_Instance,
                           iBase_EntitySetHandle entity_set_handle,
                           int entity_dimension,
                           int requested_array_size,
                           iGeom_EntityArrIterator* entArr_iterator,
                           int* err ){ }

void iGeom_getNextEntIter( iGeom_Instance,
                           iGeom_EntityIterator,
                           iBase_EntityHandle* entity_handle,
                           int* has_data,
                           int* err ){ }

void iGeom_getNextEntArrIter( iGeom_Instance,
                              iGeom_EntityArrIterator,
                              iBase_EntityHandle** entity_handles,
                              int* entity_handles_allocated,
                              int* entity_handles_size,
                              int* has_data,
                              int* err ){ }

void iGeom_resetEntIter( iGeom_Instance,
                         iGeom_EntityIterator,
                         int* err ){ }

void iGeom_resetEntArrIter( iGeom_Instance,
                            iGeom_EntityArrIterator,
                            int* err ){ }

void iGeom_endEntIter( iGeom_Instance, iGeom_EntityIterator, int* err ){ }

void iGeom_endEntArrIter( iGeom_Instance, iGeom_EntityArrIterator, int* err ){ }

void iGeom_copyEnt( iGeom_Instance,
                    iBase_EntityHandle source,
                    iBase_EntityHandle* copy,
                    int* err ){ }

void iGeom_sweepEntAboutAxis( iGeom_Instance,
			      iBase_EntityHandle geom_entity,
			      double angle,
			      double axis_normal_x,
			      double axis_normal_y,
			      double axis_normal_z,
			      iBase_EntityHandle* geom_entity2,
			      int* err )
{
}

void iGeom_deleteAll(iGeom_Instance instance, int* err )
{
  MBErrorCode rval;
  for (int i = 0; i < 4; i++) {
    rval = MBI->delete_entities(_my_gsets[i]);
    MBERRORR("Failed to delete entities in iGeom_deleteAll.");
    _my_gsets[i].clear();
  }

  RETURN(iBase_SUCCESS);
}

void iGeom_deleteEnt( iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      int* err )
{
  int type;
  iGeom_getEntType(instance, entity_handle, &type, err);
  ERRORR("Failed to get entity type in iGeom_deleteEnt.");

  MBRange::iterator iter = _my_gsets[type].find(MBH_cast(entity_handle));
  if (iter == _my_gsets[type].end()) {
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  _my_gsets[type].erase(iter);

  MBEntityHandle this_entity = MBH_cast(entity_handle);
  MBErrorCode rval = MBI->delete_entities(&this_entity, 1);
  MBERRORR("Failed to delete entity.");
}

void iGeom_createSphere( iGeom_Instance,
			 double radius,
			 iBase_EntityHandle* sphere_handle_out,
			 int* err ){ }

void iGeom_createPrism( iGeom_Instance,
                   double height,
		   int n_sides,
		   double major_rad,
	           double minor_rad,
	           iBase_EntityHandle* prism_handle_out,
                   int* err ){ }


void iGeom_createBrick( iGeom_Instance,
			double x,
			double y,
			double z,
			iBase_EntityHandle* geom_entity,
			int* err ){ }

void iGeom_createCylinder( iGeom_Instance,
			   double height,
			   double major_rad,
			   double minor_rad,
			   iBase_EntityHandle* geom_entity,
			   int* err ){ }

void iGeom_createCone( iGeom_Instance,
		       double height,
		       double major_rad_base,
		       double minor_rad_base,
		       double rad_top,
		       iBase_EntityHandle* geom_entity,
		       int* err ){ }

  void iGeom_createTorus( iGeom_Instance,
                          double major_rad,
                          double minor_rad,
                          iBase_EntityHandle* geom_entity,
                          int* err ){ }

  void iGeom_moveEnt( iGeom_Instance,
                      iBase_EntityHandle geom_entity,
                      double x,
                      double y,
                      double z,
                      int* err ){ }

  void iGeom_rotateEnt( iGeom_Instance,
                        iBase_EntityHandle geom_entity,
                        double angle,
                        double axis_normal_x,
                        double axis_normal_y,
                        double axis_normal_z,
                        int* err ){ }


  void iGeom_reflectEnt( iGeom_Instance,
                         iBase_EntityHandle geom_entity,
                         double plane_normal_x,
                         double plane_normal_y,
                         double plane_normal_z,
                         int* err ){ }

  void iGeom_scaleEnt( iGeom_Instance,
                       iBase_EntityHandle geom_entity,
                       double scale_x,
                       double scale_y,
                       double scale_z,
                       int* err ){ }

  void iGeom_uniteEnts( iGeom_Instance,
                        iBase_EntityHandle const* geom_entities,
                        int geom_entities_size,
                        iBase_EntityHandle* geom_entity,
                        int* err ){ }

  void iGeom_subtractEnts( iGeom_Instance,
                           iBase_EntityHandle blank,
                           iBase_EntityHandle tool,
                           iBase_EntityHandle* geom_entity,
                           int* err ){ }

  void iGeom_intersectEnts( iGeom_Instance,
                            iBase_EntityHandle entity2,
			    iBase_EntityHandle entity1,
			    iBase_EntityHandle* geom_entity,
			    int* err ){ }

  void iGeom_sectionEnt( iGeom_Instance,
                         iBase_EntityHandle geom_entity,
                         double plane_normal_x,
                         double plane_normal_y,
                         double plane_normal_z,
                         double offset,
                         int reverse,
                         iBase_EntityHandle* geom_entity2,
                         int* err ){ }

  void iGeom_imprintEnts( iGeom_Instance,
                          iBase_EntityHandle const* geom_entities,
                          int geom_entities_size,
                          int* err ){ }

  void iGeom_mergeEnts( iGeom_Instance,
                        iBase_EntityHandle const* geom_entities,
                        int geom_entities_size,
                        double tolerance,
                        int* err ){ }

void iGeom_createEntSet(iGeom_Instance instance,
                        int isList,
                        iBase_EntitySetHandle* entity_set_created, 
                        int *err)
{
  iMesh_createEntSet(IMESH_INSTANCE(instance), isList,
		     entity_set_created, err);
  ERRORR("Failed to get create entity set.");
}

void iGeom_destroyEntSet(iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set, 
                         int *err){ }

void iGeom_isList(iGeom_Instance instance,
                  iBase_EntitySetHandle entity_set,
                  int *is_list, 
                  int *err){ }

void iGeom_getNumEntSets(iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set_handle,
                         int num_hops,
                         int *num_sets, 
                         int *err)
{
  iMesh_getNumEntSets(IMESH_INSTANCE(instance), entity_set_handle,
		      num_hops, num_sets, err);
  ERRORR("Failed to get number of entity sets.");
}

void iGeom_getEntSets(iGeom_Instance instance,
                      iBase_EntitySetHandle entity_set_handle,
                      int num_hops,
                      iBase_EntitySetHandle** contained_set_handles,
                      int* contained_set_handles_allocated,
                      int* contained_set_handles_size, 
                      int *err)
{
  iMesh_getEntSets(IMESH_INSTANCE(instance), entity_set_handle,
		   num_hops, contained_set_handles,
		   contained_set_handles_allocated,
		   contained_set_handles_size, err);
  ERRORR("Failed to get entity sets.");
}

void iGeom_addEntToSet(iGeom_Instance instance,
                       iBase_EntityHandle entity_handle,
                       iBase_EntitySetHandle entity_set, 
                       int *err)
{
  iMesh_addEntToSet(IMESH_INSTANCE(instance), entity_handle,
		    entity_set, err);
  ERRORR("Failed to add entity to set.");
}

void iGeom_rmvEntFromSet(iGeom_Instance instance,
                         iBase_EntityHandle entity_handle,
                         iBase_EntitySetHandle entity_set, 
                         int *err)
{
  iMesh_rmvEntFromSet(IMESH_INSTANCE(instance), entity_handle,
		      entity_set, err);
  ERRORR("Failed to remove entity from set.");
}

void iGeom_addEntArrToSet(iGeom_Instance instance,
                          const iBase_EntityHandle* entity_handles,
                          int entity_handles_size,
                          iBase_EntitySetHandle entity_set, 
                          int *err)
{
  iMesh_addEntArrToSet(IMESH_INSTANCE(instance), entity_handles,
		       entity_handles_size, entity_set, err);
  ERRORR("Failed to add entities to set.");
}


void iGeom_rmvEntArrFromSet(iGeom_Instance instance,
                            const iBase_EntityHandle* entity_handles,
                            int entity_handles_size,
                            iBase_EntitySetHandle entity_set,
                            int *err)
{
  iMesh_rmvEntArrFromSet(IMESH_INSTANCE(instance), entity_handles,
			 entity_handles_size, entity_set, err);
  ERRORR("Failed to remove entities from set.");
}

void iGeom_addEntSet(iGeom_Instance instance,
                     iBase_EntitySetHandle entity_set_to_add,
                     iBase_EntitySetHandle entity_set_handle, 
                     int *err)
{
  iMesh_addEntSet(IMESH_INSTANCE(instance), entity_set_to_add,
		  entity_set_handle, err);
  ERRORR("Failed to add entity set to entity set.");
}


void iGeom_rmvEntSet(iGeom_Instance instance,
                     iBase_EntitySetHandle entity_set_to_remove,
                     iBase_EntitySetHandle entity_set_handle, 
                     int *err)
{
  iMesh_rmvEntSet(IMESH_INSTANCE(instance), entity_set_to_remove,
		  entity_set_handle, err);
  ERRORR("Failed to remove entity set from entity set.");
}

void iGeom_isEntContained(iGeom_Instance instance,
                          iBase_EntitySetHandle containing_entity_set,
                          iBase_EntityHandle contained_entity,
                          int *is_contained, 
                          int *err)
{
  iMesh_isEntContained(IMESH_INSTANCE(instance), containing_entity_set,
		       contained_entity, is_contained, err);
  ERRORR("Failed to check if entity is contained to entity set.");
}

void iGeom_isEntArrContained( iGeom_Instance instance,
                              iBase_EntitySetHandle containing_set,
                              const iBase_EntityHandle* entity_handles,
                              int num_entity_handles,
                              int** is_contained,
                              int* is_contained_allocated,
                              int* is_contained_size,
                              int* err )
{
  iMesh_isEntArrContained(IMESH_INSTANCE(instance),
			  containing_set, entity_handles,
			  num_entity_handles, is_contained,
			  is_contained_allocated, is_contained_size,
			  err );
  ERRORR("Failed to check if entities are contained to entity sets.");
 }

void iGeom_isEntSetContained(iGeom_Instance instance,
                             iBase_EntitySetHandle containing_entity_set,
                             iBase_EntitySetHandle contained_entity_set,
                             int *is_contained, 
                             int *err)
{
  iMesh_isEntSetContained (IMESH_INSTANCE(instance),
			   containing_entity_set, contained_entity_set,
			   is_contained, err);
  ERRORR("Failed to check if entity set is contained to entity set.");
}

void iGeom_addPrntChld(iGeom_Instance instance,
                       iBase_EntitySetHandle parent_entity_set,
                       iBase_EntitySetHandle child_entity_set, 
                       int *err)
{
  iMesh_addPrntChld(IMESH_INSTANCE(instance),
		    parent_entity_set, child_entity_set, err);
  ERRORR("Failed to add parent and child relation of geometry sets.");
}

void iGeom_rmvPrntChld(iGeom_Instance instance,
                       iBase_EntitySetHandle parent_entity_set,
                       iBase_EntitySetHandle child_entity_set, 
                       int *err)
{
  iMesh_rmvPrntChld(IMESH_INSTANCE(instance),
		    parent_entity_set, child_entity_set, err);
  ERRORR("Failed to remove parent and child relation of geometry sets.");
}

void iGeom_isChildOf(iGeom_Instance instance,
                     iBase_EntitySetHandle parent_entity_set,
                     iBase_EntitySetHandle child_entity_set,
                     int *is_child, 
                     int *err)
{
  iMesh_isChildOf(IMESH_INSTANCE(instance),
		  parent_entity_set, child_entity_set, is_child, err);
  ERRORR("Failed to check if there is a parent/child relation.");
 }

void iGeom_getNumChld(iGeom_Instance instance,
                      iBase_EntitySetHandle entity_set,
                      int num_hops,
                      int *num_child, 
                      int *err)
{
  iMesh_getNumChld(IMESH_INSTANCE(instance), entity_set,
		   num_hops, num_child, err);
  ERRORR("Failed to get number of children.");
}

void iGeom_getNumPrnt(iGeom_Instance instance,
                      iBase_EntitySetHandle entity_set,
                      int num_hops,
                      int *num_parent, 
                      int *err)
{
  iMesh_getNumPrnt(IMESH_INSTANCE(instance), entity_set,
		   num_hops, num_parent, err);
  ERRORR("Failed to get number of parents.");
}

void iGeom_getChldn(iGeom_Instance instance,
                    iBase_EntitySetHandle from_entity_set,
                    int num_hops,
                    iBase_EntitySetHandle** entity_set_handles,
                    int* entity_set_handles_allocated,
                    int* entity_set_handles_size, 
                    int *err)
{
  iMesh_getChldn(IMESH_INSTANCE(instance), from_entity_set,
		 num_hops, entity_set_handles, entity_set_handles_allocated,
		 entity_set_handles_size, err);
  ERRORR("Failed to get children.");
}

void iGeom_getPrnts(iGeom_Instance instance,
                    iBase_EntitySetHandle from_entity_set,
                    int num_hops,
                    iBase_EntitySetHandle** entity_set_handles,
                    int* entity_set_handles_allocated,
                    int* entity_set_handles_size, 
                    int *err)
{
  iMesh_getPrnts(IMESH_INSTANCE(instance), from_entity_set,
		 num_hops, entity_set_handles, entity_set_handles_allocated,
		 entity_set_handles_size, err);
  ERRORR("Failed to get parents.");
}

void iGeom_createTag(iGeom_Instance instance,
                     const char* tag_name,
                     int tag_size,
                     int tag_type,
                     iBase_TagHandle* tag_handle, 
                     int *err,
                     int tag_name_len )
{
  
  iMesh_createTag(IMESH_INSTANCE(instance), tag_name, tag_size,
		  tag_type, tag_handle, err, tag_name_len);
  ERRORR("Failure to create tag.");
}

void iGeom_destroyTag(iGeom_Instance instance,
                      iBase_TagHandle tag_handle,
                      int forced, 
                      int *err)
{
  iMesh_destroyTag(IMESH_INSTANCE(instance), tag_handle,
		   forced, err);
  ERRORR("Failure to destroy tag.");
 }

void iGeom_getTagName(iGeom_Instance instance,
                      iBase_TagHandle tag_handle,
                      char *name, 
                      int* err,
                      int name_len)
{
  iMesh_getTagName(IMESH_INSTANCE(instance), tag_handle, name,
		   err, name_len);
  ERRORR("Failure to get tag name.");
}

void iGeom_getTagSizeValues(iGeom_Instance instance,
                            iBase_TagHandle tag_handle,
                            int *tag_size, 
                            int *err)
{
  iMesh_getTagSizeValues(IMESH_INSTANCE(instance), tag_handle,
			 tag_size, err);
  ERRORR("Failure to get numbers tag size.");
}

void iGeom_getTagSizeBytes(iGeom_Instance instance,
                           iBase_TagHandle tag_handle,
                           int *tag_size, 
                           int *err)
{
  iMesh_getTagSizeBytes(IMESH_INSTANCE(instance), tag_handle,
			tag_size, err);
  ERRORR("Failure to get byte tag size.");
}

void iGeom_getTagHandle(iGeom_Instance instance,
                        const char* tag_name,
                        iBase_TagHandle *tag_handle, 
                        int *err,
                        int tag_name_len)
{
  iMesh_getTagHandle(IMESH_INSTANCE(instance), tag_name,
		     tag_handle, err, tag_name_len);
  ERRORR("Failure to get tag name.");
}

void iGeom_getTagType(iGeom_Instance instance,
                      iBase_TagHandle tag_handle,
                      int *tag_type, 
                      int *err)
{
  iMesh_getTagType(IMESH_INSTANCE(instance), tag_handle,
		   tag_type, err);
  ERRORR("Failure to get tag type.");
}

void iGeom_setEntSetData(iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set_handle,
                         iBase_TagHandle tag_handle,
                         const char* tag_value,
                         int tag_value_size, 
                         int *err)
{
  iMesh_setEntSetData(IMESH_INSTANCE(instance), entity_set_handle,
		      tag_handle, tag_value, tag_value_size, err);
  ERRORR("Failure to set arbitrary data to entity set.");
}

void iGeom_setEntSetIntData(iGeom_Instance instance,
                            iBase_EntitySetHandle entity_set,
                            iBase_TagHandle tag_handle,
                            int tag_value, 
                            int *err)
{
  iMesh_setEntSetIntData(IMESH_INSTANCE(instance), entity_set,
			 tag_handle, tag_value, err);
  ERRORR("Failure to set integer data to entity set.");
}


void iGeom_setEntSetDblData(iGeom_Instance instance,
                            iBase_EntitySetHandle entity_set,
                            iBase_TagHandle tag_handle,
                            double tag_value, 
                            int *err)
{
  iMesh_setEntSetDblData(IMESH_INSTANCE(instance), entity_set,
			 tag_handle, tag_value, err);
  ERRORR("Failure to set double data to entity set.");
}


void iGeom_setEntSetEHData(iGeom_Instance instance,
                           iBase_EntitySetHandle entity_set,
                           iBase_TagHandle tag_handle,
                           iBase_EntityHandle tag_value, 
                           int *err)
{
  iMesh_setEntSetEHData(IMESH_INSTANCE(instance), entity_set,
			tag_handle, tag_value, err);
  ERRORR("Failure to set entity handle data to entity set.");
}


void iGeom_getEntSetData(iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set_handle,
                         iBase_TagHandle tag_handle,
                         char** tag_value,
                         int* tag_value_allocated,
                         int* tag_value_size, 
                         int *err)
{
  iMesh_getEntSetData (IMESH_INSTANCE(instance), entity_set_handle,
		       tag_handle, tag_value, tag_value_allocated,
		       tag_value_size, err);
  ERRORR("Failure to get arbitrary data from entity set.");
}

void iGeom_getEntSetIntData(iGeom_Instance instance,
                            iBase_EntitySetHandle entity_set,
                            iBase_TagHandle tag_handle,
                            int *out_data, 
                            int *err)
{
  iMesh_getEntSetIntData(IMESH_INSTANCE(instance), entity_set,
			 tag_handle, out_data, err);
  ERRORR("Failure to get integer data from entity set.");
}

void iGeom_getEntSetDblData(iGeom_Instance instance,
                            iBase_EntitySetHandle entity_set,
                            iBase_TagHandle tag_handle,
                            double *out_data, 
                            int *err)
{
  iMesh_getEntSetDblData (IMESH_INSTANCE(instance), entity_set,
			  tag_handle, out_data, err);
  ERRORR("Failure to get double data from entity set.");
}

void iGeom_getEntSetEHData(iGeom_Instance instance,
                           iBase_EntitySetHandle entity_set,
                           iBase_TagHandle tag_handle,
                           iBase_EntityHandle *out_data, 
                           int *err)
{
  iMesh_getEntSetEHData(IMESH_INSTANCE(instance), entity_set,
			tag_handle, out_data, err);
  ERRORR("Failure to get double data from entity set.");
}

void iGeom_getAllEntSetTags(iGeom_Instance instance,
                            iBase_EntitySetHandle entity_set_handle,
                            iBase_TagHandle** tag_handles,
                            int* tag_handles_allocated,
                            int* tag_handles_size, 
                            int *err)
{
  iMesh_getAllEntSetTags(IMESH_INSTANCE(instance), entity_set_handle,
			 tag_handles, tag_handles_allocated,
			 tag_handles_size, err);
  ERRORR("Failure to get double data from entity set.");
}

void iGeom_rmvEntSetTag(iGeom_Instance instance,
                        iBase_EntitySetHandle entity_set_handle,
                        iBase_TagHandle tag_handle, 
                        int *err)
{
  iMesh_rmvEntSetTag (IMESH_INSTANCE(instance), entity_set_handle,
		      tag_handle, err);
  ERRORR("Failure to remove entity set tag.");
}


void iGeom_getArrData(iGeom_Instance instance,
                      const iBase_EntityHandle* entity_handles,
                      int entity_handles_size,
                      iBase_TagHandle tag_handle,
                      char** tag_values,
                      int* tag_values_allocated,
                      int* tag_values_size, 
                      int *err)
{
  iMesh_getArrData (IMESH_INSTANCE(instance), entity_handles,
		    entity_handles_size, tag_handle,
		    tag_values, tag_values_allocated,
		    tag_values_size, err);
  ERRORR("Failure to get tag values of arbitrary type on an array of entities.");
}

void iGeom_getIntArrData(iGeom_Instance instance,
                         const iBase_EntityHandle* entity_handles,
                         int entity_handles_size,
                         iBase_TagHandle tag_handle,
                         int** tag_values,
                         int* tag_values_allocated,
                         int* tag_values_size, 
                         int *err)
{
  iMesh_getIntArrData(IMESH_INSTANCE(instance), entity_handles,
		      entity_handles_size, tag_handle,
		      tag_values, tag_values_allocated,
		      tag_values_size, err);
  ERRORR("Failure to get integer tag values on an array of entities.");
}

void iGeom_getDblArrData(iGeom_Instance instance,
                         const iBase_EntityHandle* entity_handles,
                         int entity_handles_size,
                         iBase_TagHandle tag_handle,
                         double** tag_values,
                         int* tag_values_allocated,
                         int* tag_values_size, 
                         int *err)
{
  iMesh_getDblArrData(IMESH_INSTANCE(instance), entity_handles,
		      entity_handles_size, tag_handle,
		      tag_values, tag_values_allocated,
		      tag_values_size, err);
  ERRORR("Failure to get double tag values on an array of entities.");
}

void iGeom_getEHArrData(iGeom_Instance instance,
                        const iBase_EntityHandle* entity_handles,
                        int entity_handles_size,
                        iBase_TagHandle tag_handle,
                        iBase_EntityHandle** tag_value,
                        int* tag_value_allocated,
                        int* tag_value_size, 
                        int *err)
{
  iMesh_getEHArrData(IMESH_INSTANCE(instance), entity_handles,
		      entity_handles_size, tag_handle,
		      tag_value, tag_value_allocated,
		      tag_value_size, err);
  ERRORR("Failure to get entity handle tag values on an array of entities.");
}

void iGeom_setArrData(iGeom_Instance instance,
                      const iBase_EntityHandle* entity_handles,
                      int entity_handles_size,
                      iBase_TagHandle tag_handle,
                      const char* tag_values,
                      int tag_values_size, 
                      int *err)
{
  iMesh_setArrData(IMESH_INSTANCE(instance), entity_handles,
		   entity_handles_size, tag_handle,
		   tag_values, tag_values_size, err);
  ERRORR("Failure to set tag values of arbitrary type on an array of entities.");
}

void iGeom_setIntArrData(iGeom_Instance instance,
                         const iBase_EntityHandle* entity_handles,
                         int entity_handles_size,
                         iBase_TagHandle tag_handle,
                         const int* tag_values,
                         int tag_values_size, 
                         int *err)
{
  iMesh_setIntArrData(IMESH_INSTANCE(instance), entity_handles,
		      entity_handles_size, tag_handle,
		      tag_values, tag_values_size, err);
  ERRORR("Failure to set interger tag values on an array of entities.");
}

void iGeom_setDblArrData(iGeom_Instance instance,
                         const iBase_EntityHandle* entity_handles,
                         int entity_handles_size,
                         iBase_TagHandle tag_handle,
                         const double* tag_values,
                         const int tag_values_size, 
                         int *err)
{
  iMesh_setDblArrData(IMESH_INSTANCE(instance), entity_handles,
		      entity_handles_size, tag_handle,
		      tag_values, tag_values_size, err);
  ERRORR("Failure to set double tag values on an array of entities.");
}

void iGeom_setEHArrData(iGeom_Instance instance,
                        const iBase_EntityHandle* entity_handles,
                        int entity_handles_size,
                        iBase_TagHandle tag_handle,
                        const iBase_EntityHandle* tag_values,
                        int tag_values_size, 
                        int *err)
{
  iMesh_setEHArrData(IMESH_INSTANCE(instance), entity_handles,
		     entity_handles_size, tag_handle,
		     tag_values, tag_values_size, err);
  ERRORR("Failure to set entity handle tag values on an array of entities.");
}

void iGeom_rmvArrTag(iGeom_Instance instance,
                     const iBase_EntityHandle* entity_handles,
                     int entity_handles_size,
                     iBase_TagHandle tag_handle, 
                     int *err)
{
  iMesh_rmvArrTag(IMESH_INSTANCE(instance), entity_handles,
		  entity_handles_size, tag_handle, err);
  ERRORR("Failure to remove tag values on an array of entities.");
}

void iGeom_getData(iGeom_Instance instance,
                   iBase_EntityHandle entity_handle,
                   iBase_TagHandle tag_handle,
                   char** tag_value,
                   int *tag_value_allocated,
                   int *tag_value_size, 
                   int *err)
{
  iMesh_getData(IMESH_INSTANCE(instance), entity_handle,
		tag_handle, tag_value, tag_value_allocated,
		tag_value_size, err);
  ERRORR("Failure to get tag values of an entity.");
}

void iGeom_getIntData(iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      iBase_TagHandle tag_handle,
                      int *out_data, 
                      int *err)
{
  iMesh_getIntData(IMESH_INSTANCE(instance), entity_handle,
		   tag_handle, out_data, err);
  ERRORR("Failure to get integer tag values of an entity.");
}

void iGeom_getDblData(iGeom_Instance instance,
                      const iBase_EntityHandle entity_handle,
                      const iBase_TagHandle tag_handle,
                      double *out_data, int *err)
{
  iMesh_getDblData(IMESH_INSTANCE(instance), entity_handle,
		   tag_handle, out_data, err);
  ERRORR("Failure to get double tag values of an entity.");
}

void iGeom_getEHData(iGeom_Instance instance,
                     iBase_EntityHandle entity_handle,
                     iBase_TagHandle tag_handle,
                     iBase_EntityHandle *out_data, 
                     int *err)
{
  iMesh_getEHData(IMESH_INSTANCE(instance), entity_handle,
		  tag_handle, out_data, err);
  ERRORR("Failure to get entity handle tag values of an entity.");
}

void iGeom_setData(iGeom_Instance instance,
                   iBase_EntityHandle entity_handle,
                   iBase_TagHandle tag_handle,
                   const char* tag_value,
                   int tag_value_size, 
                   int *err)
{
  iMesh_setData(IMESH_INSTANCE(instance), entity_handle,
		tag_handle, tag_value, tag_value_size, err);
  ERRORR("Failure to set tag values of an entity.");
}

void iGeom_setIntData(iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      iBase_TagHandle tag_handle,
                      int tag_value, 
                      int *err)
{
  iMesh_setIntData(IMESH_INSTANCE(instance), entity_handle,
		   tag_handle, tag_value, err);
  ERRORR("Failure to set integer tag values of an entity.");
}

void iGeom_setDblData(iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      iBase_TagHandle tag_handle,
                      double tag_value, 
                      int *err)
{
  iMesh_setDblData(IMESH_INSTANCE(instance), entity_handle,
		   tag_handle, tag_value, err);
  ERRORR("Failure to set double tag values of an entity.");
}

void iGeom_setEHData(iGeom_Instance instance,
                     iBase_EntityHandle entity_handle,
                     iBase_TagHandle tag_handle,
                     iBase_EntityHandle tag_value, 
                     int *err)
{
  iMesh_setEHData(IMESH_INSTANCE(instance), entity_handle,
		  tag_handle, tag_value, err);
  ERRORR("Failure to set entity handle tag values of an entity.");
}

void iGeom_getAllTags(iGeom_Instance instance,
                      iBase_EntityHandle entity_handle,
                      iBase_TagHandle** tag_handles,
                      int* tag_handles_allocated,
                      int* tag_handles_size, 
                      int *err)
{
  iMesh_getAllTags(IMESH_INSTANCE(instance), entity_handle,
		   tag_handles, tag_handles_allocated,
		   tag_handles_size, err);
  ERRORR("Failure to get all tags.");
}

void iGeom_rmvTag(iGeom_Instance instance,
                  iBase_EntityHandle entity_handle,
                  iBase_TagHandle tag_handle, 
                  int *err)
{
  iMesh_rmvTag(IMESH_INSTANCE(instance), entity_handle,
	       tag_handle, err);
  ERRORR("Failure to remove tag.");
}

void iGeom_subtract(iGeom_Instance instance,
                    iBase_EntitySetHandle entity_set_1,
                    iBase_EntitySetHandle entity_set_2,
                    iBase_EntitySetHandle* result_entity_set, 
                    int *err){ }

void iGeom_intersect(iGeom_Instance instance,
                     iBase_EntitySetHandle entity_set_1,
                     iBase_EntitySetHandle entity_set_2,
                     iBase_EntitySetHandle* result_entity_set, 
                     int *err){ }

void iGeom_unite(iGeom_Instance instance,
                 iBase_EntitySetHandle entity_set_1,
                 iBase_EntitySetHandle entity_set_2,
                 iBase_EntitySetHandle* result_entity_set, 
                 int *err){ }

static inline void
iGeom_processError( iBase_ErrorType code, const char* desc ) 
{
  std::strncpy( iGeom_LAST_ERROR.description, desc,
                sizeof(iGeom_LAST_ERROR.description) );
  iGeom_LAST_ERROR.error_type = code;
}

static void
iGeom_get_adjacent_entities(iGeom_Instance instance,
			    const MBEntityHandle from, 
			    const int to_dim,
			    MBRange &adjs,
			    //std::vector<MBEntityHandle>& adjs,
			    int* err)
{
  int this_dim = -1;
  for (int i = 0; i < 4; i++) {
    if (_my_gsets[i].find(from) != _my_gsets[i].end()) {
      this_dim = i;
      break;
    }
  }
  
  // check target dimension
  if (-1 == this_dim) {
    iGeom_processError(iBase_FAILURE, "Entity not a geometry entity.");
    RETURN(iBase_FAILURE);
  }
  else if (0 > to_dim || 3 < to_dim) {
    iGeom_processError(iBase_FAILURE, "To dimension must be between 0 and 3.");
    RETURN(iBase_FAILURE);
  }
  else if (to_dim == this_dim) {
    iGeom_processError(iBase_FAILURE, "To dimension must be different from entity dimension.");
    RETURN(iBase_FAILURE);
  }
  
  MBErrorCode rval;
  if (to_dim > this_dim) {
    rval = MBI->get_parent_meshsets(from, adjs, to_dim - this_dim);
  }
  else {
    rval = MBI->get_child_meshsets(from, adjs, to_dim - this_dim);
  }

  RETURN(iBase_SUCCESS);
}


