#include "GenerateHierarchy.hpp"
#include <iostream>

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
GenerateHierarchy::GenerateHierarchy(Interface *impl, ErrorCode &return_value)
{
 
  if (NULL == impl)
    {
      MBI = new Core();
    }

  MBI = impl;

  return_value = setup();
  
}

//---------------------------------------------------------------------------//
// destructor
//---------------------------------------------------------------------------//
GenerateHierarchy::~GenerateHierarchy()
{
}

ErrorCode GenerateHierarchy::setup()
{
  ErrorCode rval;              

  DAG = new DagMC( MBI );

  // init_obb                                                                 
  rval = DAG->load_existing_contents();         
  MB_CHK_SET_ERR(rval, "Failed to load DAG existing contents.");

  // build obbs
  rval = DAG->setup_obbs();
  MB_CHK_SET_ERR(rval, "Failed to load set up obbs.");
 
  // setup indices
  rval = DAG->setup_indices();
  MB_CHK_SET_ERR(rval, "Failed to load set up indices.");

  // get all tag handles
  rval = get_all_handles();
  MB_CHK_SET_ERR(rval, "Failed to get all tag handles.");

  return MB_SUCCESS;
}

ErrorCode GenerateHierarchy::tear_down()
{
  ErrorCode rval;
  OrientedBoxTreeTool *obbTree = new OrientedBoxTreeTool(MBI, "OBB", false);


  // delete root meshset 
  rval = MBI->delete_entities( &(root), 1); MB_CHK_ERR(rval);

  //delete OBB Trees
  for ( int i = 1 ; i <= DAG->num_entities(3) ; i++ )
    {
      //get obb tree root node
      moab::EntityHandle obb_root;
      rval = DAG->get_root(DAG->entity_by_index(3, i), obb_root); MB_CHK_ERR(rval);
      
      //delete tree
      rval = obbTree->delete_tree(obb_root); MB_CHK_ERR(rval);
    }

  //delete obb tag
  rval = MBI->tag_delete(obb_tree_tag); MB_CHK_ERR(rval);

  delete obbTree;
  
  return MB_SUCCESS;
}

ErrorCode GenerateHierarchy::build_hierarchy() 
{                                                                               
  ErrorCode rval;

  // create root meshset-- this will be top of tree
  std::string meshset_name = "root";
  rval = MBI->create_meshset( MESHSET_SET, root); MB_CHK_ERR(rval);
  rval = MBI->tag_set_data( name_tag, &root, 1, meshset_name.c_str()); MB_CHK_ERR(rval);
  rval = MBI->tag_set_data( obj_name_tag, &root, 1,"ROOT"); MB_CHK_ERR(rval);

  
  // loop over volumes and insert each into tree
  for ( int i = 1 ; i <= DAG->num_entities(3) ; i++ )
   {
     // get eh for vol with index i 
     EntityHandle volume = DAG->entity_by_index(3,i);
     rval = insert_in_tree(volume);
     MB_CHK_SET_ERR(rval, "Failed to insert volume into tree.");
   }      

  return MB_SUCCESS;                                                          
}


ErrorCode GenerateHierarchy::get_all_handles()
{
  ErrorCode rval;

  rval = MBI->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
                                name_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  MB_CHK_ERR(rval);
  
  rval = MBI->tag_get_handle( "OBJECT_NAME", 32, MB_TYPE_OPAQUE,
                               obj_name_tag, MB_TAG_SPARSE|MB_TAG_CREAT); 
  MB_CHK_ERR(rval);
  
  rval = MBI->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
                                geom_tag, MB_TAG_SPARSE|MB_TAG_CREAT,&negone);
  MB_CHK_ERR(rval);
  

  rval = MBI->tag_get_handle( MB_OBB_TREE_TAG_NAME, 1, MB_TYPE_HANDLE, 
                                obb_tree_tag, MB_TAG_DENSE );
  MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

ErrorCode GenerateHierarchy::insert_in_tree(EntityHandle volume)
{
  bool inserted = false;
  EntityHandle current_volume = volume; // volume to be inserted 
  EntityHandle tree_volume; // volume already existing in the tree
  Range children; // child meshsets belonging to queried tree volume 
  Range::const_iterator current_child; 
  EntityHandle parent;// = root;  //initialize parent
  ErrorCode rval;
  Range child_volumes;
  Range::iterator child_volume; 
  Range::iterator it;

  //set first tree volume to check against as the root
  tree_volume = root; 
 
  // while not inserted in tree
  while ( !inserted )
    {
      // if current volume is insde of tree volume
      if ( is_A_in_B(current_volume, tree_volume) )
        {
          parent = tree_volume;  

          // if tree_volume has children then we must test them,
          // (tree_volume will change)
          child_volumes = get_children_by_dimension(tree_volume, 3);
          if (child_volumes.size() > 0 ) 
            {
              tree_volume = pop_next_child(child_volumes);
            }

	  // otherwise current_volume is the only child of the tree volume
          else
            { 
              rval = MBI->add_parent_child(parent, current_volume);
              MB_CHK_SET_ERR(rval, "Failed to add parent-child relationship.");

              inserted = true;
            }
	  }
      // if current volume is not in the tree volume, the converse may be true
      else 
	  {
	    // if the tree volume is inside the current volume
	    if( is_A_in_B(tree_volume, current_volume) )
	      {
	        // reverse their parentage
	        rval = MBI->remove_parent_child(parent, tree_volume);
	        MB_CHK_SET_ERR(rval, "Failed to remove parent-child relationship.");
                rval = MBI->add_parent_child(current_volume, tree_volume);
                MB_CHK_SET_ERR(rval, "Failed to add parent-child relationship.");
	      }
	  
	    if (child_volumes.size() == 0 )
	      {
	        rval = MBI->add_parent_child(parent, current_volume);
                MB_CHK_SET_ERR(rval, "Failed to add parent-child relationship.");
	        inserted = true;
	      }

            else
              {
                tree_volume = pop_next_child(child_volumes);
              }
	}   
    }
  return MB_SUCCESS;
}

// returns a valid point on the surface of volume
//ErrorCode GenerateHierarchy::point_on_surface(EntityHandle vol, std::vector<double> &coordinates)
ErrorCode GenerateHierarchy::point_on_surface(EntityHandle vol, double *coord)
{
  Range child_surfaces, triangles, vertices;
  Range::iterator it;
  ErrorCode rval;
  
  // get surface corresponding to volume, then get the triangles  
  child_surfaces = get_children_by_dimension(vol, 2);
  rval = MBI-> get_entities_by_type(*child_surfaces.begin(), MBTRI, triangles); MB_CHK_ERR(rval);

  // now get 1st triangle vertices
  rval = MBI->get_connectivity(&(*triangles.begin()),1,vertices); MB_CHK_ERR(rval);
  
  // now get coordinates of first vertex
  rval = MBI->get_coords(&(*vertices.begin()), 1 , &(coord[0])); MB_CHK_ERR(rval);
  
  return MB_SUCCESS;
}

// runs DAGMC point_in_vol and to test if vol A is inside vol B
//  returns true or false
bool GenerateHierarchy::is_A_in_B(EntityHandle volume_A, EntityHandle volume_B)
{
  // assert that all tree volumes are inside the rootset 
  if ( volume_B == root )
    {
      return true;
    }
  // A should never be the rootset
  if ( volume_A == root )
    {
      return false;
    }
  
//  std::vector<double> coordinates; // vector of x y z
  double coord[3]; // coord[0] = x, etc.
  int result; // point in vol result; 0=F, 1=T
  ErrorCode rval;
  
  // find coordinates of point on surface of A
  rval = point_on_surface(volume_A, coord); 

  // if point on A is inside vol B, return T; o.w. return F
  rval = DAG->point_in_volume(volume_B, coord, result);
  MB_CHK_SET_ERR(rval, "Failed to complete point in volume query.");
  if (result == 0) 
    {   
     return false;
    }
  else
    {  
     return true;
    }
}  


Range GenerateHierarchy::get_children_by_dimension(EntityHandle parent, int desired_dimension)
{
  Range all_children, desired_children;
  Range::iterator it;
  ErrorCode rval;
  int actual_dimension;

  all_children.clear();
  rval = MBI->get_child_meshsets(parent, all_children);

  for ( it = all_children.begin() ; it != all_children.end() ; ++it)
    {
      rval = MBI->tag_get_data(geom_tag, &(*it), 1, &actual_dimension);
      if ( actual_dimension == desired_dimension )
	  {
          desired_children.insert(*it);
        }
    }

  return desired_children;
  
}

EntityHandle GenerateHierarchy::pop_next_child(Range &child_volumes)
{
  EntityHandle next_child = *(child_volumes.begin());
  child_volumes.erase(child_volumes.begin());
  return next_child;
}

ErrorCode GenerateHierarchy::construct_topology()
{

  ErrorCode rval;
  std::map<EntityHandle,EntityHandle> volume_surface; //map of volume
                                                      // to its surface
  for ( int i = 1; i <= DAG->num_entities(3) ; i++ )
    {
      //get the EntityHandle of each volume
      EntityHandle volume = DAG->entity_by_index(3,i);
     
      //get the surface corresponding to each volume
      // at this point, each volume meshset only has one 'child' surface
      // which exactly corresponds to that volume
      Range child_surfaces;
      Range::const_iterator surface;
      child_surfaces = get_children_by_dimension(volume, 2);
      volume_surface[volume]=*child_surfaces.begin();
    }
 

  myGeomTool = new GeomTopoTool(MBI);
  //for each original volume, get its child volumes
  for ( int i = 1 ; i <= DAG->num_entities(3) ; i++)
    {
      EntityHandle volume = DAG->entity_by_index(3,i);
      Range volume_children = get_children_by_dimension(volume, 3);
      
      if (volume_children.size() !=0)
        {
          //loop over all of original volume's child volumes
          for ( Range::iterator j = volume_children.begin() ; j != volume_children.end() ; ++j )
            {
              //set the sense of the surface mapped to the child volume to REVERSE
              // wrt the parent volume
              rval = myGeomTool->set_sense(volume_surface[*j], volume, SENSE_REVERSE);
              MB_CHK_SET_ERR(rval, "Failed to set sense.");
          
              //add the child volume's surface as a child of the original volume
              // and delete the child volume as a child of original volume
              rval = MBI->add_parent_child(volume,volume_surface[*j]);
	      MB_CHK_SET_ERR(rval, "Failed to add parent-child relationship.");
              rval = MBI->remove_parent_child(volume,*j);
	      MB_CHK_SET_ERR(rval, "Failed to remove parent-child relationship.");
            }
        }

    }

  return MB_SUCCESS;

}



