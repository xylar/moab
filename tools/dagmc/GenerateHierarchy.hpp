#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "DagMC.hpp"
#define CHKERR if(MB_SUCCESS != rval){std::cout << "moab error " << rval << std::endl; return rval;} 
#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define MB_OBB_TAG_NAME "OBB"
#include <string>

using namespace moab;

//===========================================================================//
/**
 * \class GenerateHierarchy
 * \brief 

 * This class reads in an geometry file (produced by Read_OBJ) that is 
 * composed of two meshsets for each object; a surface and a volume meshset. 
 * The surface meshset is a child of the volume meshset.
 * There is no sense of hierarchy amongst the objects.
 * Each surface is composed of triangles.  Each triagle has three vertices.
 * 
 * This class has two main, public-facing functions: 
 *     * build_hierarchy
 *     * construct_topology
 *
 * Build_hierarchy will test every object in the geometry and decide where it
 * belongs in a hierarchical tree.  For example, if object A is inside 
 * object B, A will be a child of B and below it in the tree.  If A is neither
 * inside nor outside of B, A and B will share a parent and be on the 
 * same level of the tree.
 * Coming in: each volume meshset has exactly one child, it's corresponding
 * surface meshset.
 * Going out: a root meshset, that contains all volume meshsets, has been
 * created. Each volume meshset will have exactly one 
 * child surface meshset and may have several child volume meshsets.
 * 
 * After the hierarchical structure has been established, DAGMC-appropriate 
 * volumes are created by construct_topology.
 * Coming in: A hierarchical structure of the geometry exists via parent-child 
 * linkages between volume meshsets. The FORWARD surface sense is already set.
 * Each surface has two volumes associated with it; one inside and one outside. 
 * Ex: Vol A has corresponding Surf A. Setting the forward sense places
 * Vol A as one of the two volumes associated with Surf A 
 * Going out: Old parent-child linkages between volumes are deleted and new links
   exist between parent volume meshsets and the surface meshsets of their
   children.  The REVERSE surface sense has also been set.
 */
//===========================================================================//

class GenerateHierarchy {

public:

  
  // constructor will call setup
  GenerateHierarchy(Interface *impl, ErrorCode &return_value); 
 
  // default destructor
  ~GenerateHierarchy(); 

  // creates parent-child relationships
  //  between volumes based on results
  //  from DAGMC point_in_vol 
  ErrorCode build_hierarchy(); 
 
  // sets the surface sense wrt 
  //  parent volume
  ErrorCode construct_topology();
 
  
  ErrorCode tear_down(); 
  
private: // functions

  // this loads the file and sets up DAG
  ErrorCode setup();

  ErrorCode get_all_handles();

  // function that inserts a moab entity into its proper 
  //  postion in the hierarhical tree
  ErrorCode insert_in_tree(EntityHandle volume);

  // gets an entity's children based on the desired dimension
  Range get_children_by_dimension(EntityHandle parent, int desired_dimension);

  // gets the next child from a range of children, 
  //  then erases it from the range
  EntityHandle pop_next_child(Range &child_volumes);

  // runs point in vol to test if a point on vol A is inside vol B or not
  bool is_A_in_B(EntityHandle volume_A, EntityHandle volume_B);

  // returns a valid point on the surface of volume
  std::vector<double>  point_on_surface(EntityHandle vol);

  //ErrorCode point_on_surface(EntityHandle vol, std::vector<double> &coordinates);
  ErrorCode point_on_surface(EntityHandle vol, double *coord);

private: // member variables

  Interface *MBI; // MOAB interface 
  GeomTopoTool *myGeomTool; //
  DagMC *DAG ; // DAGMC instance that is created
  EntityHandle root; //root of the tree
 


  // tags are global variables
  Tag name_tag;
  Tag obj_name_tag;
  Tag geom_tag;
  Tag sense_tag;
  Tag obb_tree_tag,obb_tag; 
  int negone = -1;

};
