
#ifndef EXOII_INTERFACE_HPP
#define EXOII_INTERFACE_HPP

enum ExoIIElementType
{
  EXOII_SPHERE = 0,
  EXOII_SPRING,
  EXOII_BAR, EXOII_BAR2, EXOII_BAR3,
  EXOII_BEAM, EXOII_BEAM2, EXOII_BEAM3,
  EXOII_TRUSS, EXOII_TRUSS2, EXOII_TRUSS3,
  EXOII_TRI, EXOII_TRI3, EXOII_TRI6, EXOII_TRI7,
  EXOII_QUAD, EXOII_QUAD4, EXOII_QUAD5, EXOII_QUAD8, EXOII_QUAD9,
  EXOII_SHEL, EXOII_SHELL4, EXOII_SHELL8, EXOII_SHELL9,
  EXOII_TETRA, EXOII_TETRA4, EXOII_TETRA8, EXOII_TETRA10, EXOII_TETRA14,
  EXOII_PYRAMID, EXOII_PYRAMID5, EXOII_PYRAMID8, EXOII_PYRAMID13, EXOII_PYRAMID18,
  EXOII_WEDGE,
  EXOII_KNIFE,
  EXOII_HEX, EXOII_HEX8, EXOII_HEX9, EXOII_HEX20, EXOII_HEX27,
  EXOII_HEXSHELL,
  EXOII_MAX_ELEM_TYPE
};


#define MAX_STR_LENGTH 33
#define MAX_LINE_LENGTH 80

#include "MBInterface.hpp"

class MB_DLL_EXPORT ExoIIInterface
{
public:
  ExoIIInterface(){}
  virtual ~ExoIIInterface(){}
      
  //! given the element name, return the type
  virtual ExoIIElementType element_name_to_type(const char* name) = 0;

  //! get the element type of the entity; this entity can either be a meshset, 
  //! in which case it will be assumed to be a material set meshset, or an 
  //! individual entity.
  virtual  ExoIIElementType get_element_type(MBEntityHandle entity,
      MBTag mid_nodes_tag, MBTag geom_dimension_tag, MBEntityType indiv_entity_type = MBMAXTYPE) = 0;

  virtual int has_mid_nodes(ExoIIElementType elem_type, int dimension) = 0;
  virtual void has_mid_nodes(ExoIIElementType elem_type, int* array) = 0;

  virtual const char* element_type_name(ExoIIElementType type) = 0;

  //! return the geometric dimension of the specified element type
  virtual int geometric_dimension(const ExoIIElementType elem_type) = 0;
  
};


#endif

