#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "MBCartVect.hpp"
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

using namespace moab;

MBInterface *MBI();
int main(int argc, char *argv[])
{ 
  //  std::string input_filename="../../fng/mcnp_n_impr_fluka.h5m";
  std::string input_filename=argv[1];
  MBEntityHandle loaded_file_set;
  MBErrorCode rval;
  
  if(argc == 1)
    {
      std::cout << "please provide filename" << std::endl;
      return 1;
    }
  
  // create meshset to load file into
  rval = MBI()->create_meshset( MESHSET_SET, loaded_file_set );
  assert( rval == MB_SUCCESS );

  // load file
  rval = MBI()->load_file( input_filename.c_str(), &loaded_file_set );
  assert( rval == MB_SUCCESS );
  std::cout << "file loaded" << std::endl;

  // 
  MBRange eh_vertices;
  MBRange eh_curves;
  MBRange eh_triangles;
  MBRange eh_volumes;

  // get the vertices
  rval = MBI()->get_entities_by_dimension(loaded_file_set,0,eh_vertices);
  assert( rval == MB_SUCCESS );
  std::cout << "There are " << eh_vertices.size() << " vertices" << std::endl;

  rval = MBI()->get_entities_by_dimension(loaded_file_set,1,eh_curves);
  assert( rval == MB_SUCCESS );
  std::cout << "There are " << eh_curves.size() << " curves" << std::endl;

  rval = MBI()->get_entities_by_dimension(loaded_file_set,2,eh_triangles);
  assert( rval == MB_SUCCESS );
  std::cout << "There are " << eh_triangles.size() << " triangles" << std::endl;

  rval = MBI()->get_entities_by_dimension(loaded_file_set,3,eh_volumes);
  assert( rval == MB_SUCCESS );
  std::cout << "There are " << eh_volumes.size() << " volumes" << std::endl;

  //rval = MBI()->get_entities_by_type(MB_VERTEX,&(eh_vertices));

  return 0;
}

MBInterface *MBI() 
{
    static MBCore instance;
    return &instance;
}

