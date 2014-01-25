
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"

using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)


#ifdef MESHDIR
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.sat";
#else
static const char input_cube[] = "/io/cube.sat";
#endif

// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// List of tests in this file
void read_cube_tris_test();
void read_cube_connectivity_test();

//Function used to match triangle connectivity and verts 
void match_tri_connectivity( Range connectivity, 
                             std::vector<EntityHandle> &reference_verts);


int main(int /* argc */, char** /* argv */)
{
  int result = 0;
 
  RUN_TEST(read_cube_connectivity_test);
 
  return result;
}



void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

void read_cube_tris_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
   
  int number_of_tris;

  rval = mb->get_number_entities_by_type(0, MBTRI , number_of_tris);
  std::cout << "Number of Triangles = " << number_of_tris << std::endl;
  CHECK_ERR(rval);

  CHECK_EQUAL( 12, number_of_tris);  

}


void read_cube_vertex_pos_test()
{
  
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //First check that the correct number of vertices are present
  int number_of_verts;
  rval = mb->get_number_entities_by_type( 0, MBVERTEX, number_of_verts );
  CHECK_ERR( rval );

  CHECK_EQUAL( 8, number_of_verts );

  //Retrieve all vertex handles from the mesh
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR( rval );

  //Get the vertex coordinates
  double x[verts.size()];
  double y[verts.size()];
  double z[verts.size()];
  rval = mb-> get_coords( verts, &x[0], &y[0], &z[0] );
  CHECK_ERR( rval );

  //Check against known locations of the vertices

  std::vector<double> x_ref;
  std::vector<double> y_ref;
  std::vector<double> z_ref;

  // Vertex 1

  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 2
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 3
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 4
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 5
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );

  // Vertex 6
  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );

  // Vertex 7
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );
 
  // Vertex 8
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );
 

  std::cout << verts.size() << std::endl;
  std::cout << x_ref.size() << std::endl;
  
  for (unsigned int i=0; i<verts.size(); i++)
    {
      for (unsigned int j=0; j<x_ref.size(); j++)
	{
	  if( x[i]==x_ref[j] && y[i]==y_ref[j] && z[i]==z_ref[j])
            {
              x_ref.erase( x_ref.begin()+j );
              y_ref.erase( y_ref.begin()+j );
              z_ref.erase( z_ref.begin()+j );
              
            }
	}
    }
  
  int leftovers = x_ref.size();
  CHECK_EQUAL(0, leftovers );

}

 
void read_cube_connectivity_test()
{

   ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //Get all vertex handles from the mesh
  std::vector<EntityHandle> verts; 
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts);
  CHECK_ERR(rval);
  
  //Duplicate the vertex handles to match the number of tris
  //they should be connected to upon a correct load
  std::vector<EntityHandle>  ref_verts=verts;
  int copy_count = 0;
  while (copy_count<4)
    {
      copy(verts.begin(),verts.end(),back_inserter(ref_verts));
      copy_count++;
    }
  
  //Make sure everything duplicated without a problem
  if( ref_verts.size()!=40 ) CHECK_ERR(MB_FAILURE);

  //Get all triangles from the mesh
  Range tris;
  rval = mb->get_entities_by_type( 0, MBTRI, tris);
  CHECK_ERR(rval);

 
  //Test connectivity
  for (Range::const_iterator i=tris.begin(); i!=tris.end(); i++)
    {
     //Get the connectivity of a triangle
      Range connect;
      rval = mb->get_connectivity( &(*i), 1, connect);
      CHECK_ERR(rval);
     
      match_tri_connectivity( connect, ref_verts);
      std::cout << "ref_verts size = " << ref_verts.size() << std::endl;
    }
  
    
}

void match_tri_connectivity( Range connectivity, std::vector<EntityHandle> &reference_verts)
{

  for(Range::const_iterator i=connectivity.begin(); i!=connectivity.end(); i++)
    {

      for(unsigned int j=0; j<reference_verts.size(); j++)
	{

	  std::cout << "Triangle Vert Handle = " << *i << std::endl;
          if( *i == reference_verts[j]) reference_verts.erase(reference_verts.begin()+j);
             
	}
    }
}
