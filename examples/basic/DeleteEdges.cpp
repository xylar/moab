/** @example GetEntities.cpp
 * Description: Get entities and report non-vertex entity connectivity and vertex adjacencies.\n
 * then delete edges, and write result
 * To run: ./DeleteEdges [meshfile] [outfile]\n
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string(MESH_DIR) + string("/hex01.vtk");
string out_file = string("outFile.h5m");

int main(int argc, char **argv)
{
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }
  if (argc > 2) {
    // User has specified an output file
    out_file = argv[2];
  }

  // Instantiate & load a mesh from a file
  Core* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  Range edges;
  rval = mb->get_entities_by_dimension(0, 1, edges);MB_CHK_ERR(rval);
  rval = mb->delete_entities(edges); MB_CHK_ERR(rval);

  rval = mb->write_file(out_file.c_str()); MB_CHK_ERR(rval);
  delete mb;

  return 0;
}
