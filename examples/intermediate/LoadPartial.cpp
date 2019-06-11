/** @example LoadPartial.cpp \n
 * \brief Load a part of a file  \n
 * <b>To run</b>: LoadPartial <file> <tag_name> <val1> <val2> ...\n
 *
 * In this example, it is shown how to load only a part of one file; the file must be organized in sets.
 * (cherry-picking only the sets we want)
 * The sets to load are identified by a tag name and the tag values for the sets of interest.
 * This procedure is used  when reading in parallel, as each processor will load only
 * its part of the file, identified either by partition or by material/block sets
 *  by default, this example will load parallel partition sets
 *  with values 1, 2, and 5 from ../MeshFiles/unittest/64bricks_1khex.h5m
 *  The example will always write the output to a file name part.h5m
 */

#include <iostream>
#include <vector>

// Include header for MOAB instance and tag conventions for
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

#define NTAGVALS 3

using namespace moab;
using namespace std;

// Function to parse input parameters
ErrorCode get_file_options(int argc, char **argv,
                           string& filename,
                           string& tagName,
                           vector<int>& tagValues)
{
  // Get mesh filename
  if (argc > 1)
    filename = string(argv[1]);
  else
    filename = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");

  // Get tag selection options
  if (argc > 2)
    tagName = string(argv[2]);
  else
    tagName = "USERTAG";

  if (argc > 3) {
    tagValues.resize(argc-3, 0);
    for (int i=3; i < argc; ++i) tagValues[i-3] = atoi(argv[i]);
  }
  else {
    for (unsigned i=0; i < tagValues.size(); ++i) tagValues[i] = 2*i+1;
  }

  if (argc > 1 && argc < 4) // print usage
    cout << " usage is " << argv[0] << " <file> <tag_name> <value1> <value2> .. \n";
  return MB_SUCCESS;
}


int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  std::string filename, tagname;
  vector<int> tagvals(NTAGVALS); // Allocate for a maximum of 5 tag values
  ErrorCode rval = get_file_options(argc, argv, filename, tagname, tagvals);MB_CHK_ERR(rval);

#ifdef MOAB_HAVE_HDF5
  // This file is in the mesh files directory
  rval = mb->load_file(filename.c_str(),
          0, 0, PARALLEL_PARTITION_TAG_NAME, tagvals.data(), (int)tagvals.size());MB_CHK_SET_ERR(rval, "Failed to read");

  // If HANDLEID tag present, convert to long, and see what we read from file
  Tag handleid_tag;
  rval = mb->tag_get_handle("HANDLEID", handleid_tag);
  if (MB_SUCCESS == rval) {
    // Convert a few values for a few vertices
    Range verts;
    rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_SET_ERR(rval, "Failed to get vertices");
    vector<long> valsTag(verts.size());
    rval = mb->tag_get_data(handleid_tag, verts, &valsTag[0]);
    if (MB_SUCCESS == rval)
      cout << "First 2 long values recovered: " << valsTag[0] << " " << valsTag[1] << "\n";
  }

  rval = mb->write_file("part.h5m");MB_CHK_SET_ERR(rval, "Failed to write partial file");
  cout << "Wrote successfully part.h5m.\n";

#else
  std::cout << "Configure MOAB with HDF5 to build and use this example correctly.\n";
#endif
  delete mb;
  return 0;
}
