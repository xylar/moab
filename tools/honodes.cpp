#include <iostream>
#include <moab/Core.hpp>
#include "moab/ProgOptions.hpp"

using namespace moab;
using namespace std;
const char BRIEF_DESC[] = "Add higher order nodes to existing mesh, eg. hex8 to hex27 (default). Use available options as desired.";
std::ostringstream LONG_DESC;


int main(int argc, char *argv[])
{
  moab::Core *mb = new moab::Core();
  moab::ErrorCode rval;
  bool edge = false;
  bool face = false;
  bool volume = false;

  LONG_DESC << "mbhonodes tool reads a mesh file and adds higher order nodes." << std::endl
            << "Options to add higher order nodes on all volume or face or edge elements of the mesh are supported" << std::endl
            << "Example1: Use the following command to create hex mid volume nodes, in.h5m is a hex8 mesh" << std::endl
            << " -> mbhonodes -i in.h5m -o o.h5m -f -e" << std::endl
            << "Example2: Use the following command to create hex27 mesh o2.h5m, in.h5m is a hex8 mesh" << std::endl
            << " -> mbhonodes -i in.h5m -o o2.h5m"
            << std::endl;

  ProgOptions opts(LONG_DESC.str(), BRIEF_DESC);
  string inFileName = "";
  opts.addRequiredArg<string>("inFile,i", "Specify the output file name string", &inFileName);
#ifdef MOAB_HAVE_HDF5
  string outFileName = "outfile.h5m";
#else
  string outFileName = "outfile.vtk";
#endif
  opts.addOpt<string>("outFile,o", "Specify the output file name string (default outfile.h5m)", &outFileName);
  opts.addOpt<void>("edge,e", "DO NOT create mid nodes along edge (default=true)", &edge);
  opts.addOpt<void>("face,f", "DO NOT create face mid nodes (default=true)", &face);
  opts.addOpt<void>("volume,v", "DO NOT create volume mid nodes (default=true)", &volume);

  opts.parseCommandLine(argc, argv);

  // load the input file
  rval = mb->load_mesh(inFileName.c_str());
  MB_CHK_SET_ERR(rval, "Failed to write the mesh file");
  std::cout << "Read input mesh file: " << inFileName << std::endl;
  moab::Range entities;
  moab::EntityHandle meshset;

  rval = mb->get_entities_by_type(0, MBHEX, entities);
  MB_CHK_SET_ERR(rval, "Failed to get hex entities");
  rval = mb->create_meshset(MESHSET_SET, meshset);
  MB_CHK_SET_ERR(rval, "Failed to create meshset");
  rval = mb->add_entities(meshset, entities);
  MB_CHK_SET_ERR(rval, "Failed to add entitites to meshset");

  rval = mb->convert_entities(meshset, !edge, !face, !volume);
  MB_CHK_SET_ERR(rval, "Failed to convert to higher dimension entities");

  rval = mb->write_mesh(outFileName.c_str());
  MB_CHK_SET_ERR(rval, "Failed to write the mesh file");
  std::cout << "Wrote mesh file: " << outFileName << std::endl;
  delete mb;
}
