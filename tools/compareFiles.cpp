/*
 * compareFiles.cpp
 * this tool will take two existing h5m files, for the same mesh;
 *  they will both have the same GLOBAL_IDs for the elements, but the entity handles can be
 *  very different (depending on how the mesh was partitioned, and saved in parallel)
 *
 *  will compare then the difference between tags, and store the result on one of the files (saved again)
 *
 *
 * example of usage:
 * ./compareFiles -i file1.h5m -j file2.h5m -n <tag_name>  -o out.file
 *
 *
 * Basically, will output a new h5m file (out.file), which has an extra tags, corresponding to the
 * difference between the 2 values
 *
 */


#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace moab;
using namespace std;

int main(int argc, char* argv[])
{

  ProgOptions opts;

  std::string inputfile1("fTargetIntx.h5m"), inputfile2("ocn_proj.h5m"),  outfile("out.h5m");

  std::string tag_name("a2oTAG_proj");

  opts.addOpt<std::string>("input1,i", "input mesh filename 1", &inputfile1);
  opts.addOpt<std::string>("input2,j", "input mesh filename 2", &inputfile2);
  opts.addOpt<std::string>("tagname,n", "tag to compare", &tag_name);
  opts.addOpt<std::string>("outfile,o", "output file", &outfile);

  opts.parseCommandLine(argc, argv);

  ErrorCode rval;
  Core *mb = new Core();

  rval = mb->load_file(inputfile1.c_str()); MB_CHK_SET_ERR(rval, "can't load input file 1");

  Core *mb2 = new Core();
  rval = mb2->load_file(inputfile2.c_str()); MB_CHK_SET_ERR(rval, "can't load input file 2");

  std::cout << " opened " << inputfile1 << " and " << inputfile2 << " with initial h5m data.\n";
  // open the netcdf file, and see if it has that variable we are looking for

  Range nodes;
  rval = mb->get_entities_by_dimension(0, 0, nodes);MB_CHK_SET_ERR(rval, "can't get nodes");

  Range edges;
  rval = mb->get_entities_by_dimension(0, 1, edges);MB_CHK_SET_ERR(rval, "can't get edges");

  Range cells;
  rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_SET_ERR(rval, "can't get cells");

  std::cout << inputfile1 << " has " << nodes.size() << " vertices " << edges.size() << " edges " << cells.size() << " cells\n";

  // construct maps between global id and handles
  std::map<int, EntityHandle> vGidHandle;
  std::map<int, EntityHandle> eGidHandle;
  std::map<int, EntityHandle> cGidHandle;
  std::vector<int> gids;
  Tag gid;
  rval = mb->tag_get_handle("GLOBAL_ID", gid);MB_CHK_SET_ERR(rval, "can't get global id tag");
  gids.resize(nodes.size());
  rval =  mb->tag_get_data(gid, nodes, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on vertices");
  int i=0;
  for (Range::iterator vit=nodes.begin(); vit!=nodes.end(); vit++)
  {
    vGidHandle[gids[i++]] = *vit;
  }

  gids.resize(edges.size());
  rval =  mb->tag_get_data(gid, edges, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on edges");
  i=0;
  for (Range::iterator vit=edges.begin(); vit!=edges.end(); vit++)
  {
    eGidHandle[gids[i++]] = *vit;
  }

  gids.resize(cells.size());
  rval =  mb->tag_get_data(gid, cells, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on cells");
  i=0;
  for (Range::iterator vit=cells.begin(); vit!=cells.end(); vit++)
  {
    cGidHandle[gids[i++]] = *vit;
  }


  Range nodes2;
  rval = mb2->get_entities_by_dimension(0, 0, nodes2);MB_CHK_SET_ERR(rval, "can't get nodes2");

  Range edges2;
  rval = mb2->get_entities_by_dimension(0, 1, edges2);MB_CHK_SET_ERR(rval, "can't get edges2");

  Range cells2;
  rval = mb2->get_entities_by_dimension(0, 2, cells2);MB_CHK_SET_ERR(rval, "can't get cells2");

  std::cout << inputfile2 << " has " << nodes2.size() << " vertices " << edges2.size() << " edges " << cells2.size() << " cells\n";

  // construct maps between global id and handles
  std::map<int, EntityHandle> vGidHandle2;
  std::map<int, EntityHandle> eGidHandle2;
  std::map<int, EntityHandle> cGidHandle2;

  Tag gid2;
  rval = mb2->tag_get_handle("GLOBAL_ID", gid2);MB_CHK_SET_ERR(rval, "can't get global id tag2");
  gids.resize(nodes2.size());
  rval =  mb2->tag_get_data(gid2, nodes2, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on vertices2");

  i=0;
  for (Range::iterator vit=nodes2.begin(); vit!=nodes2.end(); vit++)
  {
    vGidHandle2[gids[i++]] = *vit;
  }

  gids.resize(edges2.size());
  rval =  mb2->tag_get_data(gid2, edges2, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on edges2");
  i=0;
  for (Range::iterator vit=edges2.begin(); vit!=edges2.end(); vit++)
  {
    eGidHandle2[gids[i++]] = *vit;
  }

  gids.resize(cells2.size());
  rval =  mb2->tag_get_data(gid2, cells2, &gids[0]);MB_CHK_SET_ERR(rval, "can't get global id on cells2");
  i=0;
  for (Range::iterator vit=cells2.begin(); vit!=cells2.end(); vit++)
  {
    cGidHandle2[gids[i++]] = *vit;
  }

  Tag tag;
  rval = mb->tag_get_handle(tag_name.c_str(), tag); MB_CHK_SET_ERR(rval, "can't get tag on file 1");

  int len_tag=0;
  rval = mb->tag_get_length(tag, len_tag); MB_CHK_SET_ERR(rval, "can't get tag length on tag");
  std::cout << "length tag : " << len_tag << "\n";


  if (cells.size() != cells2.size())
  {
    std::cout << " meshes are different between 2 files, cells.size do not agree \n";
    exit(1);
  }
  std::vector<double> vals;
  vals.resize(len_tag * cells.size());
  rval = mb->tag_get_data(tag, cells, &vals[0]); MB_CHK_SET_ERR(rval, "can't get tag data");

  Tag tag2;
  rval = mb2->tag_get_handle(tag_name.c_str(), tag2); MB_CHK_SET_ERR(rval, "can't get tag on file 2");
  std::vector<double> vals2;
  vals2.resize(len_tag * cells2.size());
  rval = mb2->tag_get_data(tag2, cells2, &vals2[0]); MB_CHK_SET_ERR(rval, "can't get tag data on file 2");


  rval = mb->delete_entities(edges); MB_CHK_SET_ERR(rval, "can't delete edges from file 1");

  std::string new_tag_name = tag_name + "_2";
  Tag newTag, newTagDiff;
  double def_val = -1000;
  rval = mb->tag_get_handle(new_tag_name.c_str(), 1, MB_TYPE_DOUBLE, newTag,
          MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define new tag");

  std::string tag_name_diff = tag_name + "_diff";
  rval = mb->tag_get_handle(tag_name_diff.c_str(), 1, MB_TYPE_DOUBLE, newTagDiff,
            MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define new tag diff");
  i=0;
  for (Range::iterator c2it = cells2.begin(); c2it != cells2.end(); c2it++)
  {
    double val2 = vals2[i];
    int id2 = gids[i];
    i++;
    EntityHandle c1 = cGidHandle[id2];

    rval = mb->tag_set_data(newTag, &c1, 1, &val2);  MB_CHK_SET_ERR(rval, "can't set new tag");
    int indx = cells.index(c1);
    double diff = vals[indx] - val2;
    rval = mb->tag_set_data(newTagDiff, &c1, 1, &diff);  MB_CHK_SET_ERR(rval, "can't set new tag");
  }


  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "can't write file");
  std::cout << " wrote file " << outfile << "\n";
  return 0;
}
