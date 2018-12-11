/*
 * addncdata.cpp
 * this tool will take an existing h5m file and add data from some chunk description files
 *  generated from e3sm runs;
 * will support mainly showing the chunks in ViSit
 *
 * example of usage:
 * ./mbaddchunk -i penta3d.h5m -n chunks_on_proc.txt -o penta3d_ch.h5m
 *
 *
 * Basically, will output a new h5m file (penta3d_ch.h5m), which has 2 extra tags, corresponding to the chunks number
 * and processors it sits on
 *
 *
 *  file  penta3d.h5m is obtained from the pentagons file, with a python script
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

  std::string inputfile("penta3d.h5m"), outfile("penta3d_ch.h5m"),  chunkfile_name, gsmapfile;

  opts.addOpt<std::string>("input,i", "input mesh filename", &inputfile);
  opts.addOpt<std::string>("chunkFile,n", "chunk file from cam run", &chunkfile_name);
  opts.addOpt<std::string>("gsMAPfile,g", "gsmap file", &gsmapfile);

  opts.addOpt<std::string>("output,o", "output mesh filename", &outfile);


  opts.parseCommandLine(argc, argv);

  ErrorCode rval;
  Core *mb = new Core();

  rval = mb->load_file(inputfile.c_str()); MB_CHK_SET_ERR(rval, "can't load input file");

  std::cout << " opened " << inputfile << " with initial h5m data.\n";
  // open the netcdf file, and see if it has that variable we are looking for

  Range nodes;
  rval = mb->get_entities_by_dimension(0, 0, nodes);MB_CHK_SET_ERR(rval, "can't get nodes");

  Range edges;
  rval = mb->get_entities_by_dimension(0, 1, edges);MB_CHK_SET_ERR(rval, "can't get edges");

  Range cells;
  rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_SET_ERR(rval, "can't get cells");

  std::cout << " it has " << nodes.size() << " vertices " << edges.size() << " edges " << cells.size() << " cells\n";

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

  if (chunkfile_name.length()>0)
  {

    // Open chunk file
    ifstream inFile;

    inFile.open(chunkfile_name.c_str() );
    if (!inFile) {
        cout << "Unable to open chunk file";
        exit(1); // terminate with error
    }
    Tag pTag, cTag;
    int def_val=-1;
    rval = mb->tag_get_handle("ProcID", 1, MB_TYPE_INTEGER, pTag,
            MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define processor tag");
    rval = mb->tag_get_handle("ChunkID", 1, MB_TYPE_INTEGER, cTag,
          MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define chunk tag");

    int proc, lcid, ncols ;
    while (inFile >> proc) {
      inFile >> lcid >> ncols;
      int Gid;
      for (i=0; i<ncols; i++)
      {
        inFile >> Gid;
        EntityHandle cell = cGidHandle[Gid];
        rval = mb->tag_set_data(pTag, &cell, 1, &proc);MB_CHK_SET_ERR(rval, "can't set proc tag");
        rval = mb->tag_set_data(cTag, &cell, 1, &lcid);MB_CHK_SET_ERR(rval, "can't set chunk tag");

      }
    }

    inFile.close();
  }

  if (gsmapfile.length()>0)
  {

    // Open chunk file
    ifstream inFile;

    inFile.open(gsmapfile.c_str() );
    if (!inFile) {
        cout << "Unable to open gsmap file";
        exit(1); // terminate with error
    }
    Tag pTag, cTag;
    int def_val=-1;
    std::string  procTagName = gsmapfile+"_proc";
    rval = mb->tag_get_handle(procTagName.c_str(), 1, MB_TYPE_INTEGER, pTag,
            MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define processor tag");
    std::string  segTagName = gsmapfile+"_seg";
    rval = mb->tag_get_handle(segTagName.c_str(), 1, MB_TYPE_INTEGER, cTag,
          MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define segment tag");

    int compid, ngseg, gsize;
    inFile >> compid >> ngseg >> gsize;
    for ( i=1; i <= ngseg; i++)
    {
      int start, len, pe;
      inFile >> start >>  len >> pe;
      int Gid;
      for (int j=0; j<len; j++)
      {
        Gid = start + j;
        EntityHandle cell = cGidHandle[Gid];
        rval = mb->tag_set_data(pTag, &cell, 1, &pe);MB_CHK_SET_ERR(rval, "can't set proc tag");
        rval = mb->tag_set_data(cTag, &cell, 1, &i);MB_CHK_SET_ERR(rval, "can't set segment tag");
      }
    }

    inFile.close();
  }

  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "can't write file");
  std::cout << " wrote file " << outfile << "\n";
  return 0;
}
