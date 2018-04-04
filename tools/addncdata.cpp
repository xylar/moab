/*
 * addncdata.cpp
 * this tool will take an existing h5m file and add data from an nc type file
 * will support mainly showing the data associated with unstructured meshes (Homme, MPAS) with Visit
 *
 */


#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include "netcdf.h"

using namespace moab;

// copy from ReadNCDF.cpp some useful macros for reading from a netcdf file (exodus?)
// ncFile is an integer initialized when opening the nc file in read mode

int ncFile;

#define INS_ID(stringvar, prefix, id) \
  sprintf(stringvar, prefix, id)

#define GET_DIM(ncdim, name, val) \
  { \
    int gdfail = nc_inq_dimid(ncFile, name, &ncdim); \
    if (NC_NOERR == gdfail) { \
      size_t tmp_val; \
      gdfail = nc_inq_dimlen(ncFile, ncdim, &tmp_val); \
      if (NC_NOERR != gdfail) { \
        MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Couldn't get dimension length"); \
      } \
      else \
        val = tmp_val; \
    } \
    else \
      val = 0; \
  }

#define GET_DIMB(ncdim, name, varname, id, val) \
  INS_ID(name, varname, id); \
  GET_DIM(ncdim, name, val);

#define GET_VAR(name, id, dims) \
  { \
    id = -1; \
    int gvfail = nc_inq_varid(ncFile, name, &id); \
    if (NC_NOERR == gvfail) { \
      int ndims; \
      gvfail = nc_inq_varndims(ncFile, id, &ndims); \
      if (NC_NOERR == gvfail) { \
        dims.resize(ndims); \
        gvfail = nc_inq_vardimid(ncFile, id, &dims[0]); \
        if (NC_NOERR != gvfail) { \
          MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Couldn't get variable dimension IDs"); \
        } \
      } \
    } \
  }

#define GET_1D_INT_VAR(name, id, vals) \
  { \
    GET_VAR(name, id, vals); \
    if (-1 != id) { \
      size_t ntmp; \
      int ivfail = nc_inq_dimlen(ncFile, vals[0], &ntmp); \
      if (NC_NOERR != ivfail) { \
        MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Couldn't get dimension length"); \
      } \
      vals.resize(ntmp); \
      size_t ntmp1 = 0; \
      ivfail = nc_get_vara_int(ncFile, id, &ntmp1, &ntmp, &vals[0]); \
      if (NC_NOERR != ivfail) { \
        MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Problem getting variable " << name); \
      } \
    } \
  }

#define GET_1D_DBL_VAR(name, id, vals) \
  { \
    std::vector<int> dum_dims; \
    GET_VAR(name, id, dum_dims); \
    if (-1 != id) { \
      size_t ntmp; \
      int dvfail = nc_inq_dimlen(ncFile, dum_dims[0], &ntmp); \
      if (NC_NOERR != dvfail) { \
        MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Couldn't get dimension length"); \
      } \
      vals.resize(ntmp); \
      size_t ntmp1 = 0; \
      dvfail = nc_get_vara_double(ncFile, id, &ntmp1, &ntmp, &vals[0]); \
      if (NC_NOERR != dvfail) { \
        MB_SET_ERR(MB_FAILURE, "ReadNCDF:: Problem getting variable " << name); \
      } \
    } \
  }


int main(int argc, char* argv[])
{

  ProgOptions opts;

  std::string inputfile, outfile("out.h5m"), netcdfFile, variable_name;

  opts.addOpt<std::string>("input,i", "input mesh filename", &inputfile);
  opts.addOpt<std::string>("netcdfFile,n", "netcdf file aligned with the mesh input file", &netcdfFile);
  opts.addOpt<std::string>("output,o", "output mesh filename", &outfile);

  opts.addOpt<std::string>("var,v", "variable to extract and add to output file", &variable_name);

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

  // Open netcdf/exodus file
  int fail = nc_open(netcdfFile.c_str(), 0, &ncFile);
  if (NC_NOWRITE != fail) {
    MB_SET_ERR(MB_FILE_DOES_NOT_EXIST, "ReadNCDF:: problem opening Netcdf II file " << netcdfFile);
  }

  std::cout << " opened " << netcdfFile << " with new data \n";
  std::vector<int> dims;
  int nc_var;

  size_t  recs;
  char recname[NC_MAX_NAME+1];

  std::cout << " looking for variable " << variable_name << "\n";
  GET_VAR(variable_name.c_str(), nc_var, dims);
  std::cout << " it has " << dims.size() << " dimensions\n";

  int dimIndex = -1; // index of the dimension of interest
  bool vertex_data = false;
  bool cell_data = false;
  for (size_t j=0; j<dims.size(); j++)
  {
    fail = nc_inq_dim(ncFile, dims[j], recname, &recs);
    std::string name_dim(recname);
    std::cout << " dimension index " << j<< " in file: " << dims[j] << " name: " << name_dim << " recs:" << recs << "\n";
    if (recs == nodes.size())
    {
      dimIndex=j;
      vertex_data = true;
    }
    if (recs == cells.size())
    {
      dimIndex=j;
      cell_data = true;
    }
  }
  int otherDim = 1-dimIndex; //used only if 2 dimensions ; could be 0 or 1;
  size_t size_tag = 1; // good for one dimension
  std::vector<int> evals; // size of size_tag
  std::vector<double> dvals; // size of size_tag
  if (( dims.size()>=1 &&dims.size()<=2)  && (vertex_data || cell_data) )
  {
    nc_type dataType;
    // read the variable, and set it to the tag
    fail = nc_inq_vartype(ncFile, nc_var, &dataType);
    DataType mbtype = MB_TYPE_DOUBLE ;
    if (NC_INT == dataType)
      mbtype = MB_TYPE_INTEGER;
    else if (NC_DOUBLE != dataType)
      MB_CHK_SET_ERR(MB_FAILURE, "unknown type");


    if (dims.size()==2)
    {
      fail = nc_inq_dim(ncFile, dims[otherDim], recname, &size_tag);
    }
    Tag newTag;
    int def_val = 0;
    rval = mb->tag_get_handle(variable_name.c_str(), (int)size_tag, mbtype, newTag,
        MB_TAG_CREAT | MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "can't define new tag");


    if (NC_INT == dataType)
    {
      std::vector<int> vals;
      if (1==dims.size())
      {
        GET_1D_INT_VAR(variable_name.c_str(), dims[0], vals);
        if (vertex_data)
        {
          for (size_t k = 0; k<vals.size(); k++)
          {
            EntityHandle vh=vGidHandle[k+1];
            rval = mb->tag_set_data(newTag, &vh, 1, &vals[k]); MB_CHK_SET_ERR(rval, "can't set tag on vertex");
          }
        }
        else // cell_data
        {
          for (size_t k = 0; k<vals.size(); k++)
          {
            EntityHandle ch=cGidHandle[k+1];  // cell handle
            rval = mb->tag_set_data(newTag, &ch, 1, &vals[k]); MB_CHK_SET_ERR(rval, "can't set tag on cell");
          }
        }
      }
      else // dims.size() == 2
      {
        // Single var for all coords
        size_t start[2] = {0, 0}, count[2] = {1,1};
        count[dimIndex] = recs;
        count[otherDim] = size_tag;
        vals.resize(recs*size_tag);
        fail = nc_get_vara_int(ncFile, nc_var, start, count, &vals[0]);
        evals.resize(size_tag);

        if (vertex_data)
        {
          for (size_t k = 0; k<recs; k++)
          {
            EntityHandle vh=vGidHandle[k+1];
            size_t start_in_vals = k*size_tag, stride = 1;
            for (size_t j=0; j<size_tag; j++)
              evals[j] = vals[start_in_vals+j*stride];
            rval = mb->tag_set_data(newTag, &vh, 1, &evals[0]); MB_CHK_SET_ERR(rval, "can't set tag on vertex");
          }
        }
        else // cell_data
        {
          for (size_t k = 0; k<recs; k++)
          {
            EntityHandle ch=cGidHandle[k+1];  // cell handle
            size_t start_in_vals = k*size_tag, stride = 1;
            for (size_t j=0; j<size_tag; j++)
              evals[j] = vals[start_in_vals+j*stride];
            rval = mb->tag_set_data(newTag, &ch, 1, &evals[0]);MB_CHK_SET_ERR(rval, "can't set tag on cell");
          }
        }
      }
    }
    else
    {
      std::vector<double> vals;
      if (1==dims.size())
      {
        GET_1D_DBL_VAR(variable_name.c_str(), dims[0], vals);
        if (vertex_data)
        {
          for (size_t k = 0; k<vals.size(); k++)
          {
            EntityHandle vh=vGidHandle[k+1];
            rval = mb->tag_set_data(newTag, &vh, 1, &vals[k]); MB_CHK_SET_ERR(rval, "can't set tag on vertex");
          }
        }
        else // cell_data
        {
          for (size_t k = 0; k<vals.size(); k++)
          {
            EntityHandle ch=cGidHandle[k+1];  // cell handle
            rval = mb->tag_set_data(newTag, &ch, 1, &vals[k]); MB_CHK_SET_ERR(rval, "can't set tag on vertex");
          }
        }
      }
      else // dims.size() == 2
      {
        // Single var for all coords
        size_t start[2] = {0, 0}, count[2] = {1,1};
        count[dimIndex] = recs;
        count[otherDim] = size_tag;
        vals.resize(recs*size_tag);
        fail = nc_get_vara_double(ncFile, nc_var, start, count, &vals[0]);
        dvals.resize(size_tag);

        if (vertex_data)
        {
          for (size_t k = 0; k<recs; k++)
          {
            EntityHandle vh=vGidHandle[k+1];
            size_t start_in_vals = k*size_tag, stride = 1;
            for (size_t j=0; j<size_tag; j++)
              dvals[j] = vals[start_in_vals+j*stride];
            rval = mb->tag_set_data(newTag, &vh, 1, &dvals[0]); MB_CHK_SET_ERR(rval, "can't set tag on vertex");
          }
        }
        else // cell_data
        {
          for (size_t k = 0; k<recs; k++)
          {
            EntityHandle ch=cGidHandle[k+1];  // cell handle
            size_t start_in_vals = k*size_tag, stride = 1;
            for (size_t j=0; j<size_tag; j++)
              dvals[j] = vals[start_in_vals+j*stride];
            rval = mb->tag_set_data(newTag, &ch, 1, &dvals[0]);MB_CHK_SET_ERR(rval, "can't set tag on cell");
          }
        }
      }
    }
  }

  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "can't write file");
  std::cout << " wrote file " << outfile << "\n";
  return 0;
}



