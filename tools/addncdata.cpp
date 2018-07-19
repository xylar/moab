/*
 * addncdata.cpp
 * this tool will take an existing h5m file and add data from an nc type file
 * will support mainly showing the data associated with unstructured meshes (Homme, MPAS) with Visit
 *
 * example of usage:
 * ./mbaddnc -i wholeFineATM.h5m -n surfdata_ne11np4_simyr1850_c160614.nc -o whole_LONGXY_surfdata.h5m -v LONGXY
 *
 * Basically, will output a new h5m file (whole_LONGXY_surfdata.h5m), which has an extra tag, corresponding to the variable
 *   LONGXY from the file surfdata_ne11np4_simyr1850_c160614.nc; matching is based on the global ids between what we think is the order
 *   on the original file (wholeFineATM.h5m) and the order of surfdata_ne11np4_simyr1850_c160614.nc
 *
 *  file  wholeFineATM.h5m is obtained from a coupled run in e3sm, with the ne 11, np 4,
 *
 *  add an option to output the nc data to the original coarse ATM file, the one that also has the GLOBAL_DOFS
 *   tag with the global DOFs of the GLL points
 */


#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include "netcdf.h"
#include <math.h>
#include <sstream>

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
        MB_SET_ERR(MB_FAILURE, "addncdata:: Couldn't get dimension length"); \
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
          MB_SET_ERR(MB_FAILURE, "addncdata:: Couldn't get variable dimension IDs"); \
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
        MB_SET_ERR(MB_FAILURE, "addncdata:: Couldn't get dimension length"); \
      } \
      vals.resize(ntmp); \
      size_t ntmp1 = 0; \
      ivfail = nc_get_vara_int(ncFile, id, &ntmp1, &ntmp, &vals[0]); \
      if (NC_NOERR != ivfail) { \
        MB_SET_ERR(MB_FAILURE, "addncdata:: Problem getting variable " << name); \
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
        MB_SET_ERR(MB_FAILURE, "addncdata:: Couldn't get dimension length"); \
      } \
      vals.resize(ntmp); \
      size_t ntmp1 = 0; \
      dvfail = nc_get_vara_double(ncFile, id, &ntmp1, &ntmp, &vals[0]); \
      if (NC_NOERR != dvfail) { \
        MB_SET_ERR(MB_FAILURE, "addncdata:: Problem getting variable " << name); \
      } \
    } \
  }

#define GET_1D_FLT_VAR(name, id, vals) \
  { \
    std::vector<int> dum_dims; \
    GET_VAR(name, id, dum_dims); \
    if (-1 != id) { \
      size_t ntmp; \
      int dvfail = nc_inq_dimlen(ncFile, dum_dims[0], &ntmp); \
      if (NC_NOERR != dvfail) { \
        MB_SET_ERR(MB_FAILURE, "addncdata:: Couldn't get dimension length"); \
      } \
      vals.resize(ntmp); \
      size_t ntmp1 = 0; \
      dvfail = nc_get_vara_float(ncFile, id, &ntmp1, &ntmp, &vals[0]); \
      if (NC_NOERR != dvfail) { \
        MB_SET_ERR(MB_FAILURE, "addncdata:: Problem getting variable " << name); \
      } \
    } \
  }

int main(int argc, char* argv[])
{

  ProgOptions opts;

  std::string inputfile, outfile("out.h5m"), netcdfFile, variable_name, sefile_name;

  opts.addOpt<std::string>("input,i", "input mesh filename", &inputfile);
  opts.addOpt<std::string>("netcdfFile,n", "netcdf file aligned with the mesh input file", &netcdfFile);
  opts.addOpt<std::string>("output,o", "output mesh filename", &outfile);

  opts.addOpt<std::string>("var,v", "variable to extract and add to output file", &variable_name);

  opts.addOpt<std::string>("sefile,s", "spectral elements file (coarse SE mesh)", &sefile_name);
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
  nc_type dataType;
  // read the variable, and set it to the tag
  fail = nc_inq_vartype(ncFile, nc_var, &dataType);
  DataType mbtype = MB_TYPE_DOUBLE ;
  bool float_var = false;
  if (NC_INT == dataType)
    mbtype = MB_TYPE_INTEGER;
  else if (NC_DOUBLE != dataType && NC_FLOAT !=dataType)
    MB_CHK_SET_ERR(MB_FAILURE, "unknown type");

  if (NC_FLOAT ==dataType)
    float_var = true;

  int time_id=-1;
  fail = nc_inq_varid(ncFile, "time", &time_id);
  std::vector<float> times;
  if (NC_NOERR == fail)
  {
    int ii;
    GET_1D_FLT_VAR("time", ii, times);
  }
  Tag newTag;

  if (( dims.size()>=1 &&dims.size()<=2)  && (vertex_data || cell_data) )
  {

    if (dims.size()==2)
    {
      fail = nc_inq_dim(ncFile, dims[otherDim], recname, &size_tag);
    }

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
  else if (( dims.size()==3)  && vertex_data && dimIndex==2 && mbtype == MB_TYPE_DOUBLE) // the last one is the vertex
  {
    // the case when the last dim is ncol (for homme type mesh))
    size_t dim0, dim1; // dim 2 is ncol..
    char recname0[NC_MAX_NAME+1], recname1[NC_MAX_NAME+1];
    fail = nc_inq_dim(ncFile, dims[0], recname0, &dim0);
    fail = nc_inq_dim(ncFile, dims[1], recname1, &dim1);
    std::string name0(recname0);
    std::string name1(recname1);
    std::string timestr("time");
    size_t start[3]={0,0,0};
    size_t count[3]={1,0,0};
    count[1]=dim1;
    count[2]=nodes.size();
    // create a few tags, with name inserted
    if (name0.compare("time")==0)
    {

      std::vector<double> dvalues(dim1*count[2]);
      std::vector<float> fvalues(dim1*count[2]);
      std::vector<double> transp(dim1*count[2]);
      // get the variable time
      for (size_t k=0; k<dim0; k++)
      {
        // create a tag for each time, and
        std::stringstream tag_name;
        tag_name<<variable_name<< "_t"<< times[k];
        std::vector<double>  defvals(dim1, 0.);
        rval = mb->tag_get_handle(tag_name.str().c_str(), (int)dim1, mbtype, newTag,
                MB_TAG_CREAT | MB_TAG_DENSE, &defvals[0]); MB_CHK_SET_ERR(rval, "can't define new tag");
        start[0]=k;

        if (float_var)
        {
          fail = nc_get_vara_float(ncFile, nc_var, start, count, &fvalues[0]);
          // now arrange them in a tag, transpose data
          for (size_t ii=0;ii<dim1; ii++)
            for (size_t j=0; j<count[2]; j++)
            {
              transp[j*dim1+ii] = fvalues[ii*count[2]+j];
            }
        }
        else // double
        {
          fail = nc_get_vara_double(ncFile, nc_var, start, count, &dvalues[0]);
          // now arrange them in a tag, transpose data
          for (size_t ii=0;ii<dim1; ii++)
            for (size_t j=0; j<count[2]; j++)
            {
              transp[j*dim1+ii] = dvalues[ii*count[2]+j];
            }
        }
        for (size_t ii=0; ii<nodes.size(); ii++)
        {
          EntityHandle vh=vGidHandle[ii+1];
          rval = mb->tag_set_data(newTag, &vh, 1, &transp[ii*dim1]); MB_CHK_SET_ERR(rval, "can't set tag on nodes");
        }
      }
    }
  }
  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "can't write file");
  std::cout << " wrote file " << outfile << "\n";

  // now, if s option, load the coarse mesh and put data on each element, according to a matrix
  if (!sefile_name.empty())
  {
    // load the file, check for GLOBAL_DOFS tag, and create a new tag with the data associated
    Core *mb2 = new Core();
    rval = mb2->load_file(sefile_name.c_str()); MB_CHK_SET_ERR(rval, "can't load spectral element file");
    std::cout << " loaded spectral file " << sefile_name << "\n";
    // look for GLOBAL_DOFS tag
    Tag gdofeTag;
    rval = mb2->tag_get_handle("GLOBAL_DOFS", gdofeTag);MB_CHK_SET_ERR(rval, "file does not have GLOBAL_DOFS file");
    int sizeTag;
    rval = mb2->tag_get_length(gdofeTag, sizeTag); MB_CHK_SET_ERR(rval, "can't get size of tag");
    int np = (int) sqrt (1.0*sizeTag);
    std::cout << " size of tag: " << sizeTag << " np = " << np << "\n";
    std::vector<int> gdofs;
    Range cells2;
    rval = mb2->get_entities_by_dimension(0, 2, cells2);MB_CHK_SET_ERR(rval, "can't get cells on spectral mesh");
    gdofs.resize(cells2.size()*sizeTag);
    rval = mb2->tag_get_data(gdofeTag, cells2, &gdofs[0]); MB_CHK_SET_ERR(rval, "can't get global dofs tag");
    // create a new tag for element data arranged as DFIELD

    std::vector<double> dfield;
    dfield.resize(sizeTag, 0.0);
    Tag newTag2;
    rval = mb2->tag_get_handle(variable_name.c_str(), (int)sizeTag, mbtype, newTag2,
         MB_TAG_CREAT | MB_TAG_DENSE, &dfield[0]); MB_CHK_SET_ERR(rval, "can't define new tag");

    int i1 = 0; // index in the gdofs array, per element

    // get the tag values from the other moab core, for newTag
    int dataTagLen;
    rval = mb->tag_get_length(newTag, dataTagLen);MB_CHK_SET_ERR(rval, "can't get size of newTag");
    //
    std::vector<double> oldData;
    oldData.resize(dataTagLen*nodes.size()); //

    // get the "old" values
    rval = mb->tag_get_data(newTag, nodes, &oldData[0]); MB_CHK_SET_ERR(rval, "can't get old values");
    for (Range::iterator it=cells2.begin(); it!=cells2.end(); it++)
    {
      EntityHandle cel = *it;
      // gdofs per element are gdofs[i:i+np*np];
      for (int k=0; k<sizeTag; k++)
      {
        int gdof = gdofs[i1+k];
        EntityHandle node = vGidHandle[gdof];
        int index = nodes.index(node);
        // get last layer of data
        double val = oldData[ index*dataTagLen+dataTagLen-1]; //
        dfield[k] = val;
      }
      i1 = i1 + sizeTag;
      rval = mb2->tag_set_data(newTag2, &cel, 1, &dfield[0]);  MB_CHK_SET_ERR(rval, "can't set new tag");

    }

    // write the appended file with the new field:
    rval = mb2->write_file("atm2.h5m");MB_CHK_SET_ERR(rval, "can't write new spectral file");
    std::cout << " wrote file atm2.h5m \n";
  }
  return 0;
}



