/*
 * process_arm.cpp
 *
 *  Created on: April 18, 2013
 */

// process the files from Mark; also, link against netcdf directly, because we will use
// netcdf calls to read data, the same as in ReadNC and ReadNCDF
//
//
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "ReadUtilIface.hpp"

#include "CslamUtils.hpp"

#include "netcdf.h"

#include <algorithm>
#include <string>
#include <assert.h>

#include <cmath>


using namespace moab;

#define INS_ID(stringvar, prefix, id) \
          sprintf(stringvar, prefix, id)

#define GET_DIM(ncdim, name, val)\
    {                            \
    int gdfail = nc_inq_dimid(ncFile, name, &ncdim);          \
    if (NC_NOERR == gdfail) {                                             \
      size_t tmp_val;                                                   \
      gdfail = nc_inq_dimlen(ncFile, ncdim, &tmp_val);                        \
      if (NC_NOERR != gdfail) {                                           \
        std::cout<<"couldn't get dimension length for" << name<< " \n"; \
        return 1;                                              \
      }                                                                 \
      else val = tmp_val;                                               \
    } else val = 0;}

#define GET_DIMB(ncdim, name, varname, id, val) \
          INS_ID(name, varname, id); \
          GET_DIM(ncdim, name, val);

#define GET_VAR(name, id, dims) \
    {                           \
    id = -1;\
    int gvfail = nc_inq_varid(ncFile, name, &id);   \
    if (NC_NOERR == gvfail) {       \
    int ndims;\
    gvfail = nc_inq_varndims(ncFile, id, &ndims);\
    if (NC_NOERR == gvfail) {\
    dims.resize(ndims);    \
    gvfail = nc_inq_vardimid(ncFile, id, &dims[0]);}}}

#define GET_1D_INT_VAR(name, id, vals) \
    {GET_VAR(name, id, vals);  \
  if (-1 != id) {\
    size_t ntmp;\
    int ivfail = nc_inq_dimlen(ncFile, vals[0], &ntmp);\
    vals.resize(ntmp);\
    size_t ntmp1 = 0;                                                           \
    ivfail = nc_get_vara_int(ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != ivfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable " <<name <<"\n";\
      return 1;}}}


#define GET_1D_DBL_VAR(name, id, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    size_t ntmp;\
    int dvfail = nc_inq_dimlen(ncFile, dum_dims[0], &ntmp);\
    vals.resize(ntmp);\
    size_t ntmp1 = 0;                                                           \
    dvfail = nc_get_vara_double(ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

/* get the variable along an index */
#define GET_2D_FLT_VAR(name, id, index, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    size_t ntmp;\
    int dvfail = nc_inq_dimlen(ncFile, dum_dims[1], &ntmp);\
    vals.resize(ntmp);\
    size_t ntmp1[2] = {index, 0}; \
    if (index>=dum_dims[0]) { \
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n"; \
      return 1; \
    } \
    size_t count[2] ={1, ntmp}; \
    dvfail = nc_get_vara_float(ncFile, id, ntmp1, count, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

int main(int argc, char ** argv)
{



  int num_el = 50000;
  if (argc>1)
  {
    num_el = atoi(argv[1]);
  }
  std::cout << "num_el=" << num_el << "\n";

  Core moab;
  Interface & mb = moab;
  //EntityHandle euler_set;
  ErrorCode rval;
  /*
  // this can be directly imported to moab
  const char *sphere_mesh = "ARM_30_x8.g";
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  if (MB_SUCCESS != rval)
    return 1;

  rval = mb.load_file(sphere_mesh, &euler_set);

  if (MB_SUCCESS != rval)
    return 1;

  mb.list_entities(&euler_set, 1);*/

  // now, open the connectivity file and read the element corners variable
  int ncFile, temp_dim;        // netcdf/exodus file

  // now, load the other file

  // now, open the data file and read the lat, lon and U850 and V850

  const char *data_file = "fc5_arm12.cam2.h2.0004-12-12-00000.nc";

  int fail = nc_open(data_file, 0, &ncFile);
  if (NC_NOWRITE != fail) {
    std::cout<<"ReadNCDF:: problem opening Netcdf/Exodus II "<<data_file <<"\n";
    return 1;
  }
  int ncol;
  GET_DIM(temp_dim, "ncol", ncol);
  std::cout << "ncol:" << ncol << "\n";

  std::vector<double> lat;
  std::vector<double> lon;
  GET_1D_DBL_VAR("lat", temp_dim, lat);

  GET_1D_DBL_VAR("lon", temp_dim, lon);

  std::cout<< " lat, lon 0" << lat[0] << " " << lon[0] << "\n";
  std::cout<< " lat, lon 1" << lat[1] << " " << lon[1] << "\n";
  std::cout << " size: " << lat.size() << "\n";


  // see if the mesh from metadata makes sense
  // create quads with the connectivity from conn array

  ReadUtilIface* readMeshIface;
  mb.query_interface( readMeshIface );
  if (NULL==readMeshIface)
    return 1;
  EntityHandle node_handle = 0;
  std::vector<double*> arrays;
  rval = readMeshIface->get_node_coords(3, ncol,
      1, node_handle, arrays);
  if (MB_SUCCESS!= rval)
    return 1;

  SphereCoords  sp;
  sp.R = 1.;
  Tag gid;
  rval = mb.tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER,
      gid, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS!= rval)
      return 1;
  std::vector<int> gids(ncol);

  double conversion_factor = M_PI/180.;
  for (int k=0; k<ncol; k++)
  {
    sp.lat=lat[k]*conversion_factor;
    sp.lon=lon[k]*conversion_factor;
    CartVect pos=spherical_to_cart (sp);
    arrays[0][k]=pos[0];
    arrays[1][k]=pos[1];
    arrays[2][k]=pos[2];
    gids[k]=k+1;
  }

  Range nodes(node_handle, node_handle+ncol-1);
  rval = mb.tag_set_data(gid, nodes, &gids[0]);
  if (MB_SUCCESS!= rval)
       return 1;

  // read the U850 and V850 variables
  std::vector<float> u850;
  GET_2D_FLT_VAR("U850", temp_dim, 0, u850);
  std::vector<float> v850;
  GET_2D_FLT_VAR("V850", temp_dim, 0, v850);

  std::cout << " U850:" << u850[0] << " " << u850[1] << " " << u850[5] << " "<< u850.size()<<"\n";
  std::cout << " V850:" << v850[0] << " " << v850[1] << " " << v850[5] << " "<< u850.size()<<"\n";
  // ok, use radius as 6371km; not needed
  /*
     *  CartVect res;
    res[0] = sc.R * cos(sc.lat)*cos(sc.lon) ; // x coordinate
    res[1] = sc.R * cos(sc.lat)*sin(sc.lon); // y
    res[2] = sc.R * sin(sc.lat);             // z
    return res;
     */
  std::vector<CartVect> velo850(ncol);
  Tag velotag;
  rval = mb.tag_get_handle("VELO", 3, MB_TYPE_DOUBLE,
        velotag, MB_TAG_DENSE|MB_TAG_CREAT);
  if (MB_SUCCESS!= rval)
    return 1;
  for (int k=0; k<ncol; k++)
  {
    double latRad=lat[k]*conversion_factor;
    double lonRad=lon[k]*conversion_factor;
    CartVect U(-sin(lonRad), cos(lonRad), 0.);
    CartVect V(-sin(latRad)*cos(lonRad), -sin(latRad)*cos(lonRad), cos(latRad));
    velo850[k]=U*u850[k] +V*v850[k];
  }
  rval = mb.tag_set_data(velotag, nodes, &(velo850[0][0]));
  if (MB_SUCCESS!= rval)
      return 1;
  EntityHandle newSet;
  rval = mb.create_meshset(MESHSET_SET, newSet);
  if (MB_SUCCESS != rval)
    return 1;

  // so the nodes will be part
  mb.add_entities(newSet, nodes);



  const char *conn_file = "ARM_30_x8-metadata.g";

  fail = nc_open(conn_file, 0, &ncFile);
  if (NC_NOWRITE != fail)
  {
    std::cout << "ReadNCDF:: problem opening Netcdf/Exodus II " << conn_file
        << "\n";
    return 1;
  }

  // read in initial variables
  int  ncenters, ncorners;
  GET_DIM(temp_dim, "ncenters", ncenters);
  GET_DIM(temp_dim, "ncorners", ncorners);
  std::cout << " ncenters, ncorners:" << ncenters << " " << ncorners << "\n";

  int nc_var;
  //int block_seq_id = 1;
  std::vector<int> dims;

  GET_VAR("element_corners", nc_var, dims);
  std::cout << " element corners dims:" << dims.size() << "\n";
  for (size_t i = 0; i < dims.size(); i++)
    std::cout << " dim " << i << " " << dims[i] << "\n";

  // tmp_ptr is of type int* and points at the same place as conn
  EntityHandle * conn = 0;
  //conn = new int[ncenters * ncorners];

  EntityHandle elh;

  readMeshIface->get_element_connect(
          num_el,
          4,
          MBQUAD,
          1,
          elh,
          conn);

  // read the connetivity into that memory,  this will take only part of the array
  // 1/2 if sizeof(EntityHandle) == 64 bits.
  size_t start[2] = { 0, 0 };
  size_t count[2] = { ncorners, ncenters };

  int* tmp_ptr =new int [ncorners*ncenters];
  fail = nc_get_vara_int(ncFile, nc_var, start, count, tmp_ptr);
  if (NC_NOERR != fail)
    return 1;
  // now move connectivity here, using start handle for new vertices
  // also, permute connectivity; it is
  for (int k=0; k<num_el; k++)
  {
    int start_el = k*4;
    for (int j=0; j<4; j++)
       conn[start_el+j] = node_handle + tmp_ptr[j*ncenters+k]-1;
  }
  std::cout << " conn:" << conn[0] << " " << conn[1] << " " << conn[3]<< "\n";

  Range erange(elh, elh+num_el-1);

  mb.add_entities(newSet, erange);
  std::vector<int> gidels(num_el);
  Tag gid2;
  rval = mb.tag_get_handle("GLOBAL_ID_EL", 1, MB_TYPE_INTEGER,
        gid2, MB_TAG_SPARSE|MB_TAG_CREAT);

  if (MB_SUCCESS != rval)
      return 1;
  for (int k=0; k<num_el; k++)
    gidels[k]=k+1;
  mb.tag_set_data(gid2, erange, &gidels[0]);

  mb.write_file("tt2.h5m", 0, 0, &newSet, 1);

  return 0;
}
