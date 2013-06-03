#include "ReadNC.hpp"
#include "NCHelper.hpp"

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <dirent.h>

#include "moab/Core.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/SpectralMeshTool.hpp"

//#include "bil.h"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

ReaderIface* ReadNC::factory(Interface* iface) {
  return new ReadNC(iface);
}

ReadNC::ReadNC(Interface* impl) :
  mbImpl(impl), CPU_WORD_SIZE(-1), IO_WORD_SIZE(-1), fileId(-1), tMin(-1), tMax(-1), iDim(-1), jDim(-1), tDim(-1), iCDim(-1),
  jCDim(-1), numUnLim(-1), mCurrentMeshHandle(0), startVertex(0), startElem(0), mGlobalIdTag(0), mpFileIdTag(NULL), max_line_length(-1),
  max_str_length(-1), vertexOffset(0), dbgOut(stderr), isParallel(false), partMethod(-1),
#ifdef USE_MPI
  myPcomm(NULL), 
#endif
  noMesh(false), noVars(false), spectralMesh(false), helper(NULL)
{
  assert(impl != NULL);

  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
    gCDims[i] = -1;
    lCDims[i] = -1;
  }

  locallyPeriodic[0] = locallyPeriodic[1] = 0;
  globallyPeriodic[0] = globallyPeriodic[1] = 0;

  impl->query_interface(readMeshIface);
}

void ReadNC::reset() {
  CPU_WORD_SIZE = -1;
  IO_WORD_SIZE = -1;
  fileId = -1;
  tMin = tMax = -1;
  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
    gCDims[i] = -1;
    lCDims[i] = -1;
  }

  iDim = jDim = tDim = -1;
  iCDim = jCDim = -1;
  numUnLim = -1;
  mCurrentMeshHandle = 0;
  startVertex = startElem = 0;
  mGlobalIdTag = 0;
  max_line_length = -1;
  max_str_length = -1;
  vertexOffset = 0;
  dbgOut = stderr;
  mCurrentMeshHandle = 0;
  vertexOffset = 0;

#ifdef USE_MPI
  myPcomm = NULL;
#endif
}

ReadNC::~ReadNC() {
  mbImpl->release_interface(readMeshIface);
  if (helper != NULL)
    delete helper;
}

ErrorCode ReadNC::load_file(const char *file_name, const EntityHandle* file_set, const FileOptions& opts,
    const ReaderIface::SubsetList* /*subset_list*/, const Tag* file_id_tag) {

  ErrorCode rval = MB_SUCCESS;

  //See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;
  /*
  if (file_id_tag)
    mGlobalIdTag = *file_id_tag;
  else {
    */
    //! get and cache predefined tag handles
  int dum_val = 0;
  rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, mGlobalIdTag, MB_TAG_DENSE | MB_TAG_CREAT, &dum_val);
  if (MB_SUCCESS != rval)
    return rval;
  mpFileIdTag = file_id_tag; // store the pointer to the tag ; if not null, set when global id tag
  // is set too, with the same data , duplicated

  std::string partition_tag_name;
  rval = parse_options(opts, var_names, tstep_nums, tstep_vals);
  ERRORR(rval, "Trouble parsing option string.");

  // Open the file
  dbgOut.tprintf(1, "Opening file %s\n", file_name);
  fileName = std::string(file_name);
  int success;

#ifdef PNETCDF_FILE
  if (isParallel)
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), file_name, 0, MPI_INFO_NULL, &fileId);
  else
    success = NCFUNC(open)(MPI_COMM_SELF, file_name, 0, MPI_INFO_NULL, &fileId);
#else
  success = NCFUNC(open)(file_name, 0, &fileId);
#endif

  ERRORS(success, "Trouble opening file.");

  // BIL data

  if (BIL_mode_enabled(file_name)) {

    rval = get_BIL_dir();
    ERRORS(rval, "Failed to find directory with BIL data.");

    dbgOut.tprintf(1, "Reading BIL data from directory: %s\n", BIL_dir.c_str());

    rval = load_BIL(BIL_dir, file_set, opts, file_id_tag);
    ERRORR(rval, "Trouble reading BIL data.");

    return rval;
  }

  // end of BIL

  // Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

  // make sure there's a file set to put things in
  EntityHandle tmp_set;
  if (noMesh && !file_set) {
    ERRORR(MB_FAILURE, "NOMESH option requires non-NULL file set on input.\n");
  }
  else if (!file_set || (file_set && *file_set == 0)) {
    rval = mbImpl->create_meshset(MESHSET_SET, tmp_set);
    ERRORR(rval, "Trouble creating file set.");
  }
  else
    tmp_set = *file_set;

  // get the scd interface
  ScdInterface *scdi = NULL;
  rval = mbImpl->query_interface(scdi);
  if (!scdi)
    return MB_FAILURE;

  // Check if CF convention is being followed
  std::map<std::string, AttData>::iterator attIt = globalAtts.find("conventions");
  if (attIt == globalAtts.end())
    attIt = globalAtts.find("Conventions");

  if (attIt == globalAtts.end()) {
    ERRORR(MB_FAILURE, "File does not have conventions global attribute.");
  }

  unsigned int sz = attIt->second.attLen;
  std::string att_data;
  att_data.resize(sz + 1);
  att_data[sz] = '\000';
  success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
  ERRORS(success, "Failed to read conventions global attribute char data.");
  if (att_data.find("CF") == std::string::npos) {
    ERRORR(MB_FAILURE, "File not following known conventions.");
  }

  if (helper != NULL)
    delete helper;

  helper = NCHelper::get_nc_helper(this, fileId, opts);
  if (helper == NULL) {
    ERRORR(MB_FAILURE, "Failed to get NCHelper class instance.");
  }

  rval = helper->init_mesh_vals(opts, tmp_set);
  if (MB_SUCCESS != rval)
    return rval;

  // Create mesh vertex/quads sequences
  Range quads;
  if (noMesh && !noVars) {
    rval = check_verts_quads(tmp_set);
    ERRORR(rval, "Mesh characteristics didn't match from last read.\n");
  }
  else if (!noMesh) {
    rval = helper->create_verts_quads(scdi, opts, tmp_set, quads);
    if (MB_SUCCESS != rval)
      return rval;
  }

  bool scdMesh = helper->is_scd_mesh();
  if (noMesh && !scdMesh)
  {
    // we need to populate localGid range with the gids of vertices from the tmp_set
    // localGid is important in reading the variable data into the nodes
    // also, for our purposes, localGid is truly the GLOBAL_ID tag data, not other
    // file_id tags that could get passed around in other scenarios for parallel reading
    // for nodal_partition, this local gid is easier, should be initialized with only
    // the owned nodes

    // we need to get all vertices from tmp_set (it is the input set in no_mesh scenario)
    Range local_verts;
    rval = mbImpl->get_entities_by_dimension(tmp_set, 0, local_verts);
    if (MB_FAILURE == rval)
      return rval;
    std::vector<int> gids(local_verts.size());
    // !IMPORTANT : this has to be the GLOBAL_ID tag
    rval=mbImpl->tag_get_data(mGlobalIdTag, local_verts, &gids[0]);
    if (MB_FAILURE == rval)
      return rval;
    // this will do a smart copy
    std::copy(gids.begin(), gids.end(), range_inserter(localGid));
  }

  // Read variables onto grid
  if (!noVars) {
    rval = read_variables(tmp_set, var_names, tstep_nums, scdMesh);
    if (MB_FAILURE == rval)
      return rval;
  }
  else {
    // read dimension variable by default, the ones that are also variables
    std::vector<std::string> filteredDimNames;
    for (unsigned int i = 0; i < dimNames.size(); i++) {
      std::map<std::string, VarData>::iterator mit = varInfo.find(dimNames[i]);
      if (mit != varInfo.end())
        filteredDimNames.push_back(dimNames[i]);
    }
    rval = read_variables(tmp_set, filteredDimNames, tstep_nums, scdMesh);
    if (MB_FAILURE == rval)
      return rval;
  }

#ifdef USE_MPI
  // create partition set, and populate with elements
  if (isParallel) {
    EntityHandle partn_set;
    rval = mbImpl->create_meshset(MESHSET_SET, partn_set);
    ERRORR(rval, "Trouble creating partition set.");

    rval = mbImpl->add_entities(partn_set,quads);
    ERRORR(rval, "Couldn't add new quads to partition set.");

    Range verts;
    rval = mbImpl->get_connectivity(quads, verts);
    ERRORR(rval, "Couldn't get verts of quads");

    rval = mbImpl->add_entities(partn_set,verts);
    ERRORR(rval, "Couldn't add new verts to partition set.");

    myPcomm->partition_sets().insert(partn_set);

    //write partition tag name on partition set
    Tag part_tag;
    rval = mbImpl->tag_get_handle( partitionTagName.c_str(), 1, MB_TYPE_INTEGER, part_tag );
    if (MB_SUCCESS != rval) {
      // fall back to the partition tag
      part_tag = myPcomm->partition_tag();
    }

    int dum_rank = myPcomm->proc_config().proc_rank();
    rval = mbImpl->tag_set_data(part_tag, &partn_set, 1, &dum_rank);
    if (MB_SUCCESS != rval) return rval;
  }
#endif

  mbImpl->release_interface(scdi);
  ERRORR(rval, "Trouble creating scd element sequence.");

  // create nc conventional tags when loading header info only
  if (noMesh && noVars) {
    rval = create_tags(scdi, tmp_set, tstep_nums);
    ERRORR(rval, "Trouble creating nc conventional tags.");
  }
  // close the file
  success = NCFUNC(close)(fileId);
  ERRORS(success, "Trouble closing file.");

  return MB_SUCCESS;
}
    
ErrorCode ReadNC::load_BIL(std::string , const EntityHandle* , const FileOptions& , const Tag* ) {
  /*
   BIL_Init( MPI_COMM_WORLD );

   void ** buffer;

   DIR * dir;
   struct dirent * ent;
   dir = opendir(dir_name.c_str());
   if (dir != NULL) {
   while ((ent = readdir(dir)) != NULL) {
   if (strlen(ent->d_name) > 3) { //filter out . and ..

   dbgOut.tprintf(1,"reading block from %s\n",ent->d_name);

   int num_dims = 3;
   int time_d = 1;
   int lev_d  = 26;
   int ncol_d = 3458;
   int block_start[3] = {0,0,0};
   int block_size[3]  = {time_d, lev_d, ncol_d};
   const char * file_name = ent->d_name;
   const char * var_name = "T";

   BIL_Add_block_nc(num_dims, block_start, block_size,
   file_name, var_name, buffer);
   }
   }
   closedir (dir);
   }

   BIL_Read();

   BIL_Finalize();
   */
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_BIL_dir() {
  std::map<std::string, AttData> dirAtt;
  ErrorCode result = get_attributes(NC_GLOBAL, 1, dirAtt);
  ERRORR(result, "Failed to get BIL_DIR attribute");

  std::string attname;
  std::map<std::string, AttData>::iterator attIt = dirAtt.find("BIL_DIR");

  unsigned int sz = attIt->second.attLen;
  char *att_data = (char *) malloc(sz + 1);
  att_data[sz] = '\000';
  int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (char*) att_data);
  ERRORS(success, "Trouble getting BIL data directory.");
  BIL_dir = std::string(att_data);

  return MB_SUCCESS;
}

bool ReadNC::BIL_mode_enabled(const char * file_name) {
  std::string file_path = std::string(file_name);
  int idx = file_path.find_last_of("/");
  std::string file = file_path.substr(idx + 1);

  return (file == "BIL_DIR.nc");
}

ErrorCode ReadNC::parse_options(const FileOptions &opts, std::vector<std::string> &var_names, std::vector<int> &tstep_nums,
                                std::vector<double> &tstep_vals) {
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NC ");
  }

  ErrorCode rval = opts.get_strs_option("VARIABLE", var_names);
  if (MB_TYPE_OUT_OF_RANGE == rval)
    noVars = true;
  else
    noVars = false;
  opts.get_ints_option("TIMESTEP", tstep_nums);
  opts.get_reals_option("TIMEVAL", tstep_vals);
  rval = opts.get_null_option("NOMESH");
  if (MB_SUCCESS == rval)
    noMesh = true;

  rval = opts.get_null_option("SPECTRAL_MESH");
  if (MB_SUCCESS == rval)
    spectralMesh = true;

  if (2 <= dbgOut.get_verbosity()) {
    if (!var_names.empty()) {
      std::cerr << "Variables requested: ";
      for (unsigned int i = 0; i < var_names.size(); i++)
        std::cerr << var_names[i];
      std::cerr << std::endl;
    }
    if (!tstep_nums.empty()) {
      std::cerr << "Timesteps requested: ";
      for (unsigned int i = 0; i < tstep_nums.size(); i++)
        std::cerr << tstep_nums[i];
      std::cerr << std::endl;
    }
    if (!tstep_vals.empty()) {
      std::cerr << "Time vals requested: ";
      for (unsigned int i = 0; i < tstep_vals.size(); i++)
        std::cerr << tstep_vals[i];
      std::cerr << std::endl;
    }
  }

#ifdef USE_MPI
  //TODO handle options better
  //rval = opts.get_option("PARTITION", partitionTagName);
  /*
   part = (rval != MB_ENTITY_NOT_FOUND);
   rval = opts.match_option("PARALLEL", "READ_PART");
   parread = (rval != MB_ENTITY_NOT_FOUND);
   rval = opts.get_null_option("TRIVIAL_PARTITION");
   partriv = (rval != MB_ENTITY_NOT_FOUND);
   */
  isParallel = (opts.match_option("PARALLEL","READ_PART") != MB_ENTITY_NOT_FOUND);

  if (!isParallel)
  // return success here, since rval still has _NOT_FOUND from not finding option
  // in this case, myPcomm will be NULL, so it can never be used; always check for isParallel 
  // before any use for myPcomm
    return MB_SUCCESS;
  
  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    readMeshIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl, MPI_COMM_WORLD);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);

  const char *part_options[] = {"alljorkori", "alljkbal", "sqij", "sqjk",
      "TRIVIAL_PARTITION" };
  int dum;
  rval = opts.match_option("PARTITION_METHOD", part_options, dum);
  if (rval == MB_FAILURE) {
    readMeshIface->report_error("Unknown partition method specified.");
    partMethod = ScdParData::ALLJORKORI;
  }
  else if (rval == MB_ENTITY_NOT_FOUND)
    partMethod = ScdParData::ALLJORKORI;
  else
    partMethod = dum;
#endif

  return MB_SUCCESS;
}

ErrorCode ReadNC::check_verts_quads(EntityHandle file_set) {
  // check parameters on this read against what was on the mesh from last read
  // get the number of vertices
  int num_verts;
  ErrorCode rval = mbImpl->get_number_entities_by_dimension(file_set, 0, num_verts);
  ERRORR(rval, "Trouble getting number of vertices.");

  // check against parameters
  //int expected_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  //if (num_verts != expected_verts)
  //ERRORR(MB_FAILURE, "Number of vertices doesn't match.");

  // check the number of elements too
  int num_elems;
  rval = mbImpl->get_number_entities_by_dimension(file_set, (-1 == lDims[2] ? 2 : 3), num_elems);
  ERRORR(rval, "Trouble getting number of elements.");

  // check against parameters
  //int expected_elems = (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2]);
  //if (num_elems != expected_elems)
  //ERRORR(MB_FAILURE, "Number of elements doesn't match.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::create_scd_verts_quads(ScdInterface *scdi, EntityHandle tmp_set, Range &quads) {
  Range tmp_range;
  ScdBox *scd_box;

  ErrorCode rval = scdi->construct_box(HomCoord(lDims[0], lDims[1], lDims[2], 1), HomCoord(lDims[3], lDims[4], lDims[5], 1), NULL,
      0, scd_box, locallyPeriodic, &parData);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

  // add box set and new vertices, elements to the file set
  tmp_range.insert(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);
  tmp_range.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements() - 1);
  tmp_range.insert(scd_box->box_set());
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");

  dbgOut.tprintf(1, "scdbox %d quads, %d vertices\n", scd_box->num_elements(), scd_box->num_vertices());

  // get a ptr to global id memory
  void *data;
  int count;
  const Range::iterator topv = tmp_range.upper_bound(tmp_range.begin(), tmp_range.end(), scd_box->start_vertex()
      + scd_box->num_vertices());
  rval = mbImpl->tag_iterate(mGlobalIdTag, tmp_range.begin(), topv, count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == scd_box->num_vertices());
  int *gid_data = (int*) data;

  // set the vertex coordinates
  double *xc, *yc, *zc;
  rval = scd_box->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl, itmp;
  int dil = lDims[3] - lDims[0] + 1;
  int djl = lDims[4] - lDims[1] + 1;
  int di = gDims[3] - gDims[0] + 1;
  int dj = gDims[4] - gDims[1] + 1;
  assert(dil == (int)ilVals.size() && djl == (int)jlVals.size() &&
      (-1 == lDims[2] || lDims[5]-lDims[2]+1 == (int)klVals.size()));
#define INDEX(i,j,k) ()
  for (kl = lDims[2]; kl <= lDims[5]; kl++) {
    k = kl - lDims[2];
    for (jl = lDims[1]; jl <= lDims[4]; jl++) {
      j = jl - lDims[1];
      for (il = lDims[0]; il <= lDims[3]; il++) {
        i = il - lDims[0];
        unsigned int pos = i + j * dil + k * dil * djl;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == lDims[2] ? 0.0 : klVals[k]);
        itmp = (!locallyPeriodic[0] && globallyPeriodic[0] && il == gDims[3] ? gDims[0] : il);
        *gid_data = (-1 != kl ? kl * di * dj : 0) + jl * di + itmp + 1;
        gid_data++;
      }
    }
  }
#undef INDEX

#ifndef NDEBUG
  int num_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  std::vector<int> gids(num_verts);
  Range verts(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);
  rval = mbImpl->tag_get_data(mGlobalIdTag, verts, &gids[0]);
  ERRORR(rval, "Trouble getting gid values.");
  int vmin = *(std::min_element(gids.begin(), gids.end())), vmax = *(std::max_element(gids.begin(), gids.end()));
  dbgOut.tprintf(1, "Vertex gids %d-%d\n", vmin, vmax);
#endif

  // add elements to the range passed in
  quads.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements() - 1);

  if (2 <= dbgOut.get_verbosity()) {
    assert(scd_box->boundary_complete());
    EntityHandle dum_ent = scd_box->start_element();
    rval = mbImpl->list_entities(&dum_ent, 1);
    ERRORR(rval, "Trouble listing first hex.");

    std::vector<EntityHandle> connect;
    rval = mbImpl->get_connectivity(&dum_ent, 1, connect);
    ERRORR(rval, "Trouble getting connectivity.");

    rval = mbImpl->list_entities(&connect[0], connect.size());
    ERRORR(rval, "Trouble listing element connectivity.");
  }

  Range edges;
  mbImpl->get_adjacencies(quads, 1, true, edges, Interface::UNION);

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variable_setup(std::vector<std::string> &var_names, std::vector<int> &tstep_nums,
    std::vector<VarData> &vdatas, std::vector<VarData> &vsetdatas, bool scdMesh) {
  std::map<std::string, VarData>::iterator mit;

  if (scdMesh) { // scd mesh
    // if empty read them all
    if (var_names.empty()) {
      for (mit = varInfo.begin(); mit != varInfo.end(); mit++) {
        VarData vd = (*mit).second;
        if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), jCDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), iCDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), iDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else
          vsetdatas.push_back(vd);
      }
    }
    else {
      for (unsigned int i = 0; i < var_names.size(); i++) {
        mit = varInfo.find(var_names[i]);
        if (mit != varInfo.end()) {
          VarData vd = (*mit).second;
          if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), jCDim) != vd.varDims.end()))
            vdatas.push_back(vd);
          else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), iCDim) != vd.varDims.end()))
            vdatas.push_back(vd);
          else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), iDim) != vd.varDims.end()))
            vdatas.push_back(vd);
          else
            vsetdatas.push_back(vd);
        }
        else ERRORR(MB_FAILURE, "Couldn't find variable.");
      }
    }
  }
  else { // ucd mesh
    if (var_names.empty()) {
      for (mit = varInfo.begin(); mit != varInfo.end(); mit++) {
        VarData vd = (*mit).second;
        if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), kDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), iDim)
            != vd.varDims.end()))
          vdatas.push_back(vd); //3d data (time, ncol, ilev) read here
        else
          vsetdatas.push_back(vd);
      }
    }
    else {
      for (unsigned int i = 0; i < var_names.size(); i++) {

        mit = varInfo.find(var_names[i]);
        if (mit != varInfo.end()) {
          VarData vd = (*mit).second;
          if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), kDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), iDim)
              != vd.varDims.end()))
            vdatas.push_back(vd); //3d data
          else
            vsetdatas.push_back(vd);
        }
        else ERRORR(MB_FAILURE, "Couldn't find variable.");
      }
    }
  }

  if (tstep_nums.empty() && -1 != tMin) {
    // no timesteps input, get them all
    for (int i = tMin; i <= tMax; i++)
      tstep_nums.push_back(i);
  }
  if (!tstep_nums.empty()) {
    for (unsigned int i = 0; i < vdatas.size(); i++) {
      vdatas[i].varTags.resize(tstep_nums.size(), 0);
      vdatas[i].varDatas.resize(tstep_nums.size());
      vdatas[i].readDims.resize(tstep_nums.size());
      vdatas[i].readCounts.resize(tstep_nums.size());
    }
    for (unsigned int i = 0; i < vsetdatas.size(); i++) {
      if ((std::find(vsetdatas[i].varDims.begin(), vsetdatas[i].varDims.end(), tDim) != vsetdatas[i].varDims.end())
          && (vsetdatas[i].varDims.size() != 1)) {
        vsetdatas[i].varTags.resize(tstep_nums.size(), 0);
        vsetdatas[i].varDatas.resize(tstep_nums.size());
        vsetdatas[i].readDims.resize(tstep_nums.size());
        vsetdatas[i].readCounts.resize(tstep_nums.size());
      }
      else {
        vsetdatas[i].varTags.resize(1, 0);
        vsetdatas[i].varDatas.resize(1);
        vsetdatas[i].readDims.resize(1);
        vsetdatas[i].readCounts.resize(1);
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variable_allocate(EntityHandle file_set, std::vector<VarData> &vdatas, std::vector<int> &tstep_nums, bool scdMesh) {
  ErrorCode rval = MB_SUCCESS;

  std::vector<EntityHandle>* ehandles = NULL;
  Range* range = NULL;
  //std::vector<EntityHandle> verts_handles;
  std::vector<EntityHandle> ns_edges_handles;
  std::vector<EntityHandle> ew_edges_handles;
  //std::vector<EntityHandle> quads_handles;

  // get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
      verts.psize() == 1);
  //verts_handles.resize(verts.size());
  //std::copy(verts.begin(), verts.end(), verts_handles.begin());

  Range edges;
  rval = mbImpl->get_entities_by_dimension(file_set, 1, edges);
  ERRORR(rval, "Trouble getting edges in set.");
  ns_edges_handles.resize(edges.size() / 2);
  ew_edges_handles.resize(edges.size() / 2);

  // get north/south edges
  // FIXME: initialize ns_edges_handles to get the right order
  //std::copy(edges.begin(), edges.end(), ns_edges_handles.begin());

  // get east/west edges in set
  // FIXME: initialize ew_edges_handles to get the right order
  //std::copy(edges.begin(), edges.end(), ew_edges_handles.begin());

  // get quads in set
  Range quads;
  rval = mbImpl->get_entities_by_dimension(file_set, 2, quads);
  ERRORR(rval, "Trouble getting quads in set.");
  assert("Should only have a single quad subrange, since they were read in one shot" &&
      quads.psize() == 1);
  //quads_handles.resize(quads.size());
  //std::copy(quads.begin(), quads.end(), quads_handles.begin());

#ifdef USE_MPI
  moab::Range quads_owned;
  if (isParallel)
  {
    rval = myPcomm->filter_pstatus(quads, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1,
      &quads_owned);
    ERRORR(rval, "Trouble getting owned quads in set.");
  }
  else
    quads_owned=quads;// not running in parallel, but still with MPI
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      std::vector<std::string>::iterator vit;
      int idx_lev = 0;
      int idx_ilev = 0;
      if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end())
        idx_lev = vit - dimNames.begin();
      if ((vit = std::find(dimNames.begin(), dimNames.end(), "ilev")) != dimNames.end())
        idx_ilev = vit - dimNames.begin();
      if (std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), idx_lev) != vdatas[i].varDims.end())
        vdatas[i].numLev = dimVals[idx_lev];
      else if (std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), idx_ilev) != vdatas[i].varDims.end())
        vdatas[i].numLev = dimVals[idx_ilev];

      // get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag(vdatas[i], tstep_nums[t], vdatas[i].varTags[t], vdatas[i].numLev);
        ERRORR(rval, "Trouble getting tag.");
      }

      // assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int) t) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");
      }
      else if (vdatas[i].varDims[0] != tDim) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Non-default timestep number given for time-independent variable.");
      }

      // set up the dimensions and counts
      // first time
      vdatas[i].readDims[t].push_back(tstep_nums[t]);
      vdatas[i].readCounts[t].push_back(1);

      // then z/y/x
      if (vdatas[i].numLev != 1) {
        vdatas[i].readDims[t].push_back(0);
        vdatas[i].readCounts[t].push_back(vdatas[i].numLev);
      }

      switch (vdatas[i].entLoc) {
        case 0:
          // vertices
          if (scdMesh) {
            // only structured mesh has j parameter that multiplies i to get total # vertices
            vdatas[i].readDims[t].push_back(lDims[1]);
            vdatas[i].readCounts[t].push_back(lDims[4] - lDims[1] + 1);
            vdatas[i].readDims[t].push_back(lDims[0]);
            vdatas[i].readCounts[t].push_back(lDims[3] - lDims[0] + 1);
          }
          else
          {
            // we will start from the first localGid, actually; we will reset that
            // later on, anyway, in a loop
            vdatas[i].readDims[t].push_back(localGid[0]-1);
            vdatas[i].readCounts[t].push_back(localGid.size());
          }

          assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
          range = &verts;
          break;
        case 1:
          // north/south edges
          /*
           vdatas[i].readDims[t].push_back(lDims[1]);
           vdatas[i].readDims[t].push_back(lCDims[0]);
           // FIXME: need to figure out what to do with extra 2 pole points
           //vdatas[i].readCounts[t].push_back(lDims[4]-lDims[1]+1);
           vdatas[i].readCounts[t].push_back(lDims[4]-lDims[1]+1-2);
           vdatas[i].readCounts[t].push_back(lCDims[3]-lCDims[0]+1);
           assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
           ehandles = &ns_edges_handles;
           range = &verts; // FIXME: should remove when edge handles are used
           */
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.")
          ;
          break;
        case 2:
          // east/west edges
          /*
           vdatas[i].readDims[t].push_back(lCDims[1]);
           vdatas[i].readDims[t].push_back(lDims[0]);
           #ifdef USE_MPI
           if ((isParallel) &&
           (myPcomm->proc_config().proc_size() != 1) &&
           ((gDims[3]-gDims[0]) == (lDims[3]-lDims[0])))
           vdatas[i].readCounts[t].push_back(lDims[3]-lDims[0]);
           else
           #endif
           vdatas[i].readCounts[t].push_back(lCDims[4]-lCDims[1]+1);
           vdatas[i].readCounts[t].push_back(lDims[3]-lDims[0]);
           assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
           ehandles = &ew_edges_handles;
           range = &verts; // FIXME: should remove when edge handles are used
           */
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.")
          ;
          break;
        case 3:
          // quads
          vdatas[i].readDims[t].push_back(lCDims[1]);
          vdatas[i].readDims[t].push_back(lCDims[0]);
          vdatas[i].readCounts[t].push_back(lCDims[4] - lCDims[1] + 1);
          vdatas[i].readCounts[t].push_back(lCDims[3] - lCDims[0] + 1);
          assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());

#ifdef USE_MPI
          range = &quads_owned;
#else
          range = &quads;
#endif
          break;
        case 4:
          // set
          break;
        default:
          ERRORR(MB_FAILURE, "Unrecoganized entity location type.")
          ;
          break;
      }

      if (!ehandles) {
        // get ptr to tag space
        void *data;
        int count;
        rval = mbImpl->tag_iterate(vdatas[i].varTags[t], range->begin(), range->end(), count, data);
        ERRORR(rval, "Failed to get tag iterator.");
        assert((unsigned)count == range->size());
        vdatas[i].varDatas[t] = data;
      }
      else {
        //FIXME should use the correct vector of edge handles
        void *data;
        int count;
        rval = mbImpl->tag_iterate(vdatas[i].varTags[t], range->begin(), range->end(), count, data);
        ERRORR(rval, "Failed to get tag iterator.");
        assert((unsigned)count == range->size());
        vdatas[i].varDatas[t] = data;
      }
    }
  }

  return rval;
}

ErrorCode ReadNC::read_variables(EntityHandle file_set, std::vector<std::string> &var_names, std::vector<int> &tstep_nums, bool scdMesh) {
  std::vector<VarData> vdatas;
  std::vector<VarData> vsetdatas;

  ErrorCode rval = read_variable_setup(var_names, tstep_nums, vdatas, vsetdatas, scdMesh);
  ERRORR(rval, "Trouble setting up read variable.");

  if (scdMesh) {
    // create COORDS tag for quads
    rval = create_quad_coordinate_tag(file_set);
    ERRORR(rval, "Trouble creating coordinate tags to entities quads");
  }

  if (!vsetdatas.empty()) {
    rval = read_variable_to_set(file_set, vsetdatas, tstep_nums, scdMesh);
    ERRORR(rval, "Trouble read variables to set.");
  }

  if (!vdatas.empty()) {
#ifdef PNETCDF_FILE
    if (!scdMesh) // in serial, we will use the old read, everything is contiguous
      // in parallel, we will use async read in pnetcdf
      // the other mechanism is not working, forget about it
      rval = read_variable_to_nonset_async(file_set, vdatas, tstep_nums);
    else
#endif
      rval = read_variable_to_nonset(file_set, vdatas, tstep_nums, scdMesh);

    ERRORR(rval, "Trouble read variables to entities verts/edges/quads.");
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variable_to_set_allocate(std::vector<VarData> &vdatas, std::vector<int> &tstep_nums) {
  ErrorCode rval = MB_SUCCESS;

  for (unsigned int i = 0; i < vdatas.size(); ++i) {
    if ((std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), tDim) != vdatas[i].varDims.end()))
      vdatas[i].has_t = true;

    for (unsigned int t = 0; t < tstep_nums.size(); ++t) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      // get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag_to_set(vdatas[i], tstep_nums[t], vdatas[i].varTags[t]);
        ERRORR(rval, "Trouble getting tag.");
      }

      // assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int) t)
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");

      // set up the dimensions and counts
      // first variable dimension is time, if it exists
      if (vdatas[i].has_t)
      {
        if (vdatas[i].varDims.size() != 1)
        {
          vdatas[i].readDims[t].push_back(tstep_nums[t]);
          vdatas[i].readCounts[t].push_back(1);
        }
        else
        {
          vdatas[i].readDims[t].push_back(0);
          vdatas[i].readCounts[t].push_back(tstep_nums.size());
        }
      }

      // set up other dimensions and counts
      if (vdatas[i].varDims.empty()) {
        // scalar variable
        vdatas[i].readDims[t].push_back(0);
        vdatas[i].readCounts[t].push_back(1);
      }
      else {
        for (unsigned int idx = 0; idx != vdatas[i].varDims.size(); ++idx){
          if (tDim != vdatas[i].varDims[idx]){
            // push other variable dimensions, except time, which was already pushed
            vdatas[i].readDims[t].push_back(0);
            vdatas[i].readCounts[t].push_back(dimVals[vdatas[i].varDims[idx]]);
          }
        }
      }
      std::size_t sz = 1;
      for (std::size_t idx = 0; idx != vdatas[i].readCounts[t].size(); ++idx)
        sz *= vdatas[i].readCounts[t][idx];
      vdatas[i].sz = sz;
      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          vdatas[i].varDatas[t] = new char[sz];
          break;
        case NC_DOUBLE:
        case NC_FLOAT:
          vdatas[i].varDatas[t] = new double[sz];
          break;
        case NC_INT:
        case NC_SHORT:
          vdatas[i].varDatas[t] = new int[sz];
          break;
        default:
          std::cerr << "Unrecognized data type for tag " << std::endl;
          rval = MB_FAILURE;
      }
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

  return rval;
}

ErrorCode ReadNC::read_variable_to_set(EntityHandle file_set, std::vector<VarData> &vdatas, std::vector<int> &tstep_nums, bool scdMesh) {
  ErrorCode rval = read_variable_to_set_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables to set.");

  // finally, read into that space
  int success;
  std::vector<int> requests(vdatas.size() * tstep_nums.size()), statuss(vdatas.size() * tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if (dummyVarNames.find(vdatas[i].varName)!= dummyVarNames.end() )
       continue;// this is a dummy one, we don't have it; we created it for the dummy tag
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void *data = vdatas[i].varDatas[t];

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          success = NCFUNCAG(_vara_text)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (char*) data NCREQ);
          ERRORS(success, "Failed to read char data.")
          ;
          break;
        case NC_DOUBLE:
          success = NCFUNCAG(_vara_double)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (double*) data NCREQ);
          ERRORS(success, "Failed to read double data.")
          ;
          break;
        case NC_FLOAT: {
          success = NCFUNCAG(_vara_float)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (float*) data NCREQ);
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT:
          success = NCFUNCAG(_vara_int)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (int*) data NCREQ);
          ERRORS(success, "Failed to read int data.")
          ;
          break;
        case NC_SHORT:
          success = NCFUNCAG(_vara_short)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (short*) data NCREQ);
          ERRORS(success, "Failed to read short data.")
          ;
          break;
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t, scdMesh);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }
  // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Setting data for variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = mbImpl->tag_set_by_ptr(vdatas[i].varTags[t], &file_set, 1, &(vdatas[i].varDatas[t]), &vdatas[i].sz);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

  return rval;
}

ErrorCode ReadNC::read_variable_to_nonset(EntityHandle file_set, std::vector<VarData> &vdatas, std::vector<int> &tstep_nums, bool scdMesh) {
  ErrorCode rval = read_variable_allocate(file_set, vdatas, tstep_nums, scdMesh);
  ERRORR(rval, "Trouble allocating read variables.");

  // finally, read into that space
  int success;
  std::vector<int> requests(vdatas.size() * tstep_nums.size()), statuss(vdatas.size() * tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void *data = vdatas[i].varDatas[t];
      std::size_t sz = 1;
      size_t ni = vdatas[i].readCounts[t][2], nj = vdatas[i].readCounts[t][3], nk = vdatas[i].readCounts[t][1];
      if (scdMesh) {
        for (std::size_t idx = 0; idx != vdatas[i].readCounts[t].size(); ++idx)
          sz *= vdatas[i].readCounts[t][idx];
      }
      else {
        sz = vdatas[i].numLev * vdatas[i].readCounts[t][2];
        nj = 1; // for ucdMesh, nj holds # quads, so here should reset to 1
      }

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          std::vector<char> tmpchardata(sz);
          success = NCFUNCAG(_vara_text)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpchardata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpchardata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpchardata.size(); ++idx)
              ((char*) data)[idx] = tmpchardata[idx];
          }
          ERRORS(success, "Failed to read char data.");
          break;
        }
        case NC_DOUBLE: {
          std::vector<double> tmpdoubledata(sz);
          success = NCFUNCAG(_vara_double)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpdoubledata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpdoubledata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); ++idx)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);

          if (scdMesh)
          {
            success = NCFUNCAG(_vara_float)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
                &tmpfloatdata[0] NCREQ);
          }
          else
          {
            // in the case of ucd mesh, and on multiple proc,
            // we need to read as many times as subranges we have in the
            // localGid range;
            // basically, we have to give a different point
            // for data to start, for every subrange :(
            size_t nbDims=vdatas[i].readDims[t].size();
            // assume that the last dimension is for the ncol,
            // node varying variable

            size_t indexInFloatArray= 0;
            size_t ic=0;
            for (
                Range::pair_iterator pair_iter = localGid.pair_begin();
                pair_iter!=localGid.pair_end();
                pair_iter++, ic++)
            {
              EntityHandle starth = pair_iter->first;
              EntityHandle endh = pair_iter->second;// inclusive
              vdatas[i].readDims[t][nbDims-1] = (NCDF_SIZE) (starth-1);
              vdatas[i].readCounts[t][nbDims-1] = (NCDF_SIZE) (endh-starth+1);

              success = NCFUNCAG(_vara_float)(fileId, vdatas[i].varId,
                  &(vdatas[i].readDims[t][0]), &(vdatas[i].readCounts[t][0]),
                              &(tmpfloatdata[indexInFloatArray]) NCREQ);
              ERRORS(success, "Failed to read float data in loop");
              // we need to increment the index in float array for the
              // next subrange
              indexInFloatArray+= (endh-starth+1)*1*vdatas[i].numLev; //
            }
            assert(ic==localGid.psize());
            //
          }
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpfloatdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); ++idx)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT: {
          std::vector<int> tmpintdata(sz);
          success = NCFUNCAG(_vara_int)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpintdata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpintdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpintdata.size(); ++idx)
              ((int*) data)[idx] = tmpintdata[idx];
          }
          ERRORS(success, "Failed to read int data.");
          break;
        }
        case NC_SHORT: {
          std::vector<short> tmpshortdata(sz);
          success = NCFUNCAG(_vara_short)(fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpshortdata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpshortdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpshortdata.size(); ++idx)
              ((short*) data)[idx] = tmpshortdata[idx];
          }
          ERRORS(success, "Failed to read short data.");
          break;
        }
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");
    }
  }

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t, scdMesh);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
    }
  }
  // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  return rval;
}

#ifdef PNETCDF_FILE
// here we know all variables are on vertices, ncol dimension, and nlevel
ErrorCode ReadNC::read_variable_to_nonset_async(EntityHandle file_set, std::vector<VarData> &vdatas,
              std::vector<int> &tstep_nums)
{
  ErrorCode rval = read_variable_allocate(file_set, vdatas, tstep_nums, false);
  ERRORR(rval, "Trouble allocating read variables.");

  // finally, read into that space
  int success;
  // MPI_offset or size_t?
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      // we will synchronize all these reads with the other processors,
      // so the wait will be inside this double loop; is it too much?
      size_t nb_reads = localGid.psize();
      std::vector<int> requests(nb_reads), statuss(nb_reads);
      size_t idxReq=0;
      void *data = vdatas[i].varDatas[t];
      std::size_t sz = 1;
      size_t ni = vdatas[i].readCounts[t][2];
      //size_t nj = vdatas[i].readCounts[t][3];
      size_t nk = vdatas[i].readCounts[t][1];

      sz = vdatas[i].numLev * vdatas[i].readCounts[t][2];
      size_t nj = 1; // for ucdMesh, nj holds # quads, so here should reset to 1

      if (sz<=0)
        continue;// nothing to read, why worry?

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_DOUBLE: {
          // copy from float case
          std::vector<double> tmpdoubledata(sz);

          // in the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGid range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims=vdatas[i].readDims[t].size();
          // assume that the last dimension is for the ncol,
          // node varying variable

          size_t indexInDoubleArray= 0;
          size_t ic=0;
          for (
              Range::pair_iterator pair_iter = localGid.pair_begin();
              pair_iter!=localGid.pair_end();
              pair_iter++, ic++)
          {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second;// inclusive
            vdatas[i].readDims[t][nbDims-1] = (NCDF_SIZE) (starth-1);
            vdatas[i].readCounts[t][nbDims-1] = (NCDF_SIZE) (endh-starth+1);

            // do a partial read, in each subrange
            // wait outside this loop
            success = NCFUNCAG2(_vara_double)(fileId, vdatas[i].varId,
                &(vdatas[i].readDims[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpdoubledata[indexInDoubleArray]) NCREQ2);
            ERRORS(success, "Failed to read double data in loop");
            // we need to increment the index in float array for the
            // next subrange
            indexInDoubleArray+= (endh-starth+1)*1*vdatas[i].numLev; //
          }
          assert(ic==localGid.psize());
          //
          success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpdoubledata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); ++idx)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);

          // in the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGid range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims=vdatas[i].readDims[t].size();
          // assume that the last dimension is for the ncol,
          // node varying variable

          size_t indexInFloatArray= 0;
          size_t ic=0;
          for (
              Range::pair_iterator pair_iter = localGid.pair_begin();
              pair_iter!=localGid.pair_end();
              pair_iter++, ic++)
          {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second;// inclusive
            vdatas[i].readDims[t][nbDims-1] = (NCDF_SIZE) (starth-1);
            vdatas[i].readCounts[t][nbDims-1] = (NCDF_SIZE) (endh-starth+1);

            // do a partial read, in each subrange
            // wait outside this loop
            success = NCFUNCAG2(_vara_float)(fileId, vdatas[i].varId,
                &(vdatas[i].readDims[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpfloatdata[indexInFloatArray]) NCREQ2);
            ERRORS(success, "Failed to read float data in loop");
            // we need to increment the index in float array for the
            // next subrange
            indexInFloatArray+= (endh-starth+1)*1*vdatas[i].numLev; //
          }
          assert(ic==localGid.psize());
          //
          success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpfloatdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); ++idx)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_SHORT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");
    }
  }

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t, false);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
    }
  }
  // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  return rval;
}
#endif

ErrorCode ReadNC::convert_variable(VarData &var_data, int tstep_num, bool scdMesh) {
  // get ptr to tag space
  void *data = var_data.varDatas[tstep_num];

  std::size_t sz = 1;
  if (scdMesh) {
    for (std::size_t idx = 0; idx != var_data.readCounts[tstep_num].size(); ++idx)
      sz *= var_data.readCounts[tstep_num][idx];
  }
  else
    sz = var_data.numLev * localGid.size();// how many nodes are we reading?

  // finally, read into that space
  int success = 0, *idata;
  double *ddata;
  float *fdata;
  short *sdata;

  switch (var_data.varDataType) {
    case NC_FLOAT:
      ddata = (double*) var_data.varDatas[tstep_num];
      fdata = (float*) var_data.varDatas[tstep_num];
      // convert in-place
      for (int i = sz - 1; i >= 0; i--)
        ddata[i] = fdata[i];
      break;
    case NC_SHORT:
      idata = (int*) var_data.varDatas[tstep_num];
      sdata = (short*) var_data.varDatas[tstep_num];
      // convert in-place
      for (int i = sz - 1; i >= 0; i--)
        idata[i] = sdata[i];
      break;
    default:
      success = 1;
  }

  if (2 <= dbgOut.get_verbosity() && !success) {
    double dmin, dmax;
    int imin, imax;
    switch (var_data.varDataType) {
      case NC_DOUBLE:
      case NC_FLOAT:
        ddata = (double*) data;
        if (sz == 0)
          break;

        dmin = dmax = ddata[0];
        for (unsigned int i = 1; i < sz; i++) {
          if (ddata[i] < dmin)
            dmin = ddata[i];
          if (ddata[i] > dmax)
            dmax = ddata[i];
        }
        dbgOut.tprintf(2, "Variable %s (double): min = %f, max = %f\n", var_data.varName.c_str(), dmin, dmax);
        break;
      case NC_INT:
      case NC_SHORT:
        idata = (int*) data;
        if (sz == 0)
          break;

        imin = imax = idata[0];
        for (unsigned int i = 1; i < sz; i++) {
          if (idata[i] < imin)
            imin = idata[i];
          if (idata[i] > imax)
            imax = idata[i];
        }
        dbgOut.tprintf(2, "Variable %s (int): min = %d, max = %d\n", var_data.varName.c_str(), imin, imax);
        break;
      case NC_NAT:
      case NC_BYTE:
      case NC_CHAR:
        break;
      default: //default case added to remove compiler warnings
        success = 1;
    }
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_tag_to_set(VarData &var_data, int tstep_num, Tag &tagh) {
  std::ostringstream tag_name;
  if ((!var_data.has_t)||( var_data.varDims.size()<=1))
    tag_name << var_data.varName;
  else if (!tstep_num) {
    std::string tmp_name = var_data.varName + "0";
    tag_name << tmp_name.c_str();
  }
  else
    tag_name << var_data.varName << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_OPAQUE, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    case NC_DOUBLE:
    case NC_FLOAT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_DOUBLE, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    case NC_INT:
    case NC_SHORT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_INTEGER, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    default:
      std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
      rval = MB_FAILURE;
  }

  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

ErrorCode ReadNC::get_tag(VarData &var_data, int tstep_num, Tag &tagh, int num_lev) {
  std::ostringstream tag_name;
  if (!tstep_num) {
    std::string tmp_name = var_data.varName + "0";
    tag_name << tmp_name.c_str();
  }
  else
    tag_name << var_data.varName << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_OPAQUE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    case NC_DOUBLE:
    case NC_FLOAT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    case NC_INT:
    case NC_SHORT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_INTEGER, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    default:
      std::cerr << "Unrecognized data type for tag " << tag_name.str() << std::endl;
      rval = MB_FAILURE;
  }

  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

void ReadNC::init_dims_with_no_cvars_info() {
  // hack: look at all dimensions, and see if we have one that does not appear in the list of varInfo names
  // right now, candidates are ncol and nbnd
  // for them, create dummy tags
  for (unsigned int i=0; i<dimNames.size(); i++)
  {
    // if there is a var with this name, skip, we are fine; if not, create a varInfo...
    if ( varInfo.find(dimNames[i])!=varInfo.end())
      continue; // we already have a variable with this dimension name

    int sizeTotalVar = varInfo.size();
    std::string var_name(dimNames[i]);
    VarData &data = varInfo[var_name];
    data.varName = std::string(var_name);
    data.varId =sizeTotalVar;
    data.varTags.resize(1, 0);
    data.varDataType = NC_DOUBLE; // could be int, actually, but we do not really need the type
    data.varDims.resize(1);
    data.varDims[0]= (int)i;
    data.numAtts=0;
    data.entLoc = ENTLOCSET;
    dbgOut.tprintf(2, "Dummy varInfo created for dimension %s\n", dimNames[i].c_str());
    dummyVarNames.insert(dimNames[i]);
  }
}

ErrorCode ReadNC::read_coordinate(const char *var_name, int lmin, int lmax, std::vector<double> &cvals) {
  std::map<std::string, VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit)
    return MB_FAILURE;

  // check to make sure it's a float or double
  int fail;
  NCDF_SIZE tmin = lmin, tcount = lmax - lmin + 1;
  NCDF_DIFF dum_stride = 1;
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    cvals.resize(tcount);
    fail = NCFUNCA(get_vars_double)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &cvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    fail = NCFUNCA(get_vars_float)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &tcvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
    std::copy(tcvals.begin(), tcvals.end(), cvals.begin());
  }
  else ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_header() {
  CPU_WORD_SIZE = sizeof(double);
  IO_WORD_SIZE = sizeof(double);

  dbgOut.tprint(1, "Reading header...\n");

  // get the global attributes
  int numgatts;
  int success;
  success = NCFUNC(inq_natts )(fileId, &numgatts);
  ERRORS(success, "Couldn't get number of global attributes.");

  // read attributes into globalAtts
  ErrorCode result = get_attributes(NC_GLOBAL, numgatts, globalAtts);
  ERRORR(result, "Getting attributes.");
  dbgOut.tprintf(1, "Read %u attributes\n", (unsigned int) globalAtts.size());

  // read in dimensions into dimVals
  result = get_dimensions(fileId, dimNames, dimVals);
  ERRORR(result, "Getting dimensions.");
  dbgOut.tprintf(1, "Read %u dimensions\n", (unsigned int) dimVals.size());

  // read in variables into varInfo
  result = get_variables();
  ERRORR(result, "Getting variables.");
  dbgOut.tprintf(1, "Read %u variables\n", (unsigned int) varInfo.size());

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_attributes(int var_id, int num_atts, std::map<std::string, AttData> &atts, const char *prefix) {

  char dum_name[120];

  for (int i = 0; i < num_atts; i++) {
    // get the name
    int success = NCFUNC(inq_attname)(fileId, var_id, i, dum_name);
    ERRORS(success, "Trouble getting attribute name.");

    AttData &data = atts[std::string(dum_name)];
    data.attName = std::string(dum_name);
    success = NCFUNC(inq_att)(fileId, var_id, dum_name, &data.attDataType, &data.attLen);
    ERRORS(success, "Trouble getting attribute info.");
    data.attVarId = var_id;

    dbgOut.tprintf(2, "%sAttribute %s: length=%u, varId=%d, type=%d\n", (prefix ? prefix : ""), data.attName.c_str(),
        (unsigned int) data.attLen, data.attVarId, data.attDataType);
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_dimensions(int file_id, std::vector<std::string> &dim_names, std::vector<int> &dim_vals) {
  // get the number of dimensions
  int num_dims;
  int success = NCFUNC(inq_ndims)(file_id, &num_dims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (num_dims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n", num_dims, (int) NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME + 1];
  NCDF_SIZE dum_len;
  dim_names.resize(num_dims);
  dim_vals.resize(num_dims);

  for (int i = 0; i < num_dims; i++) {
    success = NCFUNC(inq_dim)(file_id, i, dim_name, &dum_len);
    ERRORS(success, "Trouble getting dimension info.");

    dim_vals[i] = dum_len;
    dim_names[i] = std::string(dim_name);

    dbgOut.tprintf(2, "Dimension %s, length=%u\n", dim_name, (unsigned int) dum_len);
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_variables() {
  // first cache the number of time steps
  std::vector<std::string>::iterator vit = std::find(dimNames.begin(), dimNames.end(), "time");
  if (vit == dimNames.end())
    vit = std::find(dimNames.begin(), dimNames.end(), "t");

  int ntimes = 0;
  if (vit != dimNames.end())
    ntimes = dimVals[vit - dimNames.begin()];
  if (!ntimes)
    ntimes = 1;

  // get the number of variables
  int num_vars;
  int success = NCFUNC(inq_nvars)(fileId, &num_vars);
  ERRORS(success, "Trouble getting number of variables.");

  if (num_vars > NC_MAX_VARS) {
    readMeshIface->report_error("ReadNC: File contains %d vars but NetCDF library supports only %d\n", num_vars, (int) NC_MAX_VARS);
    return MB_FAILURE;
  }

  char var_name[NC_MAX_NAME + 1];
  int var_ndims;

  for (int i = 0; i < num_vars; i++) {
    // get the name first, so we can allocate a map iterate for this var
    success = NCFUNC(inq_varname )(fileId, i, var_name);
    ERRORS(success, "Trouble getting var name.");
    VarData &data = varInfo[std::string(var_name)];
    data.varName = std::string(var_name);
    data.varId = i;
    data.varTags.resize(ntimes, 0);

    // get the data type
    success = NCFUNC(inq_vartype)(fileId, i, &data.varDataType);
    ERRORS(success, "Trouble getting variable data type.");

    // get the number of dimensions, then the dimensions
    success = NCFUNC(inq_varndims)(fileId, i, &var_ndims);
    ERRORS(success, "Trouble getting number of dims of a variable.");
    data.varDims.resize(var_ndims);

    success = NCFUNC(inq_vardimid)(fileId, i, &data.varDims[0]);
    ERRORS(success, "Trouble getting variable dimensions.");

    // finally, get the number of attributes, then the attributes
    success = NCFUNC(inq_varnatts)(fileId, i, &data.numAtts);
    ERRORS(success, "Trouble getting number of dims of a variable.");

    // print debug info here so attribute info comes afterwards
    dbgOut.tprintf(2, "Variable %s: Id=%d, numAtts=%d, datatype=%d, num_dims=%u\n", data.varName.c_str(), data.varId, data.numAtts,
        data.varDataType, (unsigned int) data.varDims.size());

    ErrorCode rval = get_attributes(i, data.numAtts, data.varAtts, "   ");
    ERRORR(rval, "Trouble getting attributes for a variable.");

  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_tag_values(const char*, const char*, const FileOptions&, std::vector<int>&, const SubsetList*) {
  return MB_FAILURE;
}

ErrorCode ReadNC::create_tags(ScdInterface *scdi, EntityHandle file_set, const std::vector<int>& tstep_nums) {
  ErrorCode rval;
  std::string tag_name;

  // <__NUM_DIMS>
  Tag numDimsTag = 0;
  tag_name = "__NUM_DIMS";
  int numDims = dimNames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_DIMS tag.");
  rval = mbImpl->tag_set_data(numDimsTag, &file_set, 1, &numDims);
  ERRORR(rval, "Trouble setting data for __NUM_DIMS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__NUM_VARS>
  Tag numVarsTag = 0;
  tag_name = "__NUM_VARS";
  int numVars = varInfo.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numVarsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_VARS tag.");
  rval = mbImpl->tag_set_data(numVarsTag, &file_set, 1, &numVars);
  ERRORR(rval, "Trouble setting data for __NUM_VARS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_NAMES>
  Tag dimNamesTag = 0;
  tag_name = "__DIM_NAMES";
  std::string dimnames;
  unsigned int dimNamesSz = dimNames.size();
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    dimnames.append(dimNames[i]);
    dimnames.push_back('\0');
  }
  int dimnamesSz = dimnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, dimNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_NAMES tag.");
  const void* ptr = dimnames.c_str();
  rval = mbImpl->tag_set_by_ptr(dimNamesTag, &file_set, 1, &ptr, &dimnamesSz);
  ERRORR(rval, "Trouble setting data for __DIM_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_VALUES>
  Tag dimValsTag = 0;
  tag_name = "__DIM_VALUES";
  //std::vector<int> dim;
  int dimValsSz = (int)dimVals.size();

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_INTEGER, dimValsTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_VALUES tag.");
  ptr = &(dimVals[0]);
  rval = mbImpl->tag_set_by_ptr(dimValsTag, &file_set, 1, &ptr, &dimValsSz);
  ERRORR(rval, "Trouble setting data for __DIM_VALUES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__VAR_NAMES>
  Tag varNamesTag = 0;
  tag_name = "__VAR_NAMES";
  std::string varnames;
  std::map<std::string, VarData>::iterator mapIter;
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varnames.append(mapIter->first);
    varnames.push_back('\0');
  }
  int varnamesSz = varnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __VAR_NAMES tag.");
  ptr = varnames.c_str();
  rval = mbImpl->tag_set_by_ptr(varNamesTag, &file_set, 1, &ptr, &varnamesSz);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // __<dim_name>_LOC_MINMAX
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] == "time") {
      std::stringstream ss_tag_name;
      ss_tag_name << "__" << dimNames[i] << "_LOC_MINMAX";
      tag_name = ss_tag_name.str();
      Tag tagh = 0;
      std::vector<int> val(2, 0);
      val[0] = tMin;
      val[1] = tMax;
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
      ERRORR(rval, "Trouble creating __<dim_name>_LOC_MINMAX tag.");
      rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
      ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_MINMAX tag.");
      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    }
  }

  // __<dim_name>_LOC_VALS
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] != "time")
      continue;
    std::vector<int> val;
    if (!tstep_nums.empty())
      val = tstep_nums;
    else {
      val.resize(tVals.size());
      for (unsigned int j = 0; j != tVals.size(); ++j)
        val[j] = j;
    }
    Tag tagh = 0;
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << dimNames[i] << "_LOC_VALS";
    tag_name = ss_tag_name.str();
    rval = mbImpl->tag_get_handle(tag_name.c_str(), val.size(), MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<dim_name>_LOC_VALS tag.");
    rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
    ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_VALS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // __<var_name>_DIMS
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    Tag varNamesDimsTag = 0;
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << mapIter->first << "_DIMS";
    tag_name = ss_tag_name.str();
    unsigned int varDimSz = varInfo[mapIter->first].varDims.size();
    if (varDimSz == 0)
      continue;
    varInfo[mapIter->first].varTags.resize(varDimSz, 0);
    for (unsigned int i = 0; i != varDimSz; ++i) {
      Tag tmptag = 0;
      std::string tmptagname = dimNames[varInfo[mapIter->first].varDims[i]];
      mbImpl->tag_get_handle(tmptagname.c_str(), 0, MB_TYPE_OPAQUE, tmptag, MB_TAG_ANY);
      varInfo[mapIter->first].varTags[i] = tmptag;
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varDimSz, MB_TYPE_HANDLE, varNamesDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_DIMS tag.");
    rval = mbImpl->tag_set_data(varNamesDimsTag, &file_set, 1, &(varInfo[mapIter->first].varTags[0]));
    ERRORR(rval, "Trouble setting data for __<var_name>_DIMS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // <PARTITION_METHOD>
  Tag part_tag = scdi->part_method_tag();
  if (!part_tag)
    ERRORR(MB_FAILURE, "Trouble getting partition method tag.");
  rval = mbImpl->tag_set_data(part_tag, &file_set, 1, &partMethod);
  ERRORR(rval, "Trouble setting data for PARTITION_METHOD tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__GLOBAL_ATTRIBS>
  tag_name = "__GLOBAL_ATTRIBS";
  Tag globalAttTag = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, globalAttTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS tag.");
  std::string gattVal;
  std::vector<int> gattLen;
  rval = create_attrib_string(globalAtts, gattVal, gattLen);
  ERRORR(rval, "Trouble creating attribute strings.");
  const void* gattptr = gattVal.c_str();
  int globalAttSz = gattVal.size();
  rval = mbImpl->tag_set_by_ptr(globalAttTag, &file_set, 1, &gattptr, &globalAttSz);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__GLOBAL_ATTRIBS_LEN>
  tag_name = "__GLOBAL_ATTRIBS_LEN";
  Tag globalAttLenTag = 0;
  if (gattLen.size() == 0)
    gattLen.push_back(0);
  rval = mbImpl->tag_get_handle(tag_name.c_str(), gattLen.size(), MB_TYPE_INTEGER, globalAttLenTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS_LEN tag.");
  rval = mbImpl->tag_set_data(globalAttLenTag, &file_set, 1, &gattLen[0]);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS_LEN tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // __<var_name>_ATTRIBS and __<var_name>_ATTRIBS_LEN
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    std::stringstream ssTagName;
    ssTagName << "__" << mapIter->first << "_ATTRIBS";
    tag_name = ssTagName.str();
    Tag varAttTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varAttTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS tag.");
    std::string varAttVal;
    std::vector<int> varAttLen;
    rval = create_attrib_string(mapIter->second.varAtts, varAttVal, varAttLen);
    ERRORR(rval, "Trouble creating attribute strings.");
    const void* varAttPtr = varAttVal.c_str();
    int varAttSz = varAttVal.size();
    rval = mbImpl->tag_set_by_ptr(varAttTag, &file_set, 1, &varAttPtr, &varAttSz);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    if (varAttLen.size() == 0)
      varAttLen.push_back(0);
    ssTagName << "_LEN";
    tag_name = ssTagName.str();
    Tag varAttLenTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varAttLen.size(), MB_TYPE_INTEGER, varAttLenTag, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS_LEN tag.");
    rval = mbImpl->tag_set_data(varAttLenTag, &file_set, 1, &varAttLen[0]);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS_LEN tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // <__VAR_NAMES_LOCATIONS>
  tag_name = "__VAR_NAMES_LOCATIONS";
  Tag varNamesLocsTag = 0;
  std::vector<int> varNamesLocs(varInfo.size());
  rval = mbImpl->tag_get_handle(tag_name.c_str(), varNamesLocs.size(), MB_TYPE_INTEGER, varNamesLocsTag, MB_TAG_CREAT
      | MB_TAG_SPARSE);
  ERRORR(rval, "Trouble creating __VAR_NAMES_LOCATIONS tag.");
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varNamesLocs[std::distance(varInfo.begin(), mapIter)] = mapIter->second.entLoc;
  }
  rval = mbImpl->tag_set_data(varNamesLocsTag, &file_set, 1, &varNamesLocs[0]);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES_LOCATIONS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__MESH_TYPE>
  Tag meshTypeTag = 0;
  tag_name = "__MESH_TYPE";
  std::string meshTypeName = helper->get_mesh_type_name();

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, meshTypeTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __MESH_TYPE tag.");
  ptr = meshTypeName.c_str();
  int leng= meshTypeName.size();
  rval = mbImpl->tag_set_by_ptr(meshTypeTag, &file_set, 1, &ptr, &leng);
  ERRORR(rval, "Trouble setting data for __MESH_TYPE tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  return MB_SUCCESS;
}

ErrorCode ReadNC::create_attrib_string(const std::map<std::string, AttData>& attMap, std::string& attVal, std::vector<int>& attLen) {
  int success;
  std::stringstream ssAtt;
  unsigned int sz = 0;
  std::map<std::string, AttData>::const_iterator attIt = attMap.begin();
  for (; attIt != attMap.end(); ++attIt) {
    ssAtt << attIt->second.attName;
    ssAtt << '\0';
    void* attData = NULL;
    switch (attIt->second.attDataType) {
      case NC_BYTE:
      case NC_CHAR:
        sz = attIt->second.attLen;
        attData = (char *) malloc(sz);
        success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (char*) attData);
        ERRORS(success, "Failed to read attribute char data.")
        ;
        ssAtt << "char;";
        break;
      case NC_DOUBLE:
        sz = attIt->second.attLen * sizeof(double);
        attData = (double *) malloc(sz);
        success = NCFUNC(get_att_double)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (double*) attData);
        ERRORS(success, "Failed to read attribute double data.")
        ;
        ssAtt << "double;";
        break;
      case NC_FLOAT:
        sz = attIt->second.attLen * sizeof(float);
        attData = (float *) malloc(sz);
        success = NCFUNC(get_att_float)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (float*) attData);
        ERRORS(success, "Failed to read attribute float data.")
        ;
        ssAtt << "float;";
        break;
      case NC_INT:
        sz = attIt->second.attLen * sizeof(int);
        attData = (int *) malloc(sz);
        success = NCFUNC(get_att_int)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (int*) attData);
        ERRORS(success, "Failed to read attribute int data.")
        ;
        ssAtt << "int;";
        break;
      case NC_SHORT:
        sz = attIt->second.attLen * sizeof(short);
        attData = (short *) malloc(sz);
        success = NCFUNC(get_att_short)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (short*) attData);
        ERRORS(success, "Failed to read attribute short data.")
        ;
        ssAtt << "short;";
        break;
      default:
        success = 1;
    }
    char* tmpc = (char *) attData;
    for (unsigned int counter = 0; counter != sz; ++counter)
      ssAtt << tmpc[counter];
    free(attData);
    ssAtt << ';';
    attLen.push_back(ssAtt.str().size() - 1);
  }
  attVal = ssAtt.str();

  return MB_SUCCESS;
}

ErrorCode ReadNC::create_quad_coordinate_tag(EntityHandle file_set) {
  Range ents;
  ErrorCode rval = mbImpl->get_entities_by_type(file_set, moab::MBQUAD, ents);
  ERRORR(rval, "Trouble getting QUAD entity.");

  std::size_t numOwnedEnts = 0;
#ifdef USE_MPI
  Range ents_owned;
  if (isParallel){
    rval = myPcomm->filter_pstatus(ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ents_owned);
    ERRORR(rval, "Trouble getting owned QUAD entity.");
    numOwnedEnts = ents_owned.size();
  }
  else
  {
    numOwnedEnts = ents.size();
    ents_owned = ents;
  }
#else
  numOwnedEnts = ents.size();
#endif

  if (numOwnedEnts == 0)
    return MB_SUCCESS;

  assert(numOwnedEnts == ilCVals.size() * jlCVals.size());
  std::vector<double> coords(numOwnedEnts * 3);
  std::size_t pos = 0;
  for (std::size_t j = 0; j != jlCVals.size(); ++j) {
    for (std::size_t i = 0; i != ilCVals.size(); ++i) {
      pos = j * ilCVals.size() * 3 + i * 3;
      coords[pos] = ilCVals[i];
      coords[pos + 1] = jlCVals[j];
      coords[pos + 2] = 0.0;
    }
  }
  std::string tag_name = "COORDS";
  Tag tagh = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating COORDS tag.");

  void *data;
  int count;
#ifdef USE_MPI
  rval = mbImpl->tag_iterate(tagh, ents_owned.begin(), ents_owned.end(), count, data);
#else
  rval = mbImpl->tag_iterate(tagh, ents.begin(), ents.end(), count, data);
#endif
  ERRORR(rval, "Failed to get COORDS tag iterator.");
  assert(count == (int)numOwnedEnts);
  double* quad_data = (double*) data;
  std::copy(coords.begin(), coords.end(), quad_data);
  return MB_SUCCESS;
}

} // namespace moab
