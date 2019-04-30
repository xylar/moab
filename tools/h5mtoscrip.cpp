// 
// Usage:
// tools/mbslavepart -d 2 -m mpas/x1.2562.grid.h5m -s mpas/x1.10242.grid.h5m -o mpas_slave.h5m -e 1e-8 -b 1e-6 -O
// 
#include <iostream>
#include <exception>
#include <cmath>
#include <vector>
#include <string>

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#endif

using namespace moab;

#if 0
void WriteNetCDF4(
  const std::string & scripfile,
  const std::map<std::string, std::string> & mapAttributes,
  NcFile::FileFormat eOutputFormat
) {
  // Temporarily change rval reporting
  NcError error_temp(NcError::verbose_fatal);

  // Open an output file
  NcFile ncMap(scripfile.c_str(), NcFile::Replace, NULL, 0, eOutputFormat);
  if (!ncMap.is_valid()) {
    _EXCEPTION1("Unable to open output map file \"%s\"",
      scripfile.c_str());
  }

  // Attributes
  ncMap.add_att("Title", "TempestRemap Offline Regridding Weight Generator");

  // Map dimensions
  int nA = (int)(m_dSourceAreas.GetRows());
  int nB = (int)(m_dTargetAreas.GetRows());

  // Write output dimensions entries
  int nSrcGridDims = (int)(m_vecSourceDimSizes.size());
  int nDstGridDims = (int)(m_vecTargetDimSizes.size());

  NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
  NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

  NcVar * varSrcGridDims =
    ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
  NcVar * varDstGridDims =
    ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

  char szDim[64];
  if ((nSrcGridDims == 1) && (m_vecSourceDimSizes[0] != nA)) {
    varSrcGridDims->put(&nA, 1);
    varSrcGridDims->add_att("name0", "num_dof");

  } else {
    for (int i = 0; i < m_vecSourceDimSizes.size(); i++) {
      varSrcGridDims->set_cur(nSrcGridDims - i - 1);
      varSrcGridDims->put(&(m_vecSourceDimSizes[i]), 1);
    }

    for (int i = 0; i < m_vecSourceDimSizes.size(); i++) {
      sprintf(szDim, "name%i", i);
      varSrcGridDims->add_att(szDim,
        m_vecSourceDimNames[nSrcGridDims - i - 1].c_str());
    }
  }

  if ((nDstGridDims == 1) && (m_vecTargetDimSizes[0] != nB)) {
    varDstGridDims->put(&nB, 1);
    varDstGridDims->add_att("name0", "num_dof");

  } else {
    for (int i = 0; i < m_vecTargetDimSizes.size(); i++) {
      varDstGridDims->set_cur(nDstGridDims - i - 1);
      varDstGridDims->put(&(m_vecTargetDimSizes[i]), 1);
    }

    for (int i = 0; i < m_vecTargetDimSizes.size(); i++) {
      sprintf(szDim, "name%i", i);
      varDstGridDims->add_att(szDim,
        m_vecTargetDimNames[nDstGridDims - i - 1].c_str());
    }
  }

  // Source and Target mesh resolutions
  NcDim * dimNA = ncMap.add_dim("n_a", nA);
  NcDim * dimNB = ncMap.add_dim("n_b", nB);

  // Number of nodes per Face
  int nSourceNodesPerFace = m_dSourceVertexLon.GetColumns();
  int nTargetNodesPerFace = m_dTargetVertexLon.GetColumns();

  NcDim * dimNVA = ncMap.add_dim("nv_a", nSourceNodesPerFace);
  NcDim * dimNVB = ncMap.add_dim("nv_b", nTargetNodesPerFace);

  // Write coordinates
  NcVar * varYCA = ncMap.add_var("yc_a", ncDouble, dimNA);
  NcVar * varYCB = ncMap.add_var("yc_b", ncDouble, dimNB);

  NcVar * varXCA = ncMap.add_var("xc_a", ncDouble, dimNA);
  NcVar * varXCB = ncMap.add_var("xc_b", ncDouble, dimNB);

  NcVar * varYVA = ncMap.add_var("yv_a", ncDouble, dimNA, dimNVA);
  NcVar * varYVB = ncMap.add_var("yv_b", ncDouble, dimNB, dimNVB);

  NcVar * varXVA = ncMap.add_var("xv_a", ncDouble, dimNA, dimNVA);
  NcVar * varXVB = ncMap.add_var("xv_b", ncDouble, dimNB, dimNVB);

  varYCA->add_att("units", "degrees");
  varYCB->add_att("units", "degrees");

  varXCA->add_att("units", "degrees");
  varXCB->add_att("units", "degrees");

  varYVA->add_att("units", "degrees");
  varYVB->add_att("units", "degrees");

  varXVA->add_att("units", "degrees");
  varXVB->add_att("units", "degrees");

  // Verify dimensionality
  if (m_dSourceCenterLon.GetRows() != nA) {
    _EXCEPTIONT("Mismatch between m_dSourceCenterLon and nA");
  }
  if (m_dSourceCenterLat.GetRows() != nA) {
    _EXCEPTIONT("Mismatch between m_dSourceCenterLat and nA");
  }
  if (m_dTargetCenterLon.GetRows() != nB) {
    _EXCEPTIONT("Mismatch between m_dTargetCenterLon and nB");
  }
  if (m_dTargetCenterLat.GetRows() != nB) {
    _EXCEPTIONT("Mismatch between m_dTargetCenterLat and nB");
  }
  if (m_dSourceVertexLon.GetRows() != nA) {
    _EXCEPTIONT("Mismatch between m_dSourceVertexLon and nA");
  }
  if (m_dSourceVertexLat.GetRows() != nA) {
    _EXCEPTIONT("Mismatch between m_dSourceVertexLat and nA");
  }
  if (m_dTargetVertexLon.GetRows() != nB) {
    _EXCEPTIONT("Mismatch between m_dTargetVertexLon and nB");
  }
  if (m_dTargetVertexLat.GetRows() != nB) {
    _EXCEPTIONT("Mismatch between m_dTargetVertexLat and nB");
  }

  varYCA->put(&(m_dSourceCenterLat[0]), nA);
  varYCB->put(&(m_dTargetCenterLat[0]), nB);

  varXCA->put(&(m_dSourceCenterLon[0]), nA);
  varXCB->put(&(m_dTargetCenterLon[0]), nB);

  varYVA->put(&(m_dSourceVertexLat[0][0]), nA, nSourceNodesPerFace);
  varYVB->put(&(m_dTargetVertexLat[0][0]), nB, nTargetNodesPerFace);

  varXVA->put(&(m_dSourceVertexLon[0][0]), nA, nSourceNodesPerFace);
  varXVB->put(&(m_dTargetVertexLon[0][0]), nB, nTargetNodesPerFace);

  // Write vector centers
  if ((m_dVectorTargetCenterLat.GetRows() != 0) &&
    (m_dVectorTargetCenterLon.GetRows() != 0)
  ) {
    NcDim * dimLatB =
      ncMap.add_dim("lat_b", m_dVectorTargetCenterLat.GetRows());
    NcDim * dimLonB =
      ncMap.add_dim("lon_b", m_dVectorTargetCenterLon.GetRows());

    NcVar * varLatCB = ncMap.add_var("latc_b", ncDouble, dimLatB);
    NcVar * varLonCB = ncMap.add_var("lonc_b", ncDouble, dimLonB);

    varLatCB->put(&(m_dVectorTargetCenterLat[0]), dimLatB->size());
    varLonCB->put(&(m_dVectorTargetCenterLon[0]), dimLonB->size());

    NcDim * dimBounds = ncMap.add_dim("bnds", 2);
    NcVar * varLatBounds =
      ncMap.add_var("lat_bnds", ncDouble, dimLatB, dimBounds);
    NcVar * varLonBounds =
      ncMap.add_var("lon_bnds", ncDouble, dimLonB, dimBounds);

    varLatBounds->put(&(m_dVectorTargetBoundsLat[0][0]),
      m_dVectorTargetBoundsLat.GetRows(), 2);
    varLonBounds->put(&(m_dVectorTargetBoundsLon[0][0]),
      m_dVectorTargetBoundsLon.GetRows(), 2);
  }

  // Write areas
  NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
  varAreaA->put(&(m_dSourceAreas[0]), nA);
  varAreaA->add_att("units", "steradians");

  NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
  varAreaB->put(&(m_dTargetAreas[0]), nB);
  varAreaB->add_att("units", "steradians");

  // Write masks
  if (m_iSourceMask.IsAttached()) {
    NcVar * varMaskA = ncMap.add_var("mask_a", ncInt, dimNA);
    varMaskA->put(&(m_iSourceMask[0]), nA);
    varMaskA->add_att("units", "unitless");

    if (!m_iTargetMask.IsAttached()) {
      NcVar * varMaskB = ncMap.add_var("mask_b", ncInt, dimNB);
      DataArray1D<int> iTargetMaskTemp(nB);
      for (int i = 0; i < nB; i++) {
        iTargetMaskTemp[i] = 1;
      }
      varMaskB->put(&(iTargetMaskTemp[0]), nB);
      varMaskB->add_att("units", "unitless");
    }
  }

  if (m_iTargetMask.IsAttached()) {
    if (!m_iSourceMask.IsAttached()) {
      NcVar * varMaskA = ncMap.add_var("mask_a", ncInt, dimNA);
      DataArray1D<int> iSourceMaskTemp(nA);
      for (int i = 0; i < nA; i++) {
        iSourceMaskTemp[i] = 1;
      }
      varMaskA->put(&(iSourceMaskTemp[0]), nA);
      varMaskA->add_att("units", "unitless");
    }

    NcVar * varMaskB = ncMap.add_var("mask_b", ncInt, dimNB);
    varMaskB->put(&(m_iTargetMask[0]), nB);
    varMaskB->add_att("units", "unitless");
  }

  // Write SparseMatrix entries
  DataArray1D<int> vecRow;
  DataArray1D<int> vecCol;
  DataArray1D<double> vecS;

  m_mapRemap.GetEntries(vecRow, vecCol, vecS);

  // Calculate and write fractional coverage arrays
  {
    DataArray1D<double> dFracA(nA);
    DataArray1D<double> dFracB(nB);
  
    for (int i = 0; i < vecS.GetRows(); i++) {
      dFracA[vecCol[i]] += vecS[i] / m_dSourceAreas[vecCol[i]] * m_dTargetAreas[vecRow[i]];
      dFracB[vecRow[i]] += vecS[i];
    }

    NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
    varFracA->put(&(dFracA[0]), nA);
    varFracA->add_att("name", "fraction of target coverage of source dof");
    varFracA->add_att("units", "unitless");

    NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
    varFracB->put(&(dFracB[0]), nB);
    varFracB->add_att("name", "fraction of source coverage of target dof");
    varFracB->add_att("units", "unitless");
  }

  // Increment vecRow and vecCol
  for (int i = 0; i < vecRow.GetRows(); i++) {
    vecRow[i]++;
    vecCol[i]++;
  }

  // Write out data
  int nS = vecRow.GetRows();
  NcDim * dimNS = ncMap.add_dim("n_s", nS);

  NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
  varRow->add_att("name", "sparse matrix target dof index");
  varRow->add_att("first_index", "1");

  NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
  varCol->add_att("name", "sparse matrix source dof index");
  varCol->add_att("first_index", "1");

  NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);
  varS->add_att("name", "sparse matrix coefficient");

  varRow->set_cur((long)0);
  varRow->put(&(vecRow[0]), nS);

  varCol->set_cur((long)0);
  varCol->put(&(vecCol[0]), nS);

  varS->set_cur((long)0);
  varS->put(&(vecS[0]), nS);

  // Add global attributes
  std::map<std::string, std::string>::const_iterator iterAttributes =
    mapAttributes.begin();
  for (; iterAttributes != mapAttributes.end(); iterAttributes++) {
    ncMap.add_att(
      iterAttributes->first.c_str(),
      iterAttributes->second.c_str());
  }
}
#endif

template<typename T>
ErrorCode get_vartag_data(moab::Interface* mbCore, Tag tag, moab::Range sets, int data_size, std::vector<T>& data)
{
  int* tag_sizes = new int [ sets.size() ];
  const void** tag_data = (const void**) new void* [sets.size()];

  ErrorCode rval = mbCore->tag_get_by_ptr(tag, sets, tag_data, tag_sizes );MB_CHK_SET_ERR(rval, "Getting matrix rows failed");

  data.reserve(data_size);
  int ioffset = 0, index = 0;
  for (Range::iterator it=sets.begin(); it!=sets.end(); it++, index++)
  {
    T* m_vals = (T*)tag_data[index];
    for (int k=0; k< tag_sizes[index]; k++, ioffset++)
    {
      assert(ioffset < data_size);
      data[ioffset] = m_vals[k];
    }
  }

  return moab::MB_SUCCESS;
}

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;
  int dimension=2;
  NcError error2 ( NcError::verbose_nonfatal );
  std::stringstream sstr;
  ProgOptions opts;
  std::string h5mfilename, scripfile;
  bool noMap=false;

#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

  opts.addOpt<std::string>("weights,w", "h5m remapping weights filename", &h5mfilename);
  opts.addOpt<std::string>("scrip,s", "Output SCRIP map filename", &scripfile);
  opts.addOpt<int>("dim,d", "Dimension of entities to use for partitioning", &dimension);
  opts.addOpt<void>("mesh,m", "Only convert the mesh and exclude the remap weight details", &noMap);
  opts.parseCommandLine(argc, argv);

  if (h5mfilename.empty() || scripfile.empty())
  {
    opts.printHelp();
    exit(1);
  }

  moab::Interface* mbCore = new ( std::nothrow ) moab::Core;

  if ( NULL == mbCore ) { return 1; }

  //Set the read options for parallel file loading
  const std::string partition_set_name = "PARALLEL_PARTITION";
  const std::string global_id_name = "GLOBAL_ID";

  //Load file
  rval = mbCore->load_mesh(h5mfilename.c_str());MB_CHK_ERR(rval);

  // moab::Range elems, verts;

  // rval = mbCore->get_entities_by_dimension(0, 2, elems);MB_CHK_ERR(rval);
  // rval = mbCore->get_entities_by_dimension(0, 0, verts);MB_CHK_ERR(rval);

  try {

    // Temporarily change rval reporting
    NcError error_temp(NcError::verbose_fatal);

    // Open an output file
    NcFile ncMap(scripfile.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
    if (!ncMap.is_valid()) {
      _EXCEPTION1("Unable to open output map file \"%s\"",
        scripfile.c_str());
    }

    // Attributes
    ncMap.add_att("Title", "MOAB-TempestRemap (MBTR) Offline Regridding Weight Converter (h5mtoscrip)");

    moab::EntityHandle rootset = 0;
    Tag smatMetadataTag;
    int smat_metadata_glb[7];
    rval = mbCore->tag_get_handle( "SMAT_DATA" , 7, MB_TYPE_INTEGER, smatMetadataTag, MB_TAG_SPARSE);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_data(smatMetadataTag, &rootset, 1, smat_metadata_glb);MB_CHK_ERR(rval);

    // Map dimensions
    int nA = smat_metadata_glb[0];
    int nB = smat_metadata_glb[1];
    int nVA = smat_metadata_glb[2];
    int nVB = smat_metadata_glb[3];
    int nDofB = smat_metadata_glb[4];
    int nDofA = smat_metadata_glb[5];
    int NNZ = smat_metadata_glb[6];

    Tag srcIDTag, srcAreaTag, tgtIDTag, tgtAreaTag;
    Tag smatRowdataTag, smatColdataTag, smatValsdataTag;
    rval = mbCore->tag_get_handle( "SourceGIDS" , srcIDTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SourceAreas" , srcAreaTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetGIDS" , tgtIDTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetAreas" , tgtAreaTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SMAT_ROWS" , smatRowdataTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SMAT_COLS" , smatColdataTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SMAT_VALS" , smatValsdataTag);MB_CHK_ERR(rval);

    // Get sets entities, by type
    moab::Range sets;
    // rval = mbCore->get_entities_by_type(0, MBENTITYSET, sets);MB_CHK_ERR(rval);
    rval = mbCore->get_entities_by_type_and_tag(0, MBENTITYSET, &smatRowdataTag, NULL, 1, sets, moab::Interface::UNION, true);MB_CHK_ERR(rval);

    std::vector<int> src_gids, tgt_gids;
    std::vector<double> src_areas, tgt_areas;
    rval = get_vartag_data(mbCore, srcIDTag, sets, nDofA*sets.size(), src_gids);MB_CHK_SET_ERR(rval, "Getting source mesh IDs failed");
    rval = get_vartag_data(mbCore, tgtIDTag, sets, nDofB*sets.size(), tgt_gids);MB_CHK_SET_ERR(rval, "Getting target mesh IDs failed");
    rval = get_vartag_data(mbCore, srcAreaTag, sets, nDofA*sets.size(), src_areas);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
    rval = get_vartag_data(mbCore, tgtAreaTag, sets, nDofB*sets.size(), tgt_areas);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");

    std::vector<double> src_glob_areas(nDofA), tgt_glob_areas(nDofB);
    for (unsigned i=0; i < nDofA*sets.size(); ++i)
      src_glob_areas[src_gids[i]] = src_areas[i];
    for (unsigned i=0; i < nDofB*sets.size(); ++i)
      tgt_glob_areas[tgt_gids[i]] = tgt_areas[i];

    // Output the number of sets    
    std::cout << "Number of primary sets is " << sets.size() << std::endl;

    // Write output dimensions entries
    int nSrcGridDims = 1;
    int nDstGridDims = 1;

    NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
    NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

    NcVar * varSrcGridDims =
      ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
    NcVar * varDstGridDims =
      ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

    if (nA == nDofA) {
      varSrcGridDims->put(&nA, 1);
      varSrcGridDims->add_att("name0", "num_elem");
    }
    else {
      varSrcGridDims->put(&nDofA, 1);
      varSrcGridDims->add_att("name1", "num_dof");
    }

    if (nB == nDofB) {
      varDstGridDims->put(&nB, 1);
      varDstGridDims->add_att("name0", "num_elem");
    }
    else {
      varDstGridDims->put(&nDofB, 1);
      varDstGridDims->add_att("name1", "num_dof");
    }

    // Source and Target mesh resolutions
    NcDim * dimNA = ncMap.add_dim("n_a", nDofA);
    NcDim * dimNB = ncMap.add_dim("n_b", nDofB);

    // Source and Target verticecs per elements
    NcDim * dimNVA = ncMap.add_dim("nv_a", nVA);
    NcDim * dimNVB = ncMap.add_dim("nv_b", nVB);

    // Write coordinates
    NcVar * varYCA = ncMap.add_var("yc_a", ncDouble, dimNA);
    NcVar * varYCB = ncMap.add_var("yc_b", ncDouble, dimNB);

    NcVar * varXCA = ncMap.add_var("xc_a", ncDouble, dimNA);
    NcVar * varXCB = ncMap.add_var("xc_b", ncDouble, dimNB);

    NcVar * varYVA = ncMap.add_var("yv_a", ncDouble, dimNA, dimNVA);
    NcVar * varYVB = ncMap.add_var("yv_b", ncDouble, dimNB, dimNVB);

    NcVar * varXVA = ncMap.add_var("xv_a", ncDouble, dimNA, dimNVA);
    NcVar * varXVB = ncMap.add_var("xv_b", ncDouble, dimNB, dimNVB);

    varYCA->add_att("units", "degrees");
    varYCB->add_att("units", "degrees");

    varXCA->add_att("units", "degrees");
    varXCB->add_att("units", "degrees");

    varYVA->add_att("units", "degrees");
    varYVB->add_att("units", "degrees");

    varXVA->add_att("units", "degrees");
    varXVB->add_att("units", "degrees");

    // Write areas
    NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
    varAreaA->put(&(src_glob_areas[0]), nDofA);
    varAreaA->add_att("units", "steradians");

    NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
    varAreaB->put(&(tgt_glob_areas[0]), nDofB);
    varAreaB->add_att("units", "steradians");

    std::vector<int> mat_rows, mat_cols;
    std::vector<double> mat_vals;
    rval = get_vartag_data(mbCore, smatRowdataTag, sets, NNZ, mat_rows);MB_CHK_SET_ERR(rval, "Getting source mesh IDs failed");
    rval = get_vartag_data(mbCore, smatColdataTag, sets, NNZ, mat_cols);MB_CHK_SET_ERR(rval, "Getting source mesh IDs failed");
    rval = get_vartag_data(mbCore, smatValsdataTag, sets, NNZ, mat_vals);MB_CHK_SET_ERR(rval, "Getting source mesh IDs failed");

    // Write out data
    NcDim * dimNS = ncMap.add_dim("n_s", NNZ);

    NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
    varRow->add_att("name", "sparse matrix target dof index");
    varRow->add_att("first_index", "1");

    NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
    varCol->add_att("name", "sparse matrix source dof index");
    varCol->add_att("first_index", "1");

    NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);
    varS->add_att("name", "sparse matrix coefficient");

    varRow->set_cur((long)0);
    varRow->put(&(mat_rows[0]), NNZ);

    varCol->set_cur((long)0);
    varCol->put(&(mat_cols[0]), NNZ);

    varS->set_cur((long)0);
    varS->put(&(mat_vals[0]), NNZ);

    // Add global attributes
    std::map<std::string, std::string> mapAttributes;
    mapAttributes["Command"] = "Converted with MBTR:h5mtoscrip";
    std::map<std::string, std::string>::const_iterator iterAttributes =
      mapAttributes.begin();
    for (; iterAttributes != mapAttributes.end(); iterAttributes++) {
      ncMap.add_att(
        iterAttributes->first.c_str(),
        iterAttributes->second.c_str());
    }

    // rval = mbCore->write_file(scripfile.c_str());MB_CHK_ERR(rval);
  }
  catch (std::exception & e)
  {
    std::cout << " exception caught during tree initialization " << e.what() << std::endl;
  }
  delete mbCore;

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

  exit(0);
}
