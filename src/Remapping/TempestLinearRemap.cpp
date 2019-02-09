///////////////////////////////////////////////////////////////////////////////
///
/// \file    TempestLinearRemap.cpp
/// \author  Vijay Mahadevan
/// \version Mar 08, 2017
///

#include "Announce.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
// #include "LinearRemapFV.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MeshUtilitiesFuzzy.h"
#include "MeshUtilitiesExact.h"
#include "MathHelper.h"
#include "SparseMatrix.h"
#include "OverlapMesh.h"

#include "moab/Remapping/TempestOnlineMap.hpp"
#include "DebugOutput.hpp"

#include "netcdfcpp.h"

#ifdef MOAB_HAVE_EIGEN
#include <Eigen/Dense>
#endif

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>

// #define VERBOSE


extern void BuildIntegrationArray (
    const Mesh & m_meshInput,
    const Mesh & m_meshOverlap,
    const TriangularQuadratureRule & triquadrule,
    int ixFirstFace,
    int ixOverlapBegin,
    int ixOverlapEnd,
    int nOrder,
    DataArray2D<double> & dIntArray
);

extern void InvertFitArray_Corrected (
    const DataArray1D<double> & dConstraint,
    DataArray2D<double> & dFitArray,
    DataArray1D<double> & dFitWeights,
    DataArray2D<double> & dFitArrayPlus
);

/// <summary>
///     Face index and distance metric pair.
/// </summary>
typedef std::pair<int, int> FaceDistancePair;

/// <summary>
///     Vector storing adjacent Faces.
/// </summary>
typedef std::vector<FaceDistancePair> AdjacentFaceVector;

extern void BuildFitArray (
    const Mesh & mesh,
    const TriangularQuadratureRule & triquadrule,
    int ixFirst,
    const AdjacentFaceVector & vecAdjFaces,
    int nOrder,
    int nFitWeightsExponent,
    const DataArray1D<double> & dConstraint,
    DataArray2D<double> & dFitArray,
    DataArray1D<double> & dFitWeights
);

extern void GetAdjacentFaceVectorByEdge (
    const Mesh & mesh,
    int iFaceInitial,
    int nRequiredFaceSetSize,
    AdjacentFaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapFVtoFV_Tempest_MOAB (
    int nOrder
)
{
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if ( m_meshInputCov->revnodearray.size() == 0 )
    {
        _EXCEPTIONT ( "ReverseNodeArray has not been calculated for m_meshInput" );
    }

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule ( TriQuadRuleOrder );

    // Number of coefficients needed at this order
#ifdef RECTANGULAR_TRUNCATION
    int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION
    int nCoefficients = nOrder * ( nOrder + 1 ) / 2;
#endif

    // Number of faces you need
    const int nRequiredFaceSetSize = nCoefficients;

    // Fit weight exponent
    const int nFitWeightsExponent = nOrder + 2;

    // Announcemnets
    if ( is_root )
    {
        Announce ( "[moab::TempestOnlineMap::LinearRemapFVtoFV_Tempest_MOAB] Finite Volume to Finite Volume Projection" );
        Announce ( "Triangular quadrature rule order %i", TriQuadRuleOrder );
        Announce ( "Number of coefficients: %i", nCoefficients );
        Announce ( "Required adjacency set size: %i", nRequiredFaceSetSize );
        Announce ( "Fit weights exponent: %i", nFitWeightsExponent );
    }

    // Current overlap face
    int ixOverlap = 0;

#ifdef VERBOSE
    const unsigned outputFrequency = (m_meshInputCov->faces.size()/10);
#endif
    // Loop through all faces on m_meshInput
    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Output every 1000 elements
        if ( ixFirst % outputFrequency == 0 )
        {
            Announce ( "Element %i/%i", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // This Face

        // Find the set of Faces that overlap faceFirst
        int ixOverlapBegin = ixOverlap;
        unsigned ixOverlapEnd = ixOverlapBegin;

        for ( ; ixOverlapEnd < m_meshOverlap->faces.size(); ixOverlapEnd++ )
        {
            if ( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapEnd] != 0 )
            {
                break;
            }
        }

        unsigned nOverlapFaces = ixOverlapEnd - ixOverlapBegin;
        // if ( is_root ) Announce ( "Element %i / %i :: [%i, %i]", ixFirst, m_meshInputCov->faces.size(), ixOverlapBegin, ixOverlapEnd );

        if ( nOverlapFaces == 0 ) continue;

        // Build integration array
        DataArray2D<double> dIntArray;

        BuildIntegrationArray (
            *m_meshInputCov,
            *m_meshOverlap,
            triquadrule,
            ixFirst,
            ixOverlapBegin,
            ixOverlapEnd,
            nOrder,
            dIntArray );

        // Set of Faces to use in building the reconstruction and associated
        // distance metric.
        AdjacentFaceVector vecAdjFaces;

        GetAdjacentFaceVectorByEdge (
            *m_meshInputCov,
            ixFirst,
            nRequiredFaceSetSize,
            vecAdjFaces );

        // Number of adjacent Faces
        int nAdjFaces = vecAdjFaces.size();

        // Determine the conservative constraint equation
        DataArray1D<double> dConstraint( nCoefficients );

        double dFirstArea = m_meshInputCov->vecFaceArea[ixFirst];

        for ( int p = 0; p < nCoefficients; p++ )
        {
            for ( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                dConstraint[p] += dIntArray[p][j];
            }
            dConstraint[p] /= dFirstArea;
        }

        // Build the fit array from the integration operator
        DataArray2D<double> dFitArray;
        DataArray1D<double> dFitWeights;
        DataArray2D<double> dFitArrayPlus;

        BuildFitArray (
            *m_meshInputCov,
            triquadrule,
            ixFirst,
            vecAdjFaces,
            nOrder,
            nFitWeightsExponent,
            dConstraint,
            dFitArray,
            dFitWeights
        );

        // Compute the inverse fit array
        InvertFitArray_Corrected (
            dConstraint,
            dFitArray,
            dFitWeights,
            dFitArrayPlus
        );

        // Multiply integration array and fit array
        DataArray2D<double> dComposedArray( nAdjFaces, nOverlapFaces );

        for ( int i = 0; i < nAdjFaces; i++ )
        {
            for ( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                for ( int k = 0; k < nCoefficients; k++ )
                {
                    dComposedArray[i][j] += dIntArray[k][j] * dFitArrayPlus[i][k];
                }
            }
        }

        // Put composed array into map
        for ( unsigned i = 0; i < vecAdjFaces.size(); i++ )
        {
            for ( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                int& ixFirstFaceLoc = vecAdjFaces[i].first;
                int& ixSecondFaceLoc = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];
                // int ixFirstFaceGlob = m_remapper->GetGlobalID(moab::Remapper::SourceMesh, ixFirstFaceLoc);
                // int ixSecondFaceGlob = m_remapper->GetGlobalID(moab::Remapper::TargetMesh, ixSecondFaceLoc);

                m_mapRemap ( ixSecondFaceLoc, ixFirstFaceLoc ) +=
                    dComposedArray[i][j]
                    / m_meshOutput->vecFaceArea[ixSecondFaceLoc];
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}


#ifdef MOAB_HAVE_EIGEN
void moab::TempestOnlineMap::CopyTempestSparseMat_Eigen()
{
    m_weightMatrix.resize(m_nTotDofs_Dest, m_nTotDofs_SrcCov);
    InitVectors();

#ifdef VERBOSE
    int locrows = std::max(m_mapRemap.GetRows(), m_nTotDofs_Dest);
    int loccols = std::max(m_mapRemap.GetColumns(), m_nTotDofs_SrcCov);

    std::cout << m_weightMatrix.rows() << ", " <<  locrows << ", " <<  m_weightMatrix.cols() << ", " << loccols << "\n";
    // assert(m_weightMatrix.rows() == locrows && m_weightMatrix.cols() == loccols);
#endif

    DataArray1D<int> lrows;
    DataArray1D<int> lcols;
    DataArray1D<double> lvals;
    m_mapRemap.GetEntries(lrows, lcols, lvals);
    unsigned locvals = lvals.GetRows();

    m_weightMatrix.reserve(locvals);
    for (unsigned iv=0; iv < locvals; iv++) {
        m_weightMatrix.insert(lrows[iv], lcols[iv]) = lvals[iv];
    }

    m_weightMatrix.makeCompressed();

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "tempestmatrix.txt.0000" << rank;
    std::ofstream output_file ( sstr.str(), std::ios::out );
    output_file << "0 " << locrows << " 0 " << loccols << "\n";
    for (unsigned iv=0; iv < locvals; iv++) {
        // output_file << lrows[iv] << " " << row_ldofmap[lrows[iv]] << " " << row_gdofmap[row_ldofmap[lrows[iv]]] << " " << col_gdofmap[col_ldofmap[lcols[iv]]] << " " << lvals[iv] << "\n";
        output_file << row_gdofmap[row_ldofmap[lrows[iv]]] << " " << col_gdofmap[col_ldofmap[lcols[iv]]] << " " << lvals[iv] << "\n";
        
    }
    output_file.flush(); // required here
    output_file.close();
#endif

    return;
}

///////////////////////////////////////////////////////////////////////////////

// #define IO_USE_PARALLEL_NETCDF
void moab::TempestOnlineMap::WriteParallelWeightsToFile(std::string strFilename)
{
    // m_weightMatrix.Print(filename.c_str(), 0, 0);

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiFile ncMap(MPI_COMM_WORLD, strFilename.c_str(), NcmpiFile::Replace, NcmpiFile::classic5);
#else
    NcFile ncMap(strFilename.c_str(), NcFile::Replace);
#endif

    if (!ncMap.is_valid()) {
        _EXCEPTION1("Unable to open output map file \"%s\"",
            strFilename.c_str());
    }

    // Attributes
    ncMap.add_att("Title", "MOAB-TempestRemap Online Regridding Weight Generator");

    // Map dimensions
    unsigned nA = (m_dSourceAreas.GetRows());
    unsigned nB = (m_dTargetAreas.GetRows());

    // Write output dimensions entries
    unsigned nSrcGridDims = (m_vecSourceDimSizes.size());
    unsigned nDstGridDims = (m_vecTargetDimSizes.size());

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiDim * dimSrcGridRank = ncMap.addDim("src_grid_rank", nSrcGridDims);
    NcmpiDim * dimDstGridRank = ncMap.addDim("dst_grid_rank", nDstGridDims);
#else
    NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
    NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);
#endif

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varSrcGridDims =
        ncMap.addVar("src_grid_dims", ncmpiInt, dimSrcGridRank);
    NcmpiVar * varDstGridDims =
        ncMap.addVar("dst_grid_dims", ncmpiInt, dimDstGridRank);
#else
    NcVar * varSrcGridDims =
        ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
    NcVar * varDstGridDims =
        ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);
#endif

    char szDim[64];
    if ((nSrcGridDims == 1) && (m_vecSourceDimSizes[0] != (int)nA)) {
        int tmp = (int)(nA);
        varSrcGridDims->put(&tmp, 1);
        varSrcGridDims->add_att("name0", "num_dof");

    } else {
        for (unsigned i = 0; i < m_vecSourceDimSizes.size(); i++) {
            varSrcGridDims->set_cur(nSrcGridDims - i - 1);
            varSrcGridDims->put(&(m_vecSourceDimSizes[i]), 1);
        }

        for (unsigned i = 0; i < m_vecSourceDimSizes.size(); i++) {
            sprintf(szDim, "name%i", i);
            varSrcGridDims->add_att(szDim,
                m_vecSourceDimNames[nSrcGridDims - i - 1].c_str());
        }
    }

    if ((nDstGridDims == 1) && (m_vecTargetDimSizes[0] != (int)nB)) {
        int tmp = (int)(nB);
        varDstGridDims->put(&tmp, 1);
        varDstGridDims->add_att("name0", "num_dof");
    } else {
        for (unsigned i = 0; i < m_vecTargetDimSizes.size(); i++) {
            varDstGridDims->set_cur(nDstGridDims - i - 1);
            varDstGridDims->put(&(m_vecTargetDimSizes[i]), 1);
        }

        for (unsigned i = 0; i < m_vecTargetDimSizes.size(); i++) {
            sprintf(szDim, "name%i", i);
            varDstGridDims->add_att(szDim,
                m_vecTargetDimNames[nDstGridDims - i - 1].c_str());
        }
    }

    // Source and Target mesh resolutions
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiDim * dimNA = ncMap.addDim("n_a", nA);
    NcmpiDim * dimNB = ncMap.addDim("n_b", nB);
#else
    NcDim * dimNA = ncMap.add_dim("n_a", nA);
    NcDim * dimNB = ncMap.add_dim("n_b", nB);
#endif

    // Number of nodes per Face
    int nSourceNodesPerFace = m_dSourceVertexLon.GetColumns();
    int nTargetNodesPerFace = m_dTargetVertexLon.GetColumns();

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiDim * dimNVA = ncMap.addDim("nv_a", nSourceNodesPerFace);
    NcmpiDim * dimNVB = ncMap.addDim("nv_b", nTargetNodesPerFace);
#else
    NcDim * dimNVA = ncMap.add_dim("nv_a", nSourceNodesPerFace);
    NcDim * dimNVB = ncMap.add_dim("nv_b", nTargetNodesPerFace);
#endif

    // Write coordinates
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varYCA = ncMap.addVar("yc_a", ncmpiDouble, dimNA);
    NcmpiVar * varYCB = ncMap.addVar("yc_b", ncmpiDouble, dimNB);

    NcmpiVar * varXCA = ncMap.addVar("xc_a", ncmpiDouble, dimNA);
    NcmpiVar * varXCB = ncMap.addVar("xc_b", ncmpiDouble, dimNB);

    NcmpiVar * varYVA = ncMap.addVar("yv_a", ncmpiDouble, dimNA, dimNVA);
    NcmpiVar * varYVB = ncMap.addVar("yv_b", ncmpiDouble, dimNB, dimNVB);

    NcmpiVar * varXVA = ncMap.addVar("xv_a", ncmpiDouble, dimNA, dimNVA);
    NcmpiVar * varXVB = ncMap.addVar("xv_b", ncmpiDouble, dimNB, dimNVB);
#else
    NcVar * varYCA = ncMap.add_var("yc_a", ncDouble, dimNA);
    NcVar * varYCB = ncMap.add_var("yc_b", ncDouble, dimNB);

    NcVar * varXCA = ncMap.add_var("xc_a", ncDouble, dimNA);
    NcVar * varXCB = ncMap.add_var("xc_b", ncDouble, dimNB);

    NcVar * varYVA = ncMap.add_var("yv_a", ncDouble, dimNA, dimNVA);
    NcVar * varYVB = ncMap.add_var("yv_b", ncDouble, dimNB, dimNVB);

    NcVar * varXVA = ncMap.add_var("xv_a", ncDouble, dimNA, dimNVA);
    NcVar * varXVB = ncMap.add_var("xv_b", ncDouble, dimNB, dimNVB);
#endif

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
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varAreaA = ncMap.addVar("area_a", ncmpiDouble, dimNA);
#else
    NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
#endif
    varAreaA->put(&(m_dSourceAreas[0]), nA);

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varAreaB = ncMap.addVar("area_b", ncmpiDouble, dimNB);
#else
    NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
#endif
    varAreaB->put(&(m_dTargetAreas[0]), nB);

    // Write frac
    DataArray1D<double> dFracA(nA);
    for (unsigned i = 0; i < nA; i++) {
        dFracA[i] = 1.0;
    }
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varFracA = ncMap.addVar("frac_a", ncmpiDouble, dimNA);
#else
    NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
#endif
    varFracA->put(&(dFracA[0]), nA);

    DataArray1D<double> dFracB(nB);
    for (unsigned i = 0; i < nB; i++) {
        dFracB[i] = 1.0;
    }
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varFracB = ncMap.addVar("frac_b", ncmpiDouble, dimNB);
#else
    NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
#endif
    varFracB->put(&(dFracB[0]), nB);

    // Write SparseMatrix entries
    int nS = m_weightMatrix.nonZeros();
    DataArray1D<int> vecRow(nS);
    DataArray1D<int> vecCol(nS);
    DataArray1D<double> vecS(nS);

    for (int i = 0; i < m_weightMatrix.outerSize(); i++) {
        for (WeightMatrix::InnerIterator it(m_weightMatrix,i); it; ++it) {
            vecRow[i] = 1+it.row(); // row index
            vecCol[i] = 1+it.col(); // col index
            vecS[i] = it.value();   // value
        }
    }

    // Load in data
#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiDim * dimNS = ncMap.add_dim("n_s", nS);
#else
    NcDim * dimNS = ncMap.add_dim("n_s", nS);
#endif

#ifdef IO_USE_PARALLEL_NETCDF
    NcmpiVar * varRow = ncMap.addVar("row", ncmpiInt, dimNS);
    NcmpiVar * varCol = ncMap.addVar("col", ncmpiInt, dimNS);
    NcmpiVar * varS = ncMap.addVar("S", ncmpiDouble, dimNS);
#else
    NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
    NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
    NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);
#endif

    varRow->set_cur((long)0);
    varRow->put(&(vecRow[0]), nS);

    varCol->set_cur((long)0);
    varCol->put(&(vecCol[0]), nS);

    varS->set_cur((long)0);
    varS->put(&(vecS[0]), nS);

    // Add global attributes
    // std::map<std::string, std::string>::const_iterator iterAttributes =
    //     mapAttributes.begin();
    // for (; iterAttributes != mapAttributes.end(); iterAttributes++) {
    //     ncMap.add_att(
    //         iterAttributes->first.c_str(),
    //         iterAttributes->second.c_str());
    // }
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::ApplyWeights (std::vector<double>& srcVals, std::vector<double>& tgtVals, bool transpose)
{
    // Reset the source and target data first
    m_rowVector.setZero();
    m_colVector.setZero();

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "projection_" << rank << ".txt";
    std::ofstream output_file ( sstr.str() );
#endif
    // Perform the actual projection of weights: application of weight matrix onto the source solution vector
    if (transpose) {
        // Permute the source data first
        for (unsigned i=0; i < srcVals.size(); ++i) {
            m_rowVector(row_dofmap[i]) = srcVals[i]; // permute and set the row (source) vector properly
        }
        
        m_colVector = m_weightMatrix.adjoint() * m_rowVector;

        // Permute the resulting target data back
        for (unsigned i=0; i < tgtVals.size(); ++i) {
            tgtVals[i] = m_colVector(col_dofmap[i]); // permute and set the row (source) vector properly
        }
    }
    else {
        // Permute the source data first
#ifdef VERBOSE
        output_file << "ColVector: " << m_colVector.size() << ", SrcVals: " << srcVals.size() << ", Sizes: " << m_nTotDofs_SrcCov << ", " << col_dofmap.size() << "\n";
#endif
        for (unsigned i=0; i < srcVals.size(); ++i) {
            assert(m_colVector.size()-col_dofmap[i]>0);
            m_colVector(col_dofmap[i]) = srcVals[i]; // permute and set the row (source) vector properly
#ifdef VERBOSE
            output_file << "Col: " << i << ", " << col_dofmap[i] << ", GID: " << col_gdofmap[i] << ", Data = " << srcVals[i]  << ", " << m_colVector(col_dofmap[i]) << "\n";
#endif
        }
        
        m_rowVector = m_weightMatrix * m_colVector;

        // Permute the resulting target data back
#ifdef VERBOSE
        output_file << "RowVector: " << m_rowVector.size() << ", TgtVals:" << tgtVals.size() << ", Sizes: " << m_nTotDofs_Dest << ", " << row_dofmap.size() << "\n";
#endif
        for (unsigned i=0; i < tgtVals.size(); ++i) {
            tgtVals[i] = m_rowVector(row_dofmap[i]); // permute and set the row (source) vector properly
#ifdef VERBOSE
            output_file << "Row: " << i << ", " << row_dofmap[i] << ", GID: " << row_gdofmap[i] << ", Data = " << m_rowVector(row_dofmap[i]) << "\n";
#endif
        }
    }

#ifdef VERBOSE
    output_file.flush(); // required here
    output_file.close();
#endif

    // All done with matvec application
    return moab::MB_SUCCESS;
}

#endif


///////////////////////////////////////////////////////////////////////////////

extern void ForceConsistencyConservation3(
    const DataArray1D<double> & vecSourceArea,
    const DataArray1D<double> & vecTargetArea,
    DataArray2D<double> & dCoeff,
    bool fMonotone
);

///////////////////////////////////////////////////////////////////////////////

extern void ForceIntArrayConsistencyConservation (
    const DataArray1D<double> & vecSourceArea,
    const DataArray1D<double> & vecTargetArea,
    DataArray2D<double> & dCoeff,
    bool fMonotone
);

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapSE4_Tempest_MOAB (
    const DataArray3D<int> & dataGLLNodes,
    const DataArray3D<double> & dataGLLJacobian,
    int nMonotoneType,
    bool fContinuousIn,
    bool fNoConservation
)
{
    // Order of the polynomial interpolant
    int nP = dataGLLNodes.GetRows();

    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule ( TriQuadRuleOrder );

    int TriQuadraturePoints = triquadrule.GetPoints();

    const DataArray2D<double> & TriQuadratureG = triquadrule.GetG();

    const DataArray1D<double> & TriQuadratureW = triquadrule.GetW();

    // Sample coefficients
    DataArray2D<double> dSampleCoeff( nP, nP );

    // GLL Quadrature nodes on quadrilateral elements
    DataArray1D<double> dG;
    DataArray1D<double> dW;
    GaussLobattoQuadrature::GetPoints ( nP, 0.0, 1.0, dG, dW );

    // Announcemnets
    if ( is_root )
    {
        Announce ( "[moab::TempestOnlineMap::LinearRemapSE4_Tempest_MOAB] Finite Element to Finite Volume Projection" );
        Announce ( "Triangular quadrature rule order %i", TriQuadRuleOrder );
        Announce ( "Order of the FE polynomial interpolant: %i", nP );
    }

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix<double> & smatMap = this->GetSparseMatrix();

    // NodeVector from m_meshOverlap
    const NodeVector & nodesOverlap = m_meshOverlap->nodes;
    const NodeVector & nodesFirst   = m_meshInputCov->nodes;

    // Vector of source areas
    DataArray1D<double> vecSourceArea( nP * nP );

    DataArray1D<double> vecTargetArea;
    DataArray2D<double> dCoeff;

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "remapdata_" << rank << ".txt";
    std::ofstream output_file ( sstr.str() );
#endif

    // Current Overlap Face
    int ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = (m_meshInputCov->faces.size()/10);
#endif

    // Loop over all input Faces
    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        const Face & faceFirst = m_meshInputCov->faces[ixFirst];

        if ( faceFirst.edges.size() != 4 )
        {
            _EXCEPTIONT ( "Only quadrilateral elements allowed for SE remapping" );
        }

#ifdef VERBOSE
        // Announce computation progress
        if ( ixFirst % outputFrequency == 0 )
        {
            Announce ( "Element %i/%i", ixFirst, m_meshInputCov->faces.size() );
        }
#endif

        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so on wrt the input mesh data
        // Then the overlap_end and overlap_begin will be correct. However, the relation with MOAB and Tempest will go out of the roof

        // Determine how many overlap Faces and triangles are present
        int nOverlapFaces = 0;
        size_t ixOverlapTemp = ixOverlap;
        for ( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];
            if ( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 )
            {
                break;
            }

            nOverlapFaces++;
        }

        // No overlaps
        if ( nOverlapFaces == 0 )
        {
            continue;
        }

        // Allocate remap coefficients array for meshFirst Face
        DataArray3D<double> dRemapCoeff( nP, nP, nOverlapFaces );

        // Find the local remap coefficients
        for ( int j = 0; j < nOverlapFaces; j++ )
        {
            const Face & faceOverlap = m_meshOverlap->faces[ixOverlap + j];

// #ifdef VERBOSE
            // if ( is_root )
            //     Announce ( "\tLocal ID: %i/%i = %i, areas = %2.8e", j + ixOverlap, nOverlapFaces, m_remapper->lid_to_gid_covsrc[m_meshOverlap->vecSourceFaceIx[ixOverlap + j]], m_meshOverlap->vecFaceArea[ixOverlap + j] );
// #endif

            int nOverlapTriangles = faceOverlap.edges.size() - 2;

// #define USE_MININDEX

#ifdef USE_MININDEX
            // first find out the minimum node, start there the triangle decomposition
            int minIndex = 0;
            int nnodes = faceOverlap.edges.size();
            for (int j1=1; j1<nnodes; j1++)
            {
              if ( nodesOverlap[faceOverlap[j1]] < nodesOverlap[faceOverlap[minIndex]] )
              {
                minIndex = j1;
              }
            }
#endif

            // Loop over all sub-triangles of this Overlap Face
            for ( int k = 0; k < nOverlapTriangles; k++ )
            {
#ifdef USE_MININDEX
                // Cornerpoints of triangle, they start at the minimal Node, for consistency
                const Node & node0 = nodesOverlap[faceOverlap[minIndex]];
                const Node & node1 = nodesOverlap[faceOverlap[(minIndex + k + 1)%nnodes]];
                const Node & node2 = nodesOverlap[faceOverlap[(minIndex + k + 2)%nnodes]];

                // Calculate the area of the modified Face
                Face faceTri ( 3 );
                faceTri.SetNode ( 0, faceOverlap[minIndex] );
                faceTri.SetNode ( 1, faceOverlap[(minIndex + k + 1)%nnodes] );
                faceTri.SetNode ( 2, faceOverlap[(minIndex + k + 2)%nnodes] );
#else
                // Cornerpoints of triangle
                const Node & node0 = nodesOverlap[faceOverlap[0]];
                const Node & node1 = nodesOverlap[faceOverlap[k+1]];
                const Node & node2 = nodesOverlap[faceOverlap[k+2]];

                // Calculate the area of the modified Face
                Face faceTri(3);
                faceTri.SetNode(0, faceOverlap[0]);
                faceTri.SetNode(1, faceOverlap[k+1]);
                faceTri.SetNode(2, faceOverlap[k+2]);
#endif

                double dTriangleArea =
                    CalculateFaceArea ( faceTri, nodesOverlap );

                // Coordinates of quadrature Node
                for ( int l = 0; l < TriQuadraturePoints; l++ )
                {
                    Node nodeQuadrature;
                    nodeQuadrature.x =
                        TriQuadratureG[l][0] * node0.x
                        + TriQuadratureG[l][1] * node1.x
                        + TriQuadratureG[l][2] * node2.x;

                    nodeQuadrature.y =
                        TriQuadratureG[l][0] * node0.y
                        + TriQuadratureG[l][1] * node1.y
                        + TriQuadratureG[l][2] * node2.y;

                    nodeQuadrature.z =
                        TriQuadratureG[l][0] * node0.z
                        + TriQuadratureG[l][1] * node1.z
                        + TriQuadratureG[l][2] * node2.z;

                    double dMag = sqrt (
                                      nodeQuadrature.x * nodeQuadrature.x
                                      + nodeQuadrature.y * nodeQuadrature.y
                                      + nodeQuadrature.z * nodeQuadrature.z );

                    nodeQuadrature.x /= dMag;
                    nodeQuadrature.y /= dMag;
                    nodeQuadrature.z /= dMag;

                    // Find components of quadrature point in basis
                    // of the first Face
                    double dAlpha;
                    double dBeta;

                    ApplyInverseMap (
                        faceFirst,
                        nodesFirst,
                        nodeQuadrature,
                        dAlpha,
                        dBeta );

                    // Check inverse map value
                    if ( ( dAlpha < -1.0e-13 ) || ( dAlpha > 1.0 + 1.0e-13 ) ||
                            ( dBeta  < -1.0e-13 ) || ( dBeta  > 1.0 + 1.0e-13 )
                       )
                    {
                        _EXCEPTION2 ( "Inverse Map out of range (%1.5e %1.5e)",
                                      dAlpha, dBeta );
                    }

                    // Sample the finite element at this point
                    SampleGLLFiniteElement (
                        nMonotoneType,
                        nP,
                        dAlpha,
                        dBeta,
                        dSampleCoeff );

                    // Add sample coefficients to the map
                    for ( int p = 0; p < nP; p++ )
                    {
                        for ( int q = 0; q < nP; q++ )
                        {
                            dRemapCoeff[p][q][j] +=
                                TriQuadratureW[l]
                                * dTriangleArea
                                * dSampleCoeff[p][q]
                                / m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
            }
        }

#ifdef VERBOSE
        output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        for ( int j = 0; j < nOverlapFaces; j++ )
        {
            for ( int p = 0; p < nP; p++ )
            {
                for ( int q = 0; q < nP; q++ )
                {
                    output_file << dRemapCoeff[p][q][j] << " ";
                }
            }
        }
        output_file << std::endl;
#endif

        // Force consistency and conservation
        if ( !fNoConservation )
        {
            double dTargetArea = 0.0;
            for ( int j = 0; j < nOverlapFaces; j++ )
            {
                dTargetArea += m_meshOverlap->vecFaceArea[ixOverlap + j];
            }

            for ( int p = 0; p < nP; p++ )
            {
                for ( int q = 0; q < nP; q++ )
                {
                    vecSourceArea[p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                }
            }

            const double areaTolerance = 1e-10;
            // Source elements are completely covered by target volumes
            if ( fabs ( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea ) <= areaTolerance )
            {
                vecTargetArea.Allocate ( nOverlapFaces );
                for ( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }

                dCoeff.Allocate ( nOverlapFaces, nP * nP );

                for ( int j = 0; j < nOverlapFaces; j++ )
                {
                    for ( int p = 0; p < nP; p++ )
                    {
                        for ( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }

                // Target volumes only partially cover source elements
            }
            else if ( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea > areaTolerance )
            {
                double dExtraneousArea = m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea;

                vecTargetArea.Allocate ( nOverlapFaces + 1 );
                for ( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }
                vecTargetArea[nOverlapFaces] = dExtraneousArea;

#ifdef VERBOSE
                Announce ( "Partial volume: %i (%1.10e / %1.10e)",
                           ixFirst, dTargetArea, m_meshInputCov->vecFaceArea[ixFirst] );
#endif
                if ( dTargetArea > m_meshInputCov->vecFaceArea[ixFirst] )
                {
                    _EXCEPTIONT ( "Partial element area exceeds total element area" );
                }

                dCoeff.Allocate ( nOverlapFaces + 1, nP * nP );

                for ( int j = 0; j < nOverlapFaces; j++ )
                {
                    for ( int p = 0; p < nP; p++ )
                    {
                        for ( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }
                for ( int p = 0; p < nP; p++ )
                {
                    for ( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] =
                            dataGLLJacobian[p][q][ixFirst];
                    }
                }
                for ( int j = 0; j < nOverlapFaces; j++ )
                {
                    for ( int p = 0; p < nP; p++ )
                    {
                        for ( int q = 0; q < nP; q++ )
                        {
                            dCoeff[nOverlapFaces][p * nP + q] -=
                                dRemapCoeff[p][q][j]
                                * m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
                for ( int p = 0; p < nP; p++ )
                {
                    for ( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] /= dExtraneousArea;
                    }
                }

                // Source elements only partially cover target volumes
            }
            else
            {
                Announce ( "Coverage area: %1.10e, and target element area: %1.10e)",
                           ixFirst, m_meshInputCov->vecFaceArea[ixFirst], dTargetArea );
                _EXCEPTIONT ( "Target grid must be a subset of source grid" );
            }

            ForceConsistencyConservation3 (
                vecSourceArea,
                vecTargetArea,
                dCoeff,
                ( nMonotoneType > 0 )
                /*, m_remapper->lid_to_gid_covsrc[ixFirst]*/ );

            for ( int j = 0; j < nOverlapFaces; j++ )
            {
                for ( int p = 0; p < nP; p++ )
                {
                    for ( int q = 0; q < nP; q++ )
                    {
                        dRemapCoeff[p][q][j] = dCoeff[j][p * nP + q];
                    }
                }
            }
        }

#ifdef VERBOSE
        // output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        // for ( int j = 0; j < nOverlapFaces; j++ )
        // {
        //     for ( int p = 0; p < nP; p++ )
        //     {
        //         for ( int q = 0; q < nP; q++ )
        //         {
        //             output_file << dRemapCoeff[p][q][j] << " ";
        //         }
        //     }
        // }
        // output_file << std::endl;
#endif

        // Put these remap coefficients into the SparseMatrix map
        for ( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];
            if (ixSecondFace < 0) // signal to not participate, because it is a ghost target
              continue; // do not do anything
            for ( int p = 0; p < nP; p++ )
            {
                for ( int q = 0; q < nP; q++ )
                {
                    if ( fContinuousIn )
                    {
                        int ixFirstNode = dataGLLNodes[p][q][ixFirst] - 1;

                        smatMap ( ixSecondFace, ixFirstNode ) +=
                            dRemapCoeff[p][q][j]
                            * m_meshOverlap->vecFaceArea[ixOverlap + j]
                            / m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                    else
                    {
                        int ixFirstNode = ixFirst * nP * nP + p * nP + q;

                        smatMap ( ixSecondFace, ixFirstNode ) +=
                            dRemapCoeff[p][q][j]
                            * m_meshOverlap->vecFaceArea[ixOverlap + j]
                            / m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                }
            }
        }
        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }
#ifdef VERBOSE
    output_file.flush(); // required here
    output_file.close();
#endif

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_MOAB (
    const DataArray3D<int> & dataGLLNodesIn,
    const DataArray3D<double> & dataGLLJacobianIn,
    const DataArray3D<int> & dataGLLNodesOut,
    const DataArray3D<double> & dataGLLJacobianOut,
    const DataArray1D<double> & dataNodalAreaOut,
    int nPin,
    int nPout,
    int nMonotoneType,
    bool fContinuousIn,
    bool fContinuousOut,
    bool fNoConservation
)
{
    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule ( 8 );

    const DataArray2D<double> & dG = triquadrule.GetG();
    const DataArray1D<double> & dW = triquadrule.GetW();

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix<double> & smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D<double> dSampleCoeffIn( nPin, nPin );

    DataArray2D<double> dSampleCoeffOut( nPout, nPout );

    // Announcemnets
    if ( is_root )
    {
        Announce ( "[moab::TempestOnlineMap::LinearRemapGLLtoGLL2_MOAB] Finite Element to Finite Element Projection" );
        Announce ( "Order of the input FE polynomial interpolant: %i", nPin );
        Announce ( "Order of the output FE polynomial interpolant: %i", nPout );
    }

    // Build the integration array for each element on m_meshOverlap
    DataArray3D<double> dGlobalIntArray(
        nPin * nPin,
        m_meshOverlap->faces.size(),
        nPout * nPout );

    // Number of overlap Faces per source Face
    DataArray1D<int> nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;

    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        size_t ixOverlapTemp = ixOverlap;
        for ( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];

            if ( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 )
            {
                break;
            }

            nAllOverlapFaces[ixFirst]++;
        }

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Geometric area of each output node
    DataArray2D<double> dGeometricOutputArea(
        m_meshOutput->faces.size(), nPout * nPout );

    // Area of each overlap element in the output basis
    DataArray2D<double> dOverlapOutputArea(
        m_meshOverlap->faces.size(), nPout * nPout );

    // Loop through all faces on m_meshInput
    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = (m_meshInputCov->faces.size()/10);
#endif

    if ( is_root )
        Announce ( "Building conservative distribution maps" );
    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if ( ixFirst % outputFrequency == 0 && is_root )
        {
            Announce ( "Element %i", ixFirst );
        }
#endif

        // Quantities from the First Mesh
        const Face & faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector & nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        if ( !nOverlapFaces ) continue;
        /*
                // Calculate total element Jacobian
                double dTotalJacobian = 0.0;
                for (int s = 0; s < nPin; s++) {
                for (int t = 0; t < nPin; t++) {
                    dTotalJacobian += dataGLLJacobianIn[s][t][ixFirst];
                }
                }
        */

        // Loop through all Overlap Faces
        for ( int i = 0; i < nOverlapFaces; i++ )
        {

            // Quantities from the overlap Mesh
            const Face & faceOverlap = m_meshOverlap->faces[ixOverlap + i];

            const NodeVector & nodesOverlap = m_meshOverlap->nodes;

            int nOverlapTriangles = faceOverlap.edges.size() - 2;

            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            const NodeVector & nodesSecond = m_meshOutput->nodes;

            const Face & faceSecond = m_meshOutput->faces[ixSecond];

            // Loop over all sub-triangles of this Overlap Face
            for ( int j = 0; j < nOverlapTriangles; j++ )
            {

                // Cornerpoints of triangle
                const Node & node0 = nodesOverlap[faceOverlap[0]];
                const Node & node1 = nodesOverlap[faceOverlap[j + 1]];
                const Node & node2 = nodesOverlap[faceOverlap[j + 2]];

                // Calculate the area of the modified Face
                Face faceTri ( 3 );
                faceTri.SetNode ( 0, faceOverlap[0] );
                faceTri.SetNode ( 1, faceOverlap[j + 1] );
                faceTri.SetNode ( 2, faceOverlap[j + 2] );

                double dTriArea =
                    CalculateFaceArea ( faceTri, nodesOverlap );

                for ( int k = 0; k < triquadrule.GetPoints(); k++ )
                {
                    // Get the nodal location of this point
										double dX[3];

										dX[0] = dG(k,0) * node0.x + dG(k,1) * node1.x + dG(k,2) * node2.x;
										dX[1] = dG(k,0) * node0.y + dG(k,1) * node1.y + dG(k,2) * node2.y;
										dX[2] = dG(k,0) * node0.z + dG(k,1) * node1.z + dG(k,2) * node2.z;

										double dMag =
											sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

										dX[0] /= dMag;
										dX[1] /= dMag;
										dX[2] /= dMag;

										Node nodeQuadrature(dX[0], dX[1], dX[2]);

										// Find the components of this quadrature point in the basis
										// of the first Face.
										double dAlphaIn;
										double dBetaIn;

										ApplyInverseMap(
											faceFirst,
											nodesFirst,
											nodeQuadrature,
											dAlphaIn,
											dBetaIn);

										// Find the components of this quadrature point in the basis
										// of the second Face.
										double dAlphaOut;
										double dBetaOut;

                    ApplyInverseMap (
                        faceSecond,
                        nodesSecond,
                        nodeQuadrature,
                        dAlphaOut,
                        dBetaOut );

                    /*
                                        // Check inverse map value
                                        if ((dAlphaIn < 0.0) || (dAlphaIn > 1.0) ||
                                            (dBetaIn  < 0.0) || (dBetaIn  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaIn, dBetaIn);
                                        }

                                        // Check inverse map value
                                        if ((dAlphaOut < 0.0) || (dAlphaOut > 1.0) ||
                                            (dBetaOut  < 0.0) || (dBetaOut  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaOut, dBetaOut);
                                        }
                    */
                    // Sample the First finite element at this point
                    SampleGLLFiniteElement (
                        nMonotoneType,
                        nPin,
                        dAlphaIn,
                        dBetaIn,
                        dSampleCoeffIn );

                    // Sample the Second finite element at this point
                    SampleGLLFiniteElement (
                        nMonotoneType,
                        nPout,
                        dAlphaOut,
                        dBetaOut,
                        dSampleCoeffOut );

                    // Overlap output area
                    for ( int s = 0; s < nPout; s++ )
                    {
                        for ( int t = 0; t < nPout; t++ )
                        {
                            double dNodeArea =
                                dSampleCoeffOut[s][t]
                                * dW[k]
                                * dTriArea;

                            dOverlapOutputArea[ixOverlap + i][s * nPout + t] +=
                                dNodeArea;

                            dGeometricOutputArea[ixSecond][s * nPout + t] +=
                                dNodeArea;
                        }
                    }

                    // Compute overlap integral
                    int ixp = 0;
                    for ( int p = 0; p < nPin; p++ )
                    {
                        for ( int q = 0; q < nPin; q++ )
                        {

                            int ixs = 0;
                            for ( int s = 0; s < nPout; s++ )
                            {
                                for ( int t = 0; t < nPout; t++ )
                                {

                                    // Sample the Second finite element at this point
                                    dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
                                        dSampleCoeffOut[s][t]
                                        * dSampleCoeffIn[p][q]
                                        * dW[k]
                                        * dTriArea;

                                    ixs++;
                                }
                            }

                            ixp++;
                        }
                    }
                }
            }
        }

        // Coefficients
        DataArray2D<double> dCoeff( nOverlapFaces * nPout * nPout, nPin * nPin );

        for ( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            int ixp = 0;
            for ( int p = 0; p < nPin; p++ )
            {
                for ( int q = 0; q < nPin; q++ )
                {

                    int ixs = 0;
                    for ( int s = 0; s < nPout; s++ )
                    {
                        for ( int t = 0; t < nPout; t++ )
                        {
                            dCoeff[i * nPout * nPout + ixs][ixp] =
                                dGlobalIntArray[ixp][ixOverlap + i][ixs]
                                / dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Source areas
        DataArray1D<double> vecSourceArea ( nPin * nPin );

        for ( int p = 0; p < nPin; p++ )
        {
            for ( int q = 0; q < nPin; q++ )
            {
                vecSourceArea[p * nPin + q] =
                    dataGLLJacobianIn[p][q][ixFirst];
            }
        }

        // Target areas
        DataArray1D<double> vecTargetArea ( nOverlapFaces * nPout * nPout );

        for ( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            int ixs = 0;
            for ( int s = 0; s < nPout; s++ )
            {
                for ( int t = 0; t < nPout; t++ )
                {

                    vecTargetArea[i * nPout * nPout + ixs] =
                        dOverlapOutputArea[ixOverlap + i][nPout * s + t];

                    ixs++;
                }
            }
        }

        // Force consistency and conservation
        if ( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation (
                vecSourceArea,
                vecTargetArea,
                dCoeff,
                ( nMonotoneType != 0 ) );
        }

        // Update global coefficients
        for ( int i = 0; i < nOverlapFaces; i++ )
        {
            int ixp = 0;
            for ( int p = 0; p < nPin; p++ )
            {
                for ( int q = 0; q < nPin; q++ )
                {
                    int ixs = 0;
                    for ( int s = 0; s < nPout; s++ )
                    {
                        for ( int t = 0; t < nPout; t++ )
                        {

                            dGlobalIntArray[ixp][ixOverlap + i][ixs] =
                                dCoeff[i * nPout * nPout + ixs][ixp]
                                * dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }
        /*
                // Check column sums (conservation)
                for (int i = 0; i < nPin * nPin; i++) {
                    double dColSum = 0.0;
                    for (int j = 0; j < nOverlapFaces * nPout * nPout; j++) {
                        dColSum += dCoeff[j][i] * vecTargetArea[j];
                    }
                    printf("Col %i: %1.15e\n", i, dColSum / vecSourceArea[i]);
                }

                // Check row sums (consistency)
                for (int j = 0; j < nOverlapFaces * nPout * nPout; j++) {
                    double dRowSum = 0.0;
                    for (int i = 0; i < nPin * nPin; i++) {
                        dRowSum += dCoeff[j][i];
                    }
                    printf("Row %i: %1.15e\n", j, dRowSum);
                }
                _EXCEPTION();
        */

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Build redistribution map within target element
    Announce ( "Building redistribution maps on target mesh" );
    DataArray1D<double> dRedistSourceArea ( nPout * nPout );
    DataArray1D<double> dRedistTargetArea ( nPout * nPout );
    std::vector< DataArray2D<double> > dRedistributionMaps;
    dRedistributionMaps.resize ( m_meshOutput->faces.size() );

    for ( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        dRedistributionMaps[ixSecond].Allocate ( nPout * nPout, nPout * nPout );

        for ( int i = 0; i < nPout * nPout; i++ )
        {
            dRedistributionMaps[ixSecond][i][i] = 1.0;
        }

        for ( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistSourceArea[s] =
                dGeometricOutputArea[ixSecond][s];
        }

        for ( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistTargetArea[s] =
                dataGLLJacobianOut[s / nPout][s % nPout][ixSecond];
        }

        if ( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation (
                dRedistSourceArea,
                dRedistTargetArea,
                dRedistributionMaps[ixSecond],
                ( nMonotoneType != 0 ) );

            for ( int s = 0; s < nPout * nPout; s++ )
            {
                for ( int t = 0; t < nPout * nPout; t++ )
                {
                    dRedistributionMaps[ixSecond][s][t] *=
                        dRedistTargetArea[s] / dRedistSourceArea[t];
                }
            }
        }
    }

    // Construct the total geometric area
    DataArray1D<double> dTotalGeometricArea ( dataNodalAreaOut.GetRows() );
    for ( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        for ( int s = 0; s < nPout; s++ )
        {
            for ( int t = 0; t < nPout; t++ )
            {
                dTotalGeometricArea[dataGLLNodesOut[s][t][ixSecond] - 1]
                += dGeometricOutputArea[ixSecond][s * nPout + t];
            }
        }
    }

    // Compose the integration operator with the output map
    ixOverlap = 0;

    Announce ( "Assembling map" );

    // Map from source DOFs to target DOFs with redistribution applied
    DataArray2D<double> dRedistributedOp (
        nPin * nPin, nPout * nPout );

    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if ( ixFirst % outputFrequency == 0 )
        {
            Announce ( "Element %i", ixFirst );
        }
#endif

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        // Put composed array into map
        for ( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];

            dRedistributedOp.Zero();
            for ( int p = 0; p < nPin * nPin; p++ )
            {
                for ( int s = 0; s < nPout * nPout; s++ )
                {
                    for ( int t = 0; t < nPout * nPout; t++ )
                    {
                        dRedistributedOp[p][s] +=
                            dRedistributionMaps[ixSecondFace][s][t]
                            * dGlobalIntArray[p][ixOverlap + j][t];
                    }
                }
            }

            int ixp = 0;
            for ( int p = 0; p < nPin; p++ )
            {
                for ( int q = 0; q < nPin; q++ )
                {

                    int ixFirstNode;
                    if ( fContinuousIn )
                    {
                        ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                    }
                    else
                    {
                        ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
                    }

                    int ixs = 0;
                    for ( int s = 0; s < nPout; s++ )
                    {
                        for ( int t = 0; t < nPout; t++ )
                        {

                            int ixSecondNode;
                            if ( fContinuousOut )
                            {
                                ixSecondNode = dataGLLNodesOut[s][t][ixSecondFace] - 1;

                                if ( !fNoConservation )
                                {
                                    smatMap ( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs]
                                        / dataNodalAreaOut[ixSecondNode];
                                }
                                else
                                {
                                    smatMap ( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs]
                                        / dTotalGeometricArea[ixSecondNode];
                                }

                            }
                            else
                            {
                                ixSecondNode =
                                    ixSecondFace * nPout * nPout + s * nPout + t;

                                if ( !fNoConservation )
                                {
                                    smatMap ( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs]
                                        / dataGLLJacobianOut[s][t][ixSecondFace];
                                }
                                else
                                {
                                    smatMap ( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs]
                                        / dGeometricOutputArea[ixSecondFace][s * nPout + t];
                                }
                            }

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_Pointwise_MOAB (
    const DataArray3D<int> & dataGLLNodesIn,
    const DataArray3D<double> & /*dataGLLJacobianIn*/,
    const DataArray3D<int> & dataGLLNodesOut,
    const DataArray3D<double> & /*dataGLLJacobianOut*/,
    const DataArray1D<double> & dataNodalAreaOut,
    int nPin,
    int nPout,
    int nMonotoneType,
    bool fContinuousIn,
    bool fContinuousOut
)
{
    // Gauss-Lobatto quadrature within Faces
    DataArray1D<double> dGL;
    DataArray1D<double> dWL;

    GaussLobattoQuadrature::GetPoints ( nPout, 0.0, 1.0, dGL, dWL );

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix<double> & smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D<double> dSampleCoeffIn( nPin, nPin );

    // Announcemnets
    if ( is_root )
    {
        Announce ( "[moab::TempestOnlineMap::LinearRemapGLLtoGLL2_Pointwise_MOAB] Finite Element to Finite Element (Pointwise) Projection" );
        Announce ( "Order of the input FE polynomial interpolant: %i", nPin );
        Announce ( "Order of the output FE polynomial interpolant: %i", nPout );
    }

    // Number of overlap Faces per source Face
    DataArray1D<int> nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;

    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        size_t ixOverlapTemp = ixOverlap;
        for ( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];

            if ( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 )
                break;

            nAllOverlapFaces[ixFirst]++;
        }

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Number of times this point was found
    DataArray1D<bool> fSecondNodeFound ( dataNodalAreaOut.GetRows() );

    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = (m_meshInputCov->faces.size()/10);
#endif

    // Loop through all faces on m_meshInputCov
    for ( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if ( ixFirst % outputFrequency == 0 )
        {
            Announce ( "Element %i", ixFirst );
        }
#endif

        // Quantities from the First Mesh
        const Face & faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector & nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        // Loop through all Overlap Faces
        for ( int i = 0; i < nOverlapFaces; i++ )
        {
            // Quantities from the overlap Mesh
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlap + i];

            // const NodeVector & nodesOverlap = m_meshOverlap->nodes;

            // size_t nOverlapTriangles = faceOverlap.edges.size() - 2;

            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            const NodeVector & nodesSecond = m_meshOutput->nodes;

            const Face & faceSecond = m_meshOutput->faces[ixSecond];

            // Loop through all nodes on the second face
            for ( int s = 0; s < nPout; s++ )
            {
                for ( int t = 0; t < nPout; t++ )
                {
                    size_t ixSecondNode;
                    if ( fContinuousOut )
                    {
                        ixSecondNode = dataGLLNodesOut[s][t][ixSecond] - 1;
                    }
                    else
                    {
                        ixSecondNode =
                            ixSecond * nPout * nPout + s * nPout + t;
                    }

                    if ( ixSecondNode >= fSecondNodeFound.GetRows() )
                    {
                        _EXCEPTIONT ( "Logic error" );
                    }

                    // Check if this node has been found already
                    if ( fSecondNodeFound[ixSecondNode] )
                    {
                        continue;
                    }

                    // Check this node
                    Node node;
                    Node dDx1G;
                    Node dDx2G;

                    ApplyLocalMap (
                        faceSecond,
                        nodesSecond,
                        dGL[t],
                        dGL[s],
                        node,
                        dDx1G,
                        dDx2G );

                    // Find the components of this quadrature point in the basis
                    // of the first Face.
                    double dAlphaIn;
                    double dBetaIn;

                    ApplyInverseMap (
                        faceFirst,
                        nodesFirst,
                        node,
                        dAlphaIn,
                        dBetaIn );

                    // Check if this node is within the first Face
                    if ( ( dAlphaIn < -1.0e-10 ) || ( dAlphaIn > 1.0 + 1.0e-10 ) ||
                            ( dBetaIn  < -1.0e-10 ) || ( dBetaIn  > 1.0 + 1.0e-10 )
                       )
                    {
                        continue;
                    }

                    // Node is within the overlap region, mark as found
                    fSecondNodeFound[ixSecondNode] = true;

                    // Sample the First finite element at this point
                    SampleGLLFiniteElement (
                        nMonotoneType,
                        nPin,
                        dAlphaIn,
                        dBetaIn,
                        dSampleCoeffIn );

                    // Add to map
                    for ( int p = 0; p < nPin; p++ )
                    {
                        for ( int q = 0; q < nPin; q++ )
                        {
                            int ixFirstNode;
                            if ( fContinuousIn )
                            {
                                ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                            }
                            else
                            {
                                ixFirstNode =
                                    ixFirst * nPin * nPin + p * nPin + q;
                            }

                            smatMap ( ixSecondNode, ixFirstNode ) +=
                                dSampleCoeffIn[p][q];
                        }
                    }
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Check for missing samples
    for ( size_t i = 0; i < fSecondNodeFound.GetRows(); i++ )
    {
        if ( !fSecondNodeFound[i] )
        {
            _EXCEPTION1 ( "Can't sample point %i", i );
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

