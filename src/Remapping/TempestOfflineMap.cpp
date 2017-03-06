///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOfflineMap.cpp
///	\author  Vijay Mahadevan
///	\version Feb 28, 2017
///

#include "Announce.h"
#include "DataMatrix3D.h"
#include "FiniteElementTools.h"
#include "LinearRemapFV.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"

#include "TempestOfflineMap.hpp"
#include "DebugOutput.hpp"

#include "netcdfcpp.h"

#include <fstream>
#include <cmath>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

moab::TempestOfflineMap::TempestOfflineMap(moab::TempestRemapper* remapper) : OfflineMap(), m_remapper(remapper)
{
    // Get the references for the MOAB core objects
    mbCore = m_remapper->get_interface();
    pcomm = m_remapper->get_parallel_communicator();

    // Update the references to the meshes
    m_meshInput = remapper->GetMesh(moab::Remapper::SourceMesh);
    m_meshInputCov = remapper->GetCoveringMesh();
    m_meshOutput = remapper->GetMesh(moab::Remapper::TargetMesh);
    m_meshOverlap = remapper->GetMesh(moab::Remapper::IntersectedMesh);

    moab::DebugOutput dbgprint(std::cout, (pcomm ? pcomm->rank() : 0));

    // Compute and store the total number of source and target DoFs corresponding
    // to number of rows and columns in the mapping.
    
    // Initialize dimension information from file
    dbgprint.printf(0, "Initializing dimensions of map\n");
    dbgprint.printf(0, "Input mesh\n");
    this->InitializeSourceDimensionsFromMesh(*m_meshInputCov);
    dbgprint.printf(0, "Output mesh\n");
    this->InitializeTargetDimensionsFromMesh(*m_meshOutput);
    // dbgprint.printf(0, "----------------------------------\n");

    // Build a matrix of source and target discretization so that we know how to assign
    // the global DoFs in parallel for the mapping weights
    // For example, FV->FV: rows X cols = faces_source X faces_target
    // 
    // MPI_Scan()
}

///////////////////////////////////////////////////////////////////////////////

moab::TempestOfflineMap::~TempestOfflineMap()
{
    mbCore = NULL;
    pcomm = NULL;
    m_meshInput = NULL;
    m_meshOutput = NULL;
    m_meshOverlap = NULL;
}

///////////////////////////////////////////////////////////////////////////////

static void ParseVariableList(
	const std::string & strVariables,
	std::vector< std::string > & vecVariableStrings
) {
	unsigned iVarBegin = 0;
	unsigned iVarCurrent = 0;

	// Parse variable name
	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
			(strVariables[iVarCurrent] == ',') ||
			(strVariables[iVarCurrent] == ' ')
		) {
			if (iVarCurrent == iVarBegin) {
				if (iVarCurrent >= strVariables.length()) {
					break;
				}
				continue;
			}

			vecVariableStrings.push_back(
				strVariables.substr(iVarBegin, iVarCurrent - iVarBegin));

			iVarBegin = iVarCurrent + 1;
		}

		iVarCurrent++;
	}
}


///////////////////////////////////////////////////////////////////////////////

extern void BuildIntegrationArray(
    const Mesh & meshInput,
    const Mesh & meshOverlap,
    const TriangularQuadratureRule & triquadrule,
    int ixFirstFace,
    int ixOverlapBegin,
    int ixOverlapEnd,
    int nOrder,
    DataMatrix<double> & dIntArray
);

extern void InvertFitArray_Corrected(
    const DataVector<double> & dConstraint,
    DataMatrix<double> & dFitArray,
    DataVector<double> & dFitWeights,
    DataMatrix<double> & dFitArrayPlus
);

/// <summary>
///     Face index and distance metric pair.
/// </summary>
typedef std::pair<int, int> FaceDistancePair;

/// <summary>
///     Vector storing adjacent Faces.
/// </summary>
typedef std::vector<FaceDistancePair> AdjacentFaceVector;

extern void BuildFitArray(
    const Mesh & mesh,
    const TriangularQuadratureRule & triquadrule,
    int ixFirst,
    const AdjacentFaceVector & vecAdjFaces,
    int nOrder,
    int nFitWeightsExponent,
    const DataVector<double> & dConstraint,
    DataMatrix<double> & dFitArray,
    DataVector<double> & dFitWeights
);

extern void GetAdjacentFaceVectorByEdge(
    const Mesh & mesh,
    int iFaceInitial,
    int nRequiredFaceSetSize,
    AdjacentFaceVector & vecFaces
);

void moab::TempestOfflineMap::LinearRemapFVtoFV_Tempest_MOAB(
    int nOrder
) {
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if (m_meshInputCov->revnodearray.size() == 0) {
        _EXCEPTIONT("ReverseNodeArray has not been calculated for m_meshInput");
    }

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

    // Number of elements needed
    const int nCoefficients = nOrder * (nOrder + 1) / 2;

// #pragma message "This should be a command-line parameter"
    // Number of faces you need
    const int nRequiredFaceSetSize = nCoefficients;

    // Fit weight exponent
    const int nFitWeightsExponent = nOrder + 2;

    // Announcemnets
    if (!pcomm->rank()) {
        Announce("Triangular quadrature rule order %i", TriQuadRuleOrder);
        Announce("Number of coefficients: %i", nCoefficients);
        Announce("Required adjacency set size: %i", nRequiredFaceSetSize);
        Announce("Fit weights exponent: %i", nFitWeightsExponent);
    }
    if (!pcomm->rank()) Announce("Rank 0 - Input size: %i ", m_meshInputCov->faces.size());
    if (pcomm->rank()) Announce("Rank 1 - Input size: %i", m_meshInputCov->faces.size());
    // const bool debug = (!pcomm->rank());

    // Current overlap face
    int ixOverlap = 0;

    // Loop through all faces on m_meshInput
    for (size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++) {

        // Output every 100 elements
        if (ixFirst % 100 == 0) {
            Announce("Element %i/%i", ixFirst, m_meshInputCov->faces.size());
        }

        std::vector<std::pair<int,int> > ixFaces;
        for (unsigned ix = 0; ix < m_meshOverlap->faces.size(); ix++) {
            if (ixFirst - m_meshOverlap->vecSourceFaceIx[ix] == 0) {
                // if (!pcomm->rank()) std::cout << "[" << ix << "] SrcFaceIx = " << m_meshOverlap->vecSourceFaceIx[ix] << std::endl;
                unsigned jx;
                for (jx = 1; jx < m_meshOverlap->faces.size()-ix; jx++) {
                    // if (!pcomm->rank()) std::cout << "\t[" << jx << "] SrcFaceJx = " << m_meshOverlap->vecSourceFaceIx[jx] << std::endl;
                    if (ixFirst - m_meshOverlap->vecSourceFaceIx[ix+jx] != 0)
                        break;
                }
                ixFaces.push_back(std::make_pair<int,int>(ix,ix+jx));
                ix += jx-1;
            }
        }

        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so on wrt the input mesh data.
        // Then the overlap_end and overlap_begin will be correct. However, the relation with MOAB and Tempest will go out of the roof.

        int gnOverlapFaces = 0;
        for (unsigned ifaceIndx = 0; ifaceIndx < ixFaces.size(); ++ifaceIndx) {

            int ixOverlapBegin = ixFaces[ifaceIndx].first;
            int ixOverlapEnd = ixFaces[ifaceIndx].second;
            int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

            // Build integration array
            DataMatrix<double> dIntArray;

            BuildIntegrationArray(
                *m_meshInputCov,
                *m_meshOverlap,
                triquadrule,
                ixFirst,
                ixOverlapBegin,
                ixOverlapEnd,
                nOrder,
                dIntArray);

            // Set of Faces to use in building the reconstruction and associated
            // distance metric.
            AdjacentFaceVector vecAdjFaces;

            GetAdjacentFaceVectorByEdge(
                *m_meshInputCov,
                ixFirst,
                nRequiredFaceSetSize,
                vecAdjFaces);

            // Number of adjacent Faces
            int nAdjFaces = vecAdjFaces.size();

            // Determine the conservative constraint equation
            DataVector<double> dConstraint;

            dConstraint.Initialize(nCoefficients);

            double dFirstArea = m_meshInputCov->vecFaceArea[ixFirst];

            for (int p = 0; p < nCoefficients; p++) {
                for (int j = 0; j < nOverlapFaces; j++) {
                    dConstraint[p] += dIntArray[p][j];
                }
                dConstraint[p] /= dFirstArea;
            }

            // Build the fit array from the integration operator
            DataMatrix<double> dFitArray;
            DataVector<double> dFitWeights;
            DataMatrix<double> dFitArrayPlus;

            BuildFitArray(
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
            InvertFitArray_Corrected(
                dConstraint,
                dFitArray,
                dFitWeights,
                dFitArrayPlus
            );
        
            // if (debug) Announce("[%i] (%d,%d) vecSourceFaceIx: %i and vecTargetFaceIx: %i -- Area: %12.10f", ixFirst, nAdjFaces, nOverlapFaces, m_meshOverlap->vecSourceFaceIx[ixOverlapBegin], m_meshOverlap->vecTargetFaceIx[ixOverlapEnd], dFirstArea);

            // Multiply integration array and fit array
            DataMatrix<double> dComposedArray;
            dComposedArray.Initialize(nAdjFaces, nOverlapFaces);

            for (int i = 0; i < nAdjFaces; i++) {
                for (int j = 0; j < nOverlapFaces; j++) {
                    for (int k = 0; k < nCoefficients; k++) {
                        dComposedArray[i][j] += dIntArray[k][j] * dFitArrayPlus[i][k];
                    }
                }
            }

/*
        for (int j = 0; j < nOverlapFaces; j++) {
            dComposedArray[0][j] = meshOverlap.vecFaceArea[ixOverlap + j];
        }

        for (int i = 0; i < nAdjFaces; i++) {
        for (int j = 0; j < nOverlapFaces; j++) {
        for (int k = 1; k < nCoefficients; k++) {
            double dOverlapArea =
                meshOverlap.vecFaceArea[ixOverlap + j];

            dComposedArray[i][j] +=
                (dIntArray[k][j] - dConstraint[k] * dOverlapArea)
                    * dFitArrayPlus[i][k];
        }
        }
        }
*/

            // Put composed array into map
            for (unsigned i = 0; i < vecAdjFaces.size(); i++) {
            for (int j = 0; j < nOverlapFaces; j++) {
                int& ixFirstFaceLoc = vecAdjFaces[i].first;
                int& ixSecondFaceLoc = m_meshOverlap->vecTargetFaceIx[ixOverlapBegin + j];
                // int ixFirstFaceGlob = m_remapper->GetGlobalID(moab::Remapper::SourceMesh, ixFirstFaceLoc);
                // int ixSecondFaceGlob = m_remapper->GetGlobalID(moab::Remapper::TargetMesh, ixSecondFaceLoc);

                // if (debug) std::cout << "!"<<gnOverlapFaces<<"! (" << ixSecondFaceGlob << ", " << ixFirstFaceGlob << ") = " << dComposedArray[i][j]/m_meshOutput->vecFaceArea[ixSecondFaceLoc] << "\t";

                m_mapRemap(ixSecondFaceLoc, ixFirstFaceLoc) +=
                    dComposedArray[i][j]
                    / m_meshOutput->vecFaceArea[ixSecondFaceLoc];
            }
            }
            // if (debug) std::cout << "\n";

            gnOverlapFaces += nOverlapFaces;

        }

        // Increment the current overlap index
        ixOverlap += gnOverlapFaces;
    }
}


///////////////////////////////////////////////////////////////////////////////

extern void ForceConsistencyConservation3(
    const DataVector<double> & vecSourceArea,
    const DataVector<double> & vecTargetArea,
    DataMatrix<double> & dCoeff,
    bool fMonotone
);

void moab::TempestOfflineMap::LinearRemapSE4_Tempest_MOAB(
    const DataMatrix3D<int> & dataGLLNodes,
    const DataMatrix3D<double> & dataGLLJacobian,
    int nMonotoneType,
    bool fContinuousIn,
    bool fNoConservation
) {
    // Order of the polynomial interpolant
    int nP = dataGLLNodes.GetRows();

    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4; 

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

    int TriQuadraturePoints = triquadrule.GetPoints();

    const DataMatrix<double> & TriQuadratureG = triquadrule.GetG();

    const DataVector<double> & TriQuadratureW = triquadrule.GetW();

    // Sample coefficients
    DataMatrix<double> dSampleCoeff;
    dSampleCoeff.Initialize(nP, nP);

    // GLL Quadrature nodes on quadrilateral elements
    DataVector<double> dG;
    DataVector<double> dW;
    GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix<double> & smatMap = this->GetSparseMatrix();

    // NodeVector from meshOverlap
    const NodeVector & nodesOverlap = m_meshOverlap->nodes;
    const NodeVector & nodesFirst   = m_meshInputCov->nodes;

    // Vector of source areas
    DataVector<double> vecSourceArea;
    vecSourceArea.Initialize(nP * nP);

    DataVector<double> vecTargetArea;

    DataMatrix<double> dCoeff;

    // Current Overlap Face
    int ixOverlap = 0;

    // Loop over all input Faces
    for (size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++) {

        const Face & faceFirst = m_meshInputCov->faces[ixFirst];

        if (faceFirst.edges.size() != 4) {
            _EXCEPTIONT("Only quadrilateral elements allowed for SE remapping");
        }

        // Output every 100 elements
        if (ixFirst % 100 == 0) {
            Announce("Element %i/%i", ixFirst, m_meshInputCov->faces.size());
        }

        std::vector<std::pair<int,int> > ixFaces;
        for (unsigned ix = 0; ix < m_meshOverlap->faces.size(); ix++) {
            if (ixFirst - m_meshOverlap->vecSourceFaceIx[ix] == 0) { 
                // if (!pcomm->rank()) std::cout << "[" << ix << "] SrcFaceIx = " << m_meshOverlap->vecSourceFaceIx[ix] << std::endl;
                unsigned jx;
                int ntotTris=0;
                for (jx = 1; jx < m_meshOverlap->faces.size()-ix; jx++) {
                    // if (!pcomm->rank()) std::cout << "\t[" << jx << "] SrcFaceJx = " << m_meshOverlap->vecSourceFaceIx[jx] << std::endl;
                    if (ixFirst - m_meshOverlap->vecSourceFaceIx[ix+jx] != 0)
                        break;

                    ntotTris += m_meshOverlap->faces[ix+jx].edges.size() - 2;
                }
                ixFaces.push_back(std::make_pair<int,int>(ix,ix+jx));
                ix += jx-1;
            }
        }

        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so on wrt the input mesh data
        // Then the overlap_end and overlap_begin will be correct. However, the relation with MOAB and Tempest will go out of the roof

        // // Determine how many overlap Faces and triangles are present
        // int ixOverlapTemp = ixOverlap;
        // for (; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++) {

        //     const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];

        //     if (m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
        //         break;
        //     }

        //     nOverlapFaces++;
        //     nTotalOverlapTriangles += faceOverlap.edges.size() - 2;
        // }

        int gnOverlapFaces = 0;
        double globTargetArea = 0.0;

        for (unsigned ifaceIndx = 0; ifaceIndx < ixFaces.size(); ++ifaceIndx) {

            // Number of overlapping Faces and triangles
            int ixOverlapBegin = ixFaces[ifaceIndx].first;
            int ixOverlapEnd = ixFaces[ifaceIndx].second;
            int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

            // No overlaps
            if (nOverlapFaces == 0) {
                continue;
            }

            // Allocate remap coefficients array for meshFirst Face
            DataMatrix3D<double> dRemapCoeff;
            dRemapCoeff.Initialize(nP, nP, nOverlapFaces);

            // Find the local remap coefficients
            for (int j = 0; j < nOverlapFaces; j++) {
                const Face & faceOverlap = m_meshOverlap->faces[ixOverlapBegin + j];

                int nOverlapTriangles = faceOverlap.edges.size() - 2;

                // Loop over all sub-triangles of this Overlap Face
                for (int k = 0; k < nOverlapTriangles; k++) {

                    // Cornerpoints of triangle
                    const Node & node0 = nodesOverlap[faceOverlap[0]];
                    const Node & node1 = nodesOverlap[faceOverlap[k+1]];
                    const Node & node2 = nodesOverlap[faceOverlap[k+2]];

                    // Calculate the area of the modified Face
                    Face faceTri(3);
                    faceTri.SetNode(0, faceOverlap[0]);
                    faceTri.SetNode(1, faceOverlap[k+1]);
                    faceTri.SetNode(2, faceOverlap[k+2]);

                    double dTriangleArea =
                        CalculateFaceArea(faceTri, nodesOverlap);

                    // Coordinates of quadrature Node
                    for (int l = 0; l < TriQuadraturePoints; l++) {
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

                        double dMag = sqrt(
                              nodeQuadrature.x * nodeQuadrature.x
                            + nodeQuadrature.y * nodeQuadrature.y
                            + nodeQuadrature.z * nodeQuadrature.z);

                        nodeQuadrature.x /= dMag;
                        nodeQuadrature.y /= dMag;
                        nodeQuadrature.z /= dMag;

                        // Find components of quadrature point in basis
                        // of the first Face
                        double dAlpha;
                        double dBeta;

                        ApplyInverseMap(
                            faceFirst,
                            nodesFirst,
                            nodeQuadrature,
                            dAlpha,
                            dBeta);

                        // Check inverse map value
                        if ((dAlpha < -1.0e-13) || (dAlpha > 1.0 + 1.0e-13) ||
                            (dBeta  < -1.0e-13) || (dBeta  > 1.0 + 1.0e-13)
                        ) {
                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                dAlpha, dBeta);
                        }

                        // Sample the finite element at this point
                        SampleGLLFiniteElement(
                            nMonotoneType,
                            nP,
                            dAlpha,
                            dBeta,
                            dSampleCoeff);

                        // Add sample coefficients to the map
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nP; q++) {

                                dRemapCoeff[p][q][j] +=
                                    TriQuadratureW[l]
                                    * dTriangleArea
                                    * dSampleCoeff[p][q]
                                    / m_meshOverlap->vecFaceArea[ixOverlapBegin + j];
                            }
                        }
                    }
                }
            }

            // Force consistency and conservation
            if (!fNoConservation) {
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nP; q++) {
                        vecSourceArea[p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                    }
                }

                double dTargetArea = 0.0;
                vecTargetArea.Initialize(nOverlapFaces);
                for (int j = 0; j < nOverlapFaces; j++) {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlapBegin + j];
                    dTargetArea += m_meshOverlap->vecFaceArea[ixOverlapBegin + j];
                }
                globTargetArea += dTargetArea;

                // const double factor = dTargetArea/m_meshInputCov->vecFaceArea[ixFirst]*0+1.0/ixFaces.size();
                // for (int j = 0; j < nOverlapFaces; j++) {
                //     vecTargetArea[j] *= factor;
                // }

                if (fabs(globTargetArea - m_meshInputCov->vecFaceArea[ixFirst]) > 1.0e-10) {
                    Announce("TempestOfflineMap: Partial element: %i, areas = %f percent", ixFirst, 100*globTargetArea/m_meshInputCov->vecFaceArea[ixFirst]);
                }
                else
                {
                    dCoeff.Initialize(nOverlapFaces, nP * nP);

                    for (int j = 0; j < nOverlapFaces; j++) {
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nP; q++)
                                dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }

                    ForceConsistencyConservation3(
                        vecSourceArea,
                        vecTargetArea,
                        dCoeff,
                        (nMonotoneType != 0));

                    for (int j = 0; j < nOverlapFaces; j++) {
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nP; q++)
                                dRemapCoeff[p][q][j] = dCoeff[j][p * nP + q];
                        }
                    }
                }
            }

            // Put these remap coefficients into the SparseMatrix map
            for (int j = 0; j < nOverlapFaces; j++) {
                int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlapBegin + j];

                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nP; q++) {

                        if (fContinuousIn) {
                            int ixFirstNode = dataGLLNodes[p][q][ixFirst] - 1;

                            smatMap(ixSecondFace, ixFirstNode) +=
                                dRemapCoeff[p][q][j]
                                * m_meshOverlap->vecFaceArea[ixOverlapBegin + j]
                                / m_meshOutput->vecFaceArea[ixSecondFace];

                        } else {
                            int ixFirstNode = ixFirst * nP * nP + p * nP + q;

                            smatMap(ixSecondFace, ixFirstNode) +=
                                dRemapCoeff[p][q][j]
                                * m_meshOverlap->vecFaceArea[ixOverlapBegin + j]
                                / m_meshOutput->vecFaceArea[ixSecondFace];
                        }
                    }
                }
            }
            gnOverlapFaces += nOverlapFaces;

        }
        // Increment the current overlap index
        ixOverlap += gnOverlapFaces;
    }
}

///////////////////////////////////////////////////////////////////////////////
moab::ErrorCode moab::TempestOfflineMap::GenerateOfflineMap( std::string strInputType, std::string strOutputType,
                                            int nPin, int nPout,
                                            bool fBubble, int fMonotoneTypeID,
                                            bool fVolumetric, bool fNoConservation, bool fNoCheck,
                                            std::string strVariables, std::string strOutputMap,
                                            std::string strInputData, std::string strOutputData,
                                            std::string strNColName, bool fOutputDouble,
                                            std::string strPreserveVariables, bool fPreserveAll, double dFillValueOverride )
{
    NcError error(NcError::silent_nonfatal);

    moab::DebugOutput dbgprint(std::cout, (pcomm ? pcomm->rank() : 0));

try {

    // Input / Output types
    enum DiscretizationType {
        DiscretizationType_FV,
        DiscretizationType_CGLL,
        DiscretizationType_DGLL
    };

    // Check command line parameters (data arguments)
    if ((strInputData != "") && (strOutputData == "")) {
        _EXCEPTIONT("--in_data specified without --out_data");
    }
    if ((strInputData == "") && (strOutputData != "")) {
        _EXCEPTIONT("--out_data specified without --in_data");
    }

    // Check command line parameters (data type arguments)
    STLStringHelper::ToLower(strInputType);
    STLStringHelper::ToLower(strOutputType);

    DiscretizationType eInputType;
    DiscretizationType eOutputType;

    if (strInputType == "fv") {
        eInputType = DiscretizationType_FV;
    } else if (strInputType == "cgll") {
        eInputType = DiscretizationType_CGLL;
    } else if (strInputType == "dgll") {
        eInputType = DiscretizationType_DGLL;
    } else {
        _EXCEPTION1("Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
            strInputType.c_str());
    }

    if (strOutputType == "fv") {
        eOutputType = DiscretizationType_FV;
    } else if (strOutputType == "cgll") {
        eOutputType = DiscretizationType_CGLL;
    } else if (strOutputType == "dgll") {
        eOutputType = DiscretizationType_DGLL;
    } else {
        _EXCEPTION1("Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
            strOutputType.c_str());
    }

    // Monotonicity flags
    int nMonotoneType = fMonotoneTypeID;

/*
    // Volumetric
    if (fVolumetric && (nMonotoneType != 0)) {
        _EXCEPTIONT("--volumetric cannot be used in conjunction with --mono#");
    }
*/

    // // // Initialize dimension information from file
    // // if (created_local)
    // {
    // 	dbgprint.printf(0, "Initializing dimensions of map\n");
	   //  dbgprint.printf(0, "Input mesh\n");
	   //  this->InitializeSourceDimensionsFromMesh(m_meshInput);
	   //  dbgprint.printf(0, "Output mesh\n");
	   //  this->InitializeTargetDimensionsFromMesh(meshOutput);
	   //  // dbgprint.printf(0, "----------------------------------\n");
    // }

    // Parse variable list
    std::vector< std::string > vecVariableStrings;
    ParseVariableList(strVariables, vecVariableStrings);

    // Parse preserve variable list
    std::vector< std::string > vecPreserveVariableStrings;
    ParseVariableList(strPreserveVariables, vecPreserveVariableStrings);

    if (fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
        _EXCEPTIONT("--preserveall and --preserve cannot both be specified");
    }

    // Calculate Face areas
    if (!pcomm->rank()) dbgprint.printf(0, "Calculating input mesh Face areas\n");
    double dTotalAreaInput_loc = m_meshInput->CalculateFaceAreas();
    Real dTotalAreaInput;
    MPI_Allreduce(&dTotalAreaInput_loc, &dTotalAreaInput, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
    if (!pcomm->rank()) dbgprint.printf(0, "Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput);

    double dTotalAreaInputCov_loc = m_meshInputCov->CalculateFaceAreas();
    Real dTotalAreaInputCov;
    MPI_Allreduce(&dTotalAreaInputCov_loc, &dTotalAreaInputCov, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
    if (!pcomm->rank()) dbgprint.printf(0, "Input Mesh (coverage set) Geometric Area: %1.15e\n", dTotalAreaInputCov);

    // Input mesh areas
    if (eInputType == DiscretizationType_FV) {
        this->SetSourceAreas(m_meshInputCov->vecFaceArea);
    }

    // Calculate Face areas
    if (!pcomm->rank()) dbgprint.printf(0, "Calculating output mesh Face areas\n");
    Real dTotalAreaOutput_loc = m_meshOutput->CalculateFaceAreas();
    Real dTotalAreaOutput;
    MPI_Allreduce(&dTotalAreaOutput_loc, &dTotalAreaOutput, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
    if (!pcomm->rank()) dbgprint.printf(0, "Output Mesh Geometric Area: %1.15e\n", dTotalAreaOutput);

    // Output mesh areas
    if (eOutputType == DiscretizationType_FV) {
        this->SetTargetAreas(m_meshOutput->vecFaceArea);
    }

    // Verify that overlap mesh is in the correct order
    int ixSourceFaceMax = (-1);
    int ixTargetFaceMax = (-1);

    if (m_meshOverlap->vecSourceFaceIx.size() !=
        m_meshOverlap->vecTargetFaceIx.size()
    ) {
        _EXCEPTIONT("Invalid overlap mesh:\n"
            "    Possible mesh file corruption?");
    }

    for (unsigned i = 0; i < m_meshOverlap->vecSourceFaceIx.size(); i++) {
        if (m_meshOverlap->vecSourceFaceIx[i] + 1 > ixSourceFaceMax) {
            ixSourceFaceMax = m_meshOverlap->vecSourceFaceIx[i] + 1;
        }
        if (m_meshOverlap->vecTargetFaceIx[i] + 1 > ixTargetFaceMax) {
            ixTargetFaceMax = m_meshOverlap->vecTargetFaceIx[i] + 1;
        }
    }

    dbgprint.printf(0, "m_meshInput->faces = %lu, m_meshInputCov->faces = %lu, m_meshOutput->faces = %lu, ixSourceFaceMax = %d\n", m_meshInput->faces.size(), m_meshInputCov->faces.size(), m_meshOutput->faces.size(), ixSourceFaceMax);

    // Check for forward correspondence in overlap mesh
    if (// m_meshInputCov->faces.size() - ixSourceFaceMax == 0 //&&
        (m_meshOutput->faces.size() - ixTargetFaceMax == 0)
    ) {
        if (!pcomm->rank()) dbgprint.printf(0, "Overlap mesh forward correspondence found\n");

    // Check for reverse correspondence in overlap mesh
    } else if (
        // m_meshOutput->faces.size() - ixSourceFaceMax == 0 //&&
        (m_meshInputCov->faces.size() - ixTargetFaceMax == 0)
    ) {
        if (!pcomm->rank()) dbgprint.printf(0, "Overlap mesh reverse correspondence found (reversing)\n");

        // Reorder overlap mesh
        m_meshOverlap->ExchangeFirstAndSecondMesh();

    // No correspondence found
    } else {
        _EXCEPTION4("Invalid overlap mesh:\n"
            "    No correspondence found with input and output meshes (%i,%i) vs (%i,%i)",
            m_meshInputCov->faces.size(), m_meshOutput->faces.size(), ixSourceFaceMax, ixTargetFaceMax);
    }

    // Calculate Face areas
    if (!pcomm->rank()) dbgprint.printf(0, "Calculating overlap mesh Face areas\n");
    Real dTotalAreaOverlap_loc = m_meshOverlap->CalculateFaceAreas();
    Real dTotalAreaOverlap;
    MPI_Allreduce(&dTotalAreaOverlap_loc, &dTotalAreaOverlap, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
    if (!pcomm->rank()) dbgprint.printf(0, "Overlap Mesh Area: %1.15e\n", dTotalAreaOverlap);

    // Partial cover
    if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
        if (!fNoCheck) {
            if (!pcomm->rank()) dbgprint.printf(0, "WARNING: Significant mismatch between overlap mesh area "
                "and input mesh area.\n  Automatically enabling --nocheck\n");
            fNoCheck = true;
        }
    }

/*
    // Recalculate input mesh area from overlap mesh
    if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
        dbgprint.printf(0, "Overlap mesh only covers a sub-area of the sphere\n");
        dbgprint.printf(0, "Recalculating source mesh areas\n");
        dTotalAreaInput = m_meshInput->CalculateFaceAreasFromOverlap(m_meshOverlap);
        dbgprint.printf(0, "New Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput);
    }
*/
    // Finite volume input / Finite volume output
    if ((eInputType  == DiscretizationType_FV) &&
        (eOutputType == DiscretizationType_FV)
    ) {

        // Generate reverse node array and edge map
        m_meshInputCov->ConstructReverseNodeArray();
        m_meshInputCov->ConstructEdgeMap();

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFV(*m_meshInputCov);
        this->InitializeTargetCoordinatesFromMeshFV(*m_meshOutput);

        // Construct OfflineMap
        // dbgprint.printf(0, "Number of source (%d) and target (%d) elements and overlap (%d)\n", m_meshInput->faces.size(), m_meshOutput->faces.size(), m_meshOverlap->faces.size());
        if (!pcomm->rank()) dbgprint.printf(0, "Calculating offline map\n");
        LinearRemapFVtoFV_Tempest_MOAB(nPin);
        // LinearRemapFVtoFV(*m_meshInput, *m_meshOutput, *m_meshOverlap, nPin, *this);

    // Finite volume input / Finite element output
    } else if (eInputType == DiscretizationType_FV) {
        DataMatrix3D<int> dataGLLNodes;
        DataMatrix3D<double> dataGLLJacobian;

        if (!pcomm->rank()) dbgprint.printf(0, "Generating output mesh meta data\n");
        double dNumericalArea_loc =
            GenerateMetaData(
                *m_meshOutput,
                nPout,
                fBubble,
                dataGLLNodes,
                dataGLLJacobian);

        Real dNumericalArea;
	    MPI_Allreduce(&dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
        if (!pcomm->rank()) dbgprint.printf(0, "Output Mesh Numerical Area: %1.15e\n", dNumericalArea);

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFV(*m_meshInputCov);
        this->InitializeTargetCoordinatesFromMeshFE(
            *m_meshOutput, nPout, dataGLLNodes);

        // Generate the continuous Jacobian
        bool fContinuous = (eOutputType == DiscretizationType_CGLL);

        if (eOutputType == DiscretizationType_CGLL) {
            GenerateUniqueJacobian(
                dataGLLNodes,
                dataGLLJacobian,
                this->GetTargetAreas());

        } else {
            GenerateDiscontinuousJacobian(
                dataGLLJacobian,
                this->GetTargetAreas());
        }

        // Generate reverse node array and edge map
        m_meshInputCov->ConstructReverseNodeArray();
        m_meshInputCov->ConstructEdgeMap();

        // Generate remap weights
        if (!pcomm->rank()) dbgprint.printf(0, "Calculating offline map\n");

        if (fVolumetric) {
            LinearRemapFVtoGLL_Volumetric(
                *m_meshInputCov,
                *m_meshOutput,
                *m_meshOverlap,
                dataGLLNodes,
                dataGLLJacobian,
                this->GetTargetAreas(),
                nPin,
                *this,
                nMonotoneType,
                fContinuous,
                fNoConservation);

        } else {
            LinearRemapFVtoGLL(
                *m_meshInputCov,
                *m_meshOutput,
                *m_meshOverlap,
                dataGLLNodes,
                dataGLLJacobian,
                this->GetTargetAreas(),
                nPin,
                *this,
                nMonotoneType,
                fContinuous,
                fNoConservation);
        }

    // Finite element input / Finite volume output
    } else if (
        (eInputType != DiscretizationType_FV) &&
        (eOutputType == DiscretizationType_FV)
    ) {
        DataMatrix3D<int> dataGLLNodes;
        DataMatrix3D<double> dataGLLJacobian;

        if (!pcomm->rank()) dbgprint.printf(0, "Generating input mesh meta data\n");
        double dNumericalArea_loc =
            GenerateMetaData(
                *m_meshInputCov,
                nPin,
                fBubble,
                dataGLLNodes,
                dataGLLJacobian);

        Real dNumericalArea;
	    MPI_Allreduce(&dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
	    if (!pcomm->rank()) dbgprint.printf(0, "Input Mesh Numerical Area: %1.15e\n", dNumericalArea);

        if (fabs(dNumericalArea - dTotalAreaInput) > 1.0e-12) {
            dbgprint.printf(0, "WARNING: Significant mismatch between input mesh "
                "numerical area and geometric area\n");
        }

        if (dataGLLNodes.GetSubColumns() != m_meshInputCov->faces.size()) {
            _EXCEPTIONT("Number of element does not match between metadata and "
                "input mesh");
        }

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFE(
            *m_meshInputCov, nPin, dataGLLNodes);
        this->InitializeTargetCoordinatesFromMeshFV(*m_meshOutput);

        // Generate the continuous Jacobian for input mesh
        bool fContinuousIn = (eInputType == DiscretizationType_CGLL);

        if (eInputType == DiscretizationType_CGLL) {
            GenerateUniqueJacobian(
                dataGLLNodes,
                dataGLLJacobian,
                this->GetSourceAreas());

        } else {
            GenerateDiscontinuousJacobian(
                dataGLLJacobian,
                this->GetSourceAreas());
        }

        // Generate offline map
        if (!pcomm->rank()) dbgprint.printf(0, "Calculating offline map\n");

        if (fVolumetric) {
            _EXCEPTIONT("Unimplemented: Volumetric currently unavailable for"
                "GLL input mesh");
        }

        LinearRemapSE4_Tempest_MOAB(
            // *m_meshInputCov,
            // *m_meshOutput,
            // *m_meshOverlap,
            dataGLLNodes,
            dataGLLJacobian,
            nMonotoneType,
            fContinuousIn,
            fNoConservation
            // *this
        );

    // Finite element input / Finite element output
    } else if (
        (eInputType  != DiscretizationType_FV) &&
        (eOutputType != DiscretizationType_FV)
    ) {
        DataMatrix3D<int> dataGLLNodesIn;
        DataMatrix3D<double> dataGLLJacobianIn;

        DataMatrix3D<int> dataGLLNodesOut;
        DataMatrix3D<double> dataGLLJacobianOut;

        // Input metadata
        if (!pcomm->rank()) dbgprint.printf(0, "Generating input mesh meta data");
        double dNumericalAreaIn_loc =
            GenerateMetaData(
                *m_meshInputCov,
                nPin,
                fBubble,
                dataGLLNodesIn,
                dataGLLJacobianIn);

        Real dNumericalAreaIn;
	    MPI_Allreduce(&dNumericalAreaIn_loc, &dNumericalAreaIn, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());
        if (!pcomm->rank()) dbgprint.printf(0, "Input Mesh Numerical Area: %1.15e", dNumericalAreaIn);

        if (fabs(dNumericalAreaIn - dTotalAreaInput) > 1.0e-12) {
            dbgprint.printf(0, "WARNING: Significant mismatch between input mesh "
                "numerical area and geometric area");
        }

        // Output metadata
        if (!pcomm->rank()) dbgprint.printf(0, "Generating output mesh meta data");
        double dNumericalAreaOut_loc =
            GenerateMetaData(
                *m_meshOutput,
                nPout,
                fBubble,
                dataGLLNodesOut,
                dataGLLJacobianOut);

        Real dNumericalAreaOut;
	    MPI_Allreduce(&dNumericalAreaOut_loc, &dNumericalAreaOut, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm());

        if (!pcomm->rank()) dbgprint.printf(0, "Output Mesh Numerical Area: %1.15e", dNumericalAreaOut);

        if (fabs(dNumericalAreaOut - dTotalAreaOutput) > 1.0e-12) {
            if (!pcomm->rank()) dbgprint.printf(0, "WARNING: Significant mismatch between output mesh "
                "numerical area and geometric area");
        }

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFE(
            *m_meshInputCov, nPin, dataGLLNodesIn);
        this->InitializeTargetCoordinatesFromMeshFE(
            *m_meshOutput, nPout, dataGLLNodesOut);

        // Generate the continuous Jacobian for input mesh
        bool fContinuousIn = (eInputType == DiscretizationType_CGLL);

        if (eInputType == DiscretizationType_CGLL) {
            GenerateUniqueJacobian(
                dataGLLNodesIn,
                dataGLLJacobianIn,
                this->GetSourceAreas());

        } else {
            GenerateDiscontinuousJacobian(
                dataGLLJacobianIn,
                this->GetSourceAreas());
        }

        // Generate the continuous Jacobian for output mesh
        bool fContinuousOut = (eOutputType == DiscretizationType_CGLL);

        if (eOutputType == DiscretizationType_CGLL) {
            GenerateUniqueJacobian(
                dataGLLNodesOut,
                dataGLLJacobianOut,
                this->GetTargetAreas());

        } else {
            GenerateDiscontinuousJacobian(
                dataGLLJacobianOut,
                this->GetTargetAreas());
        }

        // Generate offline map
        if (!pcomm->rank()) dbgprint.printf(0, "Calculating offline map");

        LinearRemapGLLtoGLL2(
            *m_meshInputCov,
            *m_meshOutput,
            *m_meshOverlap,
            dataGLLNodesIn,
            dataGLLJacobianIn,
            dataGLLNodesOut,
            dataGLLJacobianOut,
            this->GetTargetAreas(),
            nPin,
            nPout,
            nMonotoneType,
            fContinuousIn,
            fContinuousOut,
            fNoConservation,
            *this
        );

    } else {
        _EXCEPTIONT("Not implemented");
    }

//#pragma warning "NOTE: VERIFICATION DISABLED"

    // Verify consistency, conservation and monotonicity
    if (!fNoCheck) {
        if (!pcomm->rank()) dbgprint.printf(0, "Verifying map");
        this->IsConsistent(1.0e-8);
        if (!fNoConservation) this->IsConservative(1.0e-8);

        if (nMonotoneType != 0) {
            this->IsMonotone(1.0e-12);
        }
    }

    // Output the Offline Map
    if (strOutputMap != "") {
        if (!pcomm->rank()) dbgprint.printf(0, "Writing offline map");
        this->Write(strOutputMap);
    }

    // Apply Offline Map to data
    if (strInputData != "") {
        if (!pcomm->rank()) dbgprint.printf(0, "Applying offline map to data\n");

        this->SetFillValueOverride(static_cast<float>(dFillValueOverride));
        this->Apply(
            strInputData,
            strOutputData,
            vecVariableStrings,
            strNColName,
            fOutputDouble,
            false);
    }

    // Copy variables from input file to output file
    if ((strInputData != "") && (strOutputData != "")) {
        if (fPreserveAll) {
            if (!pcomm->rank()) dbgprint.printf(0, "Preserving variables");
            this->PreserveAllVariables(strInputData, strOutputData);

        } else if (vecPreserveVariableStrings.size() != 0) {
            if (!pcomm->rank()) dbgprint.printf(0, "Preserving variables");
            this->PreserveVariables(
                strInputData,
                strOutputData,
                vecPreserveVariableStrings);
        }
    }

} catch(Exception & e) {
    dbgprint.printf(0, "%s", e.ToString().c_str());
    return (moab::MB_FAILURE);

} catch(...) {
    return (moab::MB_FAILURE);
}
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsConsistent(
	double dTolerance
) {
	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

    // Calculate row sums
    DataVector<double> dRowSums;
	if (pcomm->size() > 1) {
        if (m_mapRemapGlobal.GetRows() == 0 || m_mapRemapGlobal.GetColumns() == 0)
            this->GatherAllToRoot();
        m_mapRemapGlobal.GetEntries(dataRows, dataCols, dataEntries);
        dRowSums.Initialize(m_mapRemapGlobal.GetRows());
        if (pcomm->size() > 1 && pcomm->rank() > 1) return true;
    }
    else {
       m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);
       dRowSums.Initialize(m_mapRemap.GetRows());
    }

	for (unsigned i = 0; i < dataRows.GetRows(); i++) {
		dRowSums[dataRows[i]] += dataEntries[i];
	}

	// Verify all row sums are equal to 1
	bool fConsistent = true;
	for (unsigned i = 0; i < dRowSums.GetRows(); i++) {
		if (fabs(dRowSums[i] - 1.0) > dTolerance) {
			fConsistent = false;
			Announce("TempestOfflineMap is not consistent in row %i (%1.15e)",
				i, dRowSums[i]);
		}
	}

	return fConsistent;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsConservative(
	double dTolerance
) {
	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;
	const DataVector<double>& dTargetAreas = this->GetGlobalTargetAreas();
	const DataVector<double>& dSourceAreas = this->GetGlobalSourceAreas();

    // Calculate column sums
    DataVector<double> dColumnSums;
    if (pcomm->size() > 1) {
        if (m_mapRemapGlobal.GetRows() == 0 || m_mapRemapGlobal.GetColumns() == 0)
            this->GatherAllToRoot();
        if (pcomm->size() > 1 && pcomm->rank()) return true;
        m_mapRemapGlobal.GetEntries(dataRows, dataCols, dataEntries);
        dColumnSums.Initialize(m_mapRemapGlobal.GetColumns());
    }
    else {
       m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);
       dColumnSums.Initialize(m_mapRemap.GetColumns());
    }

	for (unsigned i = 0; i < dataRows.GetRows(); i++) {
		dColumnSums[dataCols[i]] +=
			dataEntries[i] * dTargetAreas[dataRows[i]];
	}

	// Verify all column sums equal the input Jacobian
	bool fConservative = true;
	for (unsigned i = 0; i < dColumnSums.GetRows(); i++) {
		if (fabs(dColumnSums[i] - dSourceAreas[i]) > dTolerance) {
			fConservative = false;
			Announce("TempestOfflineMap is not conservative in column "
				"%i (%1.15e / %1.15e)",
				i, dColumnSums[i], dSourceAreas[i]);
		}
	}

	return fConservative;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsMonotone(
	double dTolerance
) {

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

    if (pcomm->size() > 1) {
        if (m_mapRemapGlobal.GetRows() == 0 || m_mapRemapGlobal.GetColumns() == 0)
            this->GatherAllToRoot();
        m_mapRemapGlobal.GetEntries(dataRows, dataCols, dataEntries);
        if (pcomm->size() > 1 && pcomm->rank() > 1) return true;
    }
    else
	   m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Verify all entries are in the range [0,1]
	bool fMonotone = true;
	for (unsigned i = 0; i < dataRows.GetRows(); i++) {
		if ((dataEntries[i] < -dTolerance) ||
			(dataEntries[i] > 1.0 + dTolerance)
		) {
			fMonotone = false;

			Announce("TempestOfflineMap is not monotone in entry (%i): %1.15e",
				i, dataEntries[i]);
		}
	}

	return fMonotone;
}

///////////////////////////////////////////////////////////////////////////////
moab::ErrorCode moab::TempestOfflineMap::GatherAllToRoot() { // Collective

	Mesh globalMesh;
	int ierr, rootProc = 0;
    moab::ErrorCode rval;

	// Write SparseMatrix entries
	DataVector<int> vecRow;
	DataVector<int> vecCol;
	DataVector<double> vecS;

	moab::DebugOutput dbgprint(std::cout, (pcomm ? pcomm->rank() : 0));

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);
    const DataVector<double>& dSourceAreas = m_meshInputCov->vecFaceArea;
    const DataVector<double>& dTargetAreas = m_meshOutput->vecFaceArea;

    // Translate the index in Row and Col to global_id and dump it out

	// Communicate the necessary data
    int grows=0, gcols=0, gsrc=0, gtar=0;
	std::vector<int> rowcolss, rowcolsv;
    DataVector<int> rows, cols, srcelmindx, tgtelmindx;
	{	// First, accumulate the sizes of rows and columns of the matrix
		if (!pcomm->rank()) rowcolss.resize(pcomm->size()*4);

		int sendarray[4];
		sendarray[0] = vecRow.GetRows();
		sendarray[1] = vecCol.GetRows();
        sendarray[2] = dSourceAreas.GetRows();
        sendarray[3] = dTargetAreas.GetRows();
		
		ierr = MPI_Gather( sendarray, 4, MPI_INTEGER, rowcolss.data(), 4, MPI_INTEGER, rootProc, pcomm->comm());

		// dbgprint.printf(0, "[%d] Dimensions: %d, %d\n", pcomm->rank(), vecRow.GetRows(), vecCol.GetRows());

		if (!pcomm->rank()) {
			for (unsigned i=0; i < pcomm->size(); ++i) {
                grows += rowcolss[4*i];
                gcols += rowcolss[4*i+1];
                gsrc += rowcolss[4*i+2];
                gtar += rowcolss[4*i+3];
            }
			rowcolsv.resize(grows+gcols+gsrc+gtar);
            rows.Initialize(grows); // we are assuming rows = cols
            cols.Initialize(gcols); // we are assuming rows = cols
            
            m_areasSrcGlobal.Initialize(gsrc); srcelmindx.Initialize(gsrc);
            m_areasTgtGlobal.Initialize(gtar); tgtelmindx.Initialize(gtar);
            dbgprint.printf(0, "Received global dimensions: %d, %d\n", vecRow.GetRows(), rows.GetRows());
            dbgprint.printf(0, "Global: n(source) = %d, and n(target) = %d\n", gsrc, gtar);
			dbgprint.printf(0, "Operator size = %d\n", grows+gcols);
		}
	}

    {
        std::stringstream sstr;
        sstr << "rowscols_" << pcomm->rank() << ".txt";
        std::ofstream output_file(sstr.str().c_str());
        output_file << "VALUES\n";
        for (unsigned ip=0; ip < vecRow.GetRows(); ++ip) {
            output_file << ip << " (" << vecRow[ip] << ", " << vecCol[ip] << ") = " << vecS[ip] << "\n";
        }
        output_file.flush(); // required here
        output_file.close();
    }

	{
        // Next, accumulate the row and column values for the matrices (in local indexing)
		const int nR = vecRow.GetRows()+vecCol.GetRows()+dSourceAreas.GetRows()+dTargetAreas.GetRows();
		std::vector<int> sendarray(nR);
        for (unsigned ix=0; ix < vecRow.GetRows(); ++ix) {
            sendarray[ix] = m_remapper->GetGlobalID(moab::Remapper::TargetMesh, vecRow[ix]);
        }
        for (unsigned ix=0, offset=vecRow.GetRows(); ix < vecCol.GetRows(); ++ix) {
            sendarray[offset+ix] = m_remapper->GetGlobalID(moab::Remapper::CoveringMesh, vecCol[ix]);
        }
        // for (unsigned ix=0, offset=vecRow.GetRows()+vecCol.GetRows(); ix < dSourceAreas.GetRows(); ++ix) {
        //     sendarray[offset+ix] = m_remapper->lid_to_gid_src[ix];
        //     if (!pcomm->rank()) std::cout << "[ " << offset+ix << "] Source: local id " << ix << " and global ID = " << sendarray[ix+offset] << "\n";
        // }
        // for (unsigned ix=0, offset=vecRow.GetRows()+vecCol.GetRows()+dSourceAreas.GetRows(); ix < dTargetAreas.GetRows(); ++ix) {
        //     sendarray[offset+ix] = m_remapper->lid_to_gid_tgt[ix];
        // }
        {
            moab::Tag gidtag;
            rval = mbCore->tag_get_handle("GLOBAL_ID", gidtag);MB_CHK_ERR(rval);
            // rval = mbCore->tag_get_data(gidtag, m_remapper->GetEntities(moab::Remapper::SourceMesh), &sendarray[vecRow.GetRows()+vecCol.GetRows()]);MB_CHK_ERR(rval);
            rval = mbCore->tag_get_data(gidtag, m_remapper->m_covering_source_entities, &sendarray[vecRow.GetRows()+vecCol.GetRows()]);MB_CHK_ERR(rval);
            // rval = mbCore->tag_get_data(gidtag, m_remapper->m_source_entities, &sendarray[vecRow.GetRows()+vecCol.GetRows()]);MB_CHK_ERR(rval);
            rval = mbCore->tag_get_data(gidtag, m_remapper->m_target_entities, &sendarray[vecRow.GetRows()+vecCol.GetRows()+dSourceAreas.GetRows()]);MB_CHK_ERR(rval);
        }

		std::vector<int> displs, rcount;
		if (!pcomm->rank()) {
			displs.resize(pcomm->size(),0);
			rcount.resize(pcomm->size(),0);
			int gsum=0;
			for (unsigned i=0; i < pcomm->size(); ++i) {
				displs[i] = gsum;
				rcount[i] = rowcolss[4*i]+rowcolss[4*i+1]+rowcolss[4*i+2]+rowcolss[4*i+3];
				gsum += rcount[i];
			}
            assert(rowcolsv.size() - gsum == 0);
		}

		// Both rows and columns have a size of "rowsize"
		ierr = MPI_Gatherv( &sendarray[0], nR, MPI_INTEGER, &rowcolsv[0], &rcount[0], &displs[0], MPI_INTEGER, rootProc, pcomm->comm());

        if (!pcomm->rank()) {
		    std::ofstream output_file("rows-cols.txt", std::ios::out);
		    output_file << "ROWS\n";
		    for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip], iend = istart + rowcolss[4*ip];
			    for (int i = istart; i < iend; ++i, ++offset) {
			    	output_file << offset << " " << rowcolsv[i] << "\n";
                    rows[offset] = rowcolsv[i];
			    }
			}
			output_file << "COLS\n";
			for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip] + rowcolss[4*ip], iend = istart + rowcolss[4*ip+1];
                for (int i = istart; i < iend; ++i, ++offset) {
			    	output_file << offset << " " << rowcolsv[i] << "\n";
                    cols[offset] = rowcolsv[i];
			    }
			}
			output_file.flush(); // required here
			output_file.close();
            for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip] + rowcolss[4*ip] + rowcolss[4*ip+1], iend = istart + rowcolss[4*ip+2];
                for (int i = istart; i < iend; ++i, ++offset) {
                    srcelmindx[offset] = rowcolsv[i];
                }
            }
            for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip] + rowcolss[4*ip] + rowcolss[4*ip+1] + rowcolss[4*ip+2], iend = istart + rowcolss[4*ip+3];
                for (int i = istart; i < iend; ++i, ++offset) {
                    tgtelmindx[offset] = rowcolsv[i];
                }
            }
		} /* (!pcomm->rank()) */
        rowcolsv.clear();
	}

    {
        // Next, accumulate the row and column values for the matrices (in local indexing)
        const int nR = vecS.GetRows()+dSourceAreas.GetRows()+dTargetAreas.GetRows();

        std::vector<double> sendarray(nR);
        std::copy((double*)vecS, (double*)vecS+vecS.GetRows(), sendarray.begin());
        std::copy((const double*)dSourceAreas, (const double*)dSourceAreas+dSourceAreas.GetRows(), sendarray.begin()+vecS.GetRows());
        std::copy((const double*)dTargetAreas, (const double*)dTargetAreas+dTargetAreas.GetRows(), sendarray.begin()+vecS.GetRows()+dSourceAreas.GetRows());

        std::vector<int> displs, rcount;
        int gsum=0;
        if (!pcomm->rank()) {
            displs.resize(pcomm->size(),0);
            rcount.resize(pcomm->size(),0);
            for (unsigned i=0; i < pcomm->size(); ++i) {
                displs[i] = gsum;
                rcount[i] = rowcolss[4*i]+rowcolss[4*i+2]+rowcolss[4*i+3];
                gsum += rcount[i];
            }
        }

        std::vector<double> rowcolsvals(grows+gsrc+gtar);

        // Both rows and columns have a size of "rowsize"
        ierr = MPI_Gatherv( sendarray.data(), sendarray.size(), MPI_DOUBLE, rowcolsvals.data(), &rcount[0], &displs[0], MPI_DOUBLE, rootProc, pcomm->comm());

        if (!pcomm->rank()) {
            std::ofstream output_file("rows-cols.txt", std::ios::app);
            DataVector<double> globvalues(grows);
            output_file << "VALUES\n";
            for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                for (int i=displs[ip]; i < displs[ip]+rowcolss[4*ip]; ++i, ++offset) {
                    globvalues[offset] = rowcolsvals[i];
                    output_file << offset << " (" << rows[offset] << ", " << cols[offset] << ") = " << rowcolsvals[i] << "\n";
                }
            }
            output_file.flush(); // required here
            output_file.close();
            m_mapRemapGlobal.SetEntries(rows, cols, globvalues);
        }

        if (!pcomm->rank()) {
            // Store the global source and target elements areas
            std::ofstream output_file("source-target-areas.txt");
            output_file << "Source areas\n";
            for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip] + rowcolss[4*ip], iend = istart + rowcolss[4*ip+2];
                for (int i=istart; i < iend; ++i, ++offset) {
                    output_file << srcelmindx[offset] << " " << rowcolsvals[i] << "\n";
                    m_areasSrcGlobal[srcelmindx[offset]] = rowcolsvals[i];
                    // m_areasSrcGlobal[offset] = rowcolsvals[i];
                }
            }
            output_file << "Target areas\n";
            for (unsigned ip=0, offset=0; ip < pcomm->size(); ++ip) {
                int istart = displs[ip] + rowcolss[4*ip] + rowcolss[4*ip+2], iend = istart + rowcolss[4*ip+3];
                for (int i=istart; i < iend; ++i, ++offset) {
                    output_file << tgtelmindx[offset] << " " << rowcolsvals[i] << "\n";
                    m_areasTgtGlobal[tgtelmindx[offset]] = rowcolsvals[i];
                    // m_areasTgtGlobal[offset] = rowcolsvals[i];
                }
            }
            output_file.flush(); // required here
            output_file.close();
        }

    }
    rowcolss.clear();

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
