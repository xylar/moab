///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOfflineMap.cpp
///	\author  Paul Ullrich
///	\version June 29, 2015
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"
#include "DataMatrix3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"

#include "OverlapMesh.h"
#include "LinearRemapSE0.h"
#include "LinearRemapFV.h"
#include "TempestOfflineMap.hpp"

#include "DebugOutput.hpp"

#include "netcdfcpp.h"

#include <fstream>
#include <cmath>


///////////////////////////////////////////////////////////////////////////////

TempestOfflineMap::TempestOfflineMap(moab::TempestRemapper* remapper) : OfflineMap(), m_remapper(remapper)
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
    this->InitializeSourceDimensionsFromMesh(*m_meshInput);
    dbgprint.printf(0, "Output mesh\n");
    this->InitializeTargetDimensionsFromMesh(*m_meshOutput);
    // dbgprint.printf(0, "----------------------------------\n");

    // Build a matrix of source and target discretization so that we know how to assign
    // the global DoFs in parallel for the mapping weights
    // For example, FV->FV: rows X cols = faces_source X faces_target
    // 
    // int ierr;
    // MPI_Scan()
}

///////////////////////////////////////////////////////////////////////////////

TempestOfflineMap::~TempestOfflineMap()
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
#include "LinearRemapFV.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MeshUtilitiesFuzzy.h"
#include "OverlapMesh.h"

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

void TempestOfflineMap::LinearRemapFVtoFV_Tempest_MOAB(
    int nOrder
) {
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if (m_meshInput->revnodearray.size() == 0) {
        _EXCEPTIONT("ReverseNodeArray has not been calculated for m_meshInput");
    }

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

    // Number of elements needed
    const int nCoefficients = nOrder * (nOrder + 1) / 2;

#pragma message "This should be a command-line parameter"
    // Number of faces you need
    const int nRequiredFaceSetSize = nCoefficients;

    // Fit weight exponent
    const int nFitWeightsExponent = nOrder + 2;

    // Announcemnets
    Announce("Triangular quadrature rule order %i", TriQuadRuleOrder);
    Announce("Number of coefficients: %i", nCoefficients);
    Announce("Required adjacency set size: %i", nRequiredFaceSetSize);
    Announce("Fit weights exponent: %i", nFitWeightsExponent);

    // Current overlap face
    int ixOverlap = 0;

    // for (unsigned ix = 0; ix < m_meshOverlap->faces.size(); ix++) {
    //     if (!pcomm->rank()) std::cout << "[" << ix << "] SrcFaceIx = " << m_meshOverlap->vecSourceFaceIx[ix] 
    //                                   << " TargetFaceIx = " << m_meshOverlap->vecTargetFaceIx[ix] << std::endl;
    // }
    // return ;

    // Loop through all faces on m_meshInput
    for (int ixFirst = 0; ixFirst < m_meshInput->faces.size(); ixFirst++) {

        // Output every 100 elements
        if (ixFirst % 100 == 0) {
            Announce("Element %i/%i", ixFirst, m_meshInput->faces.size());
        }

        // This Face
        // const Face & faceFirst = m_meshInput->faces[ixFirst];

        std::vector<std::pair<int,int> > ixFaces;
        for (unsigned ix = 0; ix < m_meshOverlap->faces.size(); ix++) {
            if (m_meshOverlap->vecSourceFaceIx[ix] == ixFirst) {
                // if (!pcomm->rank()) std::cout << "[" << ix << "] SrcFaceIx = " << m_meshOverlap->vecSourceFaceIx[ix] << std::endl;
                unsigned jx;
                for (jx = 1; jx < m_meshOverlap->faces.size()-ix; jx++) {
                    // if (!pcomm->rank()) std::cout << "\t[" << jx << "] SrcFaceJx = " << m_meshOverlap->vecSourceFaceIx[jx] << std::endl;
                    if (m_meshOverlap->vecSourceFaceIx[ix+jx] != ixFirst)
                        break;
                }
                ixFaces.push_back(std::make_pair<int,int>(ix,ix+jx));
                ix += jx-1;
            }
        }

        // int nOverlapFaces = ixFaces.size();

        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so on wrt the input mesh data.
        // Then the overlap_end and overlap_begin will be correct. However, the relation with MOAB and Tempest will go out of the roof.

        // int ifaceIndx = 0;
        int gnOverlapFaces = 0;
        for (unsigned ifaceIndx = 0; ifaceIndx < ixFaces.size(); ++ifaceIndx) {

            int ixOverlapBegin = ixFaces[ifaceIndx].first;
            int ixOverlapEnd = ixFaces[ifaceIndx].second;
            int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

            // Build integration array
            DataMatrix<double> dIntArray;

            BuildIntegrationArray(
                *m_meshInput,
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
                *m_meshInput,
                ixFirst,
                nRequiredFaceSetSize,
                vecAdjFaces);

            // Number of adjacent Faces
            int nAdjFaces = vecAdjFaces.size();

            // Determine the conservative constraint equation
            DataVector<double> dConstraint;

            dConstraint.Initialize(nCoefficients);

            double dFirstArea = m_meshInput->vecFaceArea[ixFirst];

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
                *m_meshInput,
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
        
            // Announce("[%i] (%d,%d) ixOverlapBegin: %i and ixOverlapEnd: %i -- Area: %12.10f", ixFirst, nAdjFaces, nOverlapFaces, m_meshOverlap->vecSourceFaceIx[ixOverlapBegin], m_meshOverlap->vecTargetFaceIx[ixOverlapEnd], dFirstArea);

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
                int ixFirstFace = vecAdjFaces[i].first;
                int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlapBegin + j];

                // std::cout << "!"<<gnOverlapFaces<<"! (" << ixSecondFace << ", " << ixFirstFace << ") = " << dComposedArray[i][j]/m_meshOutput->vecFaceArea[ixSecondFace] << "\t";

                m_mapRemap(ixSecondFace, ixFirstFace) +=
                    dComposedArray[i][j]
                    / m_meshOutput->vecFaceArea[ixSecondFace];
            }
            }

            gnOverlapFaces += nOverlapFaces;

        }

        // Increment the current overlap index
        ixOverlap += gnOverlapFaces;
        // break;
    }
}


///////////////////////////////////////////////////////////////////////////////
moab::ErrorCode TempestOfflineMap::GenerateOfflineMap( std::string strInputType, std::string strOutputType,
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

    // Input mesh areas
    if (eInputType == DiscretizationType_FV) {
        this->SetSourceAreas(m_meshInput->vecFaceArea);
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

    dbgprint.printf(0, "m_meshInputCov->faces = %d, m_meshOutput->faces = %d, ixSourceFaceMax = %d\n", m_meshInput->faces.size(), m_meshOutput->faces.size(), ixSourceFaceMax);

    // Check for forward correspondence in overlap mesh
    if (m_meshInput->faces.size() - ixSourceFaceMax == 0 //&&
        //(ixTargetFaceMax == m_meshOutput.faces.size())
    ) {
        if (!pcomm->rank()) dbgprint.printf(0, "Overlap mesh forward correspondence found\n");

    // Check for reverse correspondence in overlap mesh
    } else if (
        m_meshOutput->faces.size() - ixSourceFaceMax == 0 //&&
        //(ixTargetFaceMax == m_meshInput->faces.size())
    ) {
        if (!pcomm->rank()) dbgprint.printf(0, "Overlap mesh reverse correspondence found (reversing)\n");

        // Reorder overlap mesh
        m_meshOverlap->ExchangeFirstAndSecondMesh();

    // No correspondence found
    } else {
        _EXCEPTION2("Invalid overlap mesh:\n"
            "    No correspondence found with input and output meshes (%i,%i)",
            ixSourceFaceMax, ixTargetFaceMax);
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
        m_meshInput->ConstructReverseNodeArray();
        m_meshInput->ConstructEdgeMap();

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFV(*m_meshInput);
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
        this->InitializeSourceCoordinatesFromMeshFV(*m_meshInput);
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
        m_meshInput->ConstructReverseNodeArray();
        m_meshInput->ConstructEdgeMap();

        // Generate remap weights
        if (!pcomm->rank()) dbgprint.printf(0, "Calculating offline map\n");

        if (fVolumetric) {
            LinearRemapFVtoGLL_Volumetric(
                *m_meshInput,
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
                *m_meshInput,
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
                *m_meshInput,
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

        if (dataGLLNodes.GetSubColumns() != m_meshInput->faces.size()) {
            _EXCEPTIONT("Number of element does not match between metadata and "
                "input mesh");
        }

        // Initialize coordinates for map
        this->InitializeSourceCoordinatesFromMeshFE(
            *m_meshInput, nPin, dataGLLNodes);
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

        LinearRemapSE4(
            *m_meshInput,
            *m_meshOutput,
            *m_meshOverlap,
            dataGLLNodes,
            dataGLLJacobian,
            nMonotoneType,
            fContinuousIn,
            fNoConservation,
            *this
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
                *m_meshInput,
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
            *m_meshInput, nPin, dataGLLNodesIn);
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
            *m_meshInput,
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
        this->IsConservative(1.0e-8);

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

bool TempestOfflineMap::IsConsistent(
	double dTolerance
) {

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Calculate row sums
	DataVector<double> dRowSums;
	dRowSums.Initialize(m_mapRemap.GetRows());

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

bool TempestOfflineMap::IsConservative(
	double dTolerance
) {
/*
	if (vecSourceAreas.GetRows() != m_mapRemap.GetColumns()) {
		_EXCEPTIONT("vecSourceAreas / mapRemap dimension mismatch");
	}
	if (vecTargetAreas.GetRows() != m_mapRemap.GetRows()) {
		_EXCEPTIONT("vecTargetAreas / mapRemap dimension mismatch");
	}
*/
	moab::ErrorCode rval;
	moab::Tag consTag;
	Real dum=0.0;
	rval = mbCore->tag_get_handle("CONSERVED", 1, moab::MB_TYPE_DOUBLE, 
		                           consTag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &dum);MB_CHK_SET_ERR(rval, "can't get conservation tag");

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;
	const DataVector<double>& dTargetAreas = this->GetTargetAreas();
	const DataVector<double>& dSourceAreas = this->GetSourceAreas();

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Calculate column sums
	DataVector<double> dColumnSums;
	dColumnSums.Initialize(m_mapRemap.GetColumns());

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

	rval = mbCore->tag_delete(consTag);MB_CHK_SET_ERR(rval, "can't delete conservation tag");

	return fConservative;
}

///////////////////////////////////////////////////////////////////////////////

bool TempestOfflineMap::IsMonotone(
	double dTolerance
) {

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

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

void TempestOfflineMap::Write(
	const std::string & strTarget
) {
	NcFile ncMap(strTarget.c_str(), NcFile::Replace);
	if (!ncMap.is_valid()) {
		_EXCEPTION1("Unable to open output map file \"%s\"",
			strTarget.c_str());
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
	// varSrcGridDims->DoCollectiveIO (true);
	NcVar * varDstGridDims =
		ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);
	// varDstGridDims->DoCollectiveIO (true);

	char szDim[64];
	if ((nSrcGridDims == 1) && (m_vecSourceDimSizes[0] != nA)) {
		varSrcGridDims->put(&nA, 1);
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

	if ((nDstGridDims == 1) && (m_vecSourceDimSizes[0] != nB)) {
		varDstGridDims->put(&nB, 1);
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
	if (m_dSourceCenterLon.GetRows() - nA != 0) {
		_EXCEPTIONT("Mismatch between m_dSourceCenterLon and nA");
	}
	if (m_dSourceCenterLat.GetRows() - nA != 0) {
		_EXCEPTIONT("Mismatch between m_dSourceCenterLat and nA");
	}
	if (m_dTargetCenterLon.GetRows() - nB != 0) {
		_EXCEPTIONT("Mismatch between m_dTargetCenterLon and nB");
	}
	if (m_dTargetCenterLat.GetRows() - nB != 0) {
		_EXCEPTIONT("Mismatch between m_dTargetCenterLat and nB");
	}
	if (m_dSourceVertexLon.GetRows() - nA != 0) {
		_EXCEPTIONT("Mismatch between m_dSourceVertexLon and nA");
	}
	if (m_dSourceVertexLat.GetRows() - nA != 0) {
		_EXCEPTIONT("Mismatch between m_dSourceVertexLat and nA");
	}
	if (m_dTargetVertexLon.GetRows() - nB != 0) {
		_EXCEPTIONT("Mismatch between m_dTargetVertexLon and nB");
	}
	if (m_dTargetVertexLat.GetRows() - nB != 0) {
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

	NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
	varAreaB->put(&(m_dTargetAreas[0]), nB);

	// Write frac
	DataVector<double> dFrac;

	dFrac.Initialize(nA);
	for (int i = 0; i < nA; i++) {
		dFrac[i] = 1.0;
	}
	NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
	varFracA->put(&(dFrac[0]), nA);

	dFrac.Initialize(nB);
	for (int i = 0; i < nB; i++) {
		dFrac[i] = 1.0;
	}
	NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
	varFracB->put(&(dFrac[0]), nB);

	// Write SparseMatrix entries
	DataVector<int> vecRow;
	DataVector<int> vecCol;
	DataVector<double> vecS;

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);

	// Increment vecRow and vecCol
	for (unsigned i = 0; i < vecRow.GetRows(); i++) {
		vecRow[i]++;
		vecCol[i]++;
	}

	// Load in data
	int nS = vecRow.GetRows();
	NcDim * dimNS = ncMap.add_dim("n_s", nS);

	NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
	NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
	NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);

	varRow->set_cur((long)0);
	varRow->put(&(vecRow[0]), nS);

	varCol->set_cur((long)0);
	varCol->put(&(vecCol[0]), nS);

	varS->set_cur((long)0);
	varS->put(&(vecS[0]), nS);
}

///////////////////////////////////////////////////////////////////////////////

void TempestOfflineMap::GatherAllToRoot() {

	Mesh globalMesh;
	int ierr, rootProc = 0;

	// Write SparseMatrix entries
	DataVector<int> vecRow;
	DataVector<int> vecCol;
	DataVector<double> vecS;

	moab::DebugOutput dbgprint(std::cout, (pcomm ? pcomm->rank() : 0));

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);

	// Increment vecRow and vecCol
	// for (unsigned i = 0; i < vecRow.GetRows(); i++) {
	// 	vecRow[i]++;
	// 	vecCol[i]++;
	// }

    // Translate the index in Row and Col to global_id and dump it out

	// Communicate the necessary data
	std::vector<int> rowcolss, rowcolsv;
	{	// First, accumulate the sizes of rows and columns of the matrix
		if (!pcomm->rank()) rowcolss.resize(pcomm->size()*2);

		int sendarray[2];
		sendarray[0] = vecRow.GetRows();
		sendarray[1] = vecCol.GetRows();
		
		ierr = MPI_Gather( sendarray, 2, MPI_INTEGER, &rowcolss[0], 2, MPI_INTEGER, rootProc, pcomm->comm());

		dbgprint.printf(0, "[%D] Dimensions: %D, %D\n", pcomm->rank(), vecRow.GetRows(), vecCol.GetRows());

		if (!pcomm->rank()) {
			int gsize=0;
			for (unsigned i=0; i < rowcolss.size(); ++i) gsize += rowcolss[i];
			rowcolsv.resize(gsize);
			dbgprint.printf(0, "Resizing rowscolv to %D\n", gsize);
		}
	}

	{ // Next, accumulate the row and column values for the matrices (in local indexing)
		const int nR = vecRow.GetRows()+vecCol.GetRows();
		std::vector<int> sendarray(nR);
		std::copy((int*)vecRow, (int*)vecRow+vecRow.GetRows(), sendarray.begin());
		std::copy((int*)vecCol, (int*)vecCol+vecCol.GetRows(), sendarray.begin()+vecRow.GetRows());

		std::vector<int> displs, rcount;
		if (!pcomm->rank()) {
			displs.resize(pcomm->size(),0);
			rcount.resize(pcomm->size(),0);
			int gsum=0;
			for (unsigned i=0; i < pcomm->size(); ++i) {
				displs[i] = gsum;
				rcount[i] = rowcolss[2*i]+rowcolss[2*i+1];
				gsum += rcount[i];
			}
			dbgprint.printf(0, "Received global dimensions: %D, %D\n", nR, gsum);
		}

		// Both rows and columns have a size of "rowsize"
		ierr = MPI_Gatherv( &sendarray[0], nR, MPI_INTEGER, &rowcolsv[0], &rcount[0], &displs[0], MPI_INTEGER, rootProc, pcomm->comm());

		if (!pcomm->rank()) {
		    std::ofstream output_file("rows-cols.txt", std::ios::out);
		    output_file << "ROWS\n";
            int offset=0;
		    for (unsigned ip=0; ip < pcomm->size(); ++ip) {
			    for (int i=displs[ip]; i < displs[ip]+rowcolss[2*ip]; ++i, ++offset) {
			    	output_file << offset << " " << rowcolsv[i] << "\n";
			    }
			}
			output_file << "COLS\n";
            offset=0;
			for (unsigned ip=0; ip < pcomm->size(); ++ip) {
			    for (int i=displs[ip]+rowcolss[2*ip]; i < displs[ip]+rowcolss[2*ip]+rowcolss[2*ip+1]; ++i, ++offset) {
			    	output_file << offset << " " << rowcolsv[i] << "\n";
			    }
			}
		    // std::ostream_iterator<int> output_iterator(output_file, "\n");
		    // std::copy(rowcolsv.begin(), rowcolsv.end(), output_iterator);
			output_file.flush(); // required here
			output_file.close();
		}
	}

    { // Next, accumulate the row and column values for the matrices (in local indexing)
        const int nR = vecS.GetRows();
        // std::vector<double> sendarray(nR);
        // std::copy((double*)vecS, (double*)vecS+vecS.GetRows(), sendarray.begin());
        std::vector<double> rowcolsvals;

        std::vector<int> displs, rcount;
        int gsum=0;
        if (!pcomm->rank()) {
            displs.resize(pcomm->size(),0);
            rcount.resize(pcomm->size(),0);
            for (unsigned i=0; i < pcomm->size(); ++i) {
                displs[i] = gsum;
                rcount[i] = rowcolss[2*i];
                gsum += rcount[i];
            }
            rowcolsvals.resize(gsum);
            // dbgprint.printf(0, "Received global dimensions: %D, %D\n", nR, gsum);
        }

        // Both rows and columns have a size of "rowsize"
        ierr = MPI_Gatherv( (double*)vecS, nR, MPI_DOUBLE, &rowcolsvals[0], &rcount[0], &displs[0], MPI_DOUBLE, rootProc, pcomm->comm());

        if (!pcomm->rank()) {
            std::ofstream output_file("rows-cols.txt", std::ios::app);
            output_file << "VALUES\n";
            for (int ip=0; ip < gsum; ++ip)
                output_file << ip << " " << rowcolsvals[ip] << "\n";
            // std::ostream_iterator<int> output_iterator(output_file, "\n");
            // std::copy(rowcolsv.begin(), rowcolsv.end(), output_iterator);
            output_file.flush(); // required here
            output_file.close();
        }
    }

}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
