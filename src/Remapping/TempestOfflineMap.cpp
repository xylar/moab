/*
 * =====================================================================================
 *
 *       Filename:  TempestOfflineMap.hpp
 *
 *    Description:  Interface to the TempestRemap library to compute the consistent,
 *                  and accurate high-order conservative remapping weights for overlap
 *                  grids on the sphere in climate simulations.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#include "Announce.h"
#include "DataMatrix3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"

#include "moab/Remapping/TempestOfflineMap.hpp"
#include "DebugOutput.hpp"

#include <fstream>
#include <cmath>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

// #define VERBOSE

///////////////////////////////////////////////////////////////////////////////

moab::TempestOfflineMap::TempestOfflineMap ( moab::TempestRemapper* remapper ) : OfflineMap(), m_remapper ( remapper )
{
    // Get the references for the MOAB core objects
    mbCore = m_remapper->get_interface();
    pcomm = m_remapper->get_parallel_communicator();

    // Update the references to the meshes
    m_meshInput = remapper->GetMesh ( moab::Remapper::SourceMesh );
    m_meshInputCov = remapper->GetCoveringMesh();
    m_meshOutput = remapper->GetMesh ( moab::Remapper::TargetMesh );
    m_meshOverlap = remapper->GetMesh ( moab::Remapper::IntersectedMesh );
    m_globalMapAvailable = false;

    moab::DebugOutput dbgprint ( std::cout, ( pcomm ? pcomm->rank() : 0 ) );

    // Compute and store the total number of source and target DoFs corresponding
    // to number of rows and columns in the mapping.

    // Initialize dimension information from file
    dbgprint.printf ( 0, "Initializing dimensions of map\n" );
    dbgprint.printf ( 0, "Input mesh\n" );
    this->InitializeSourceDimensionsFromMesh ( *m_meshInputCov );
    dbgprint.printf ( 0, "Output mesh\n" );
    this->InitializeTargetDimensionsFromMesh ( *m_meshOutput );

    // Build a matrix of source and target discretization so that we know how to assign
    // the global DoFs in parallel for the mapping weights
    // For example, FV->FV: rows X cols = faces_source X faces_target
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

static void ParseVariableList (
    const std::string & strVariables,
    std::vector< std::string > & vecVariableStrings
)
{
    unsigned iVarBegin = 0;
    unsigned iVarCurrent = 0;

    // Parse variable name
    for ( ;; )
    {
        if ( ( iVarCurrent >= strVariables.length() ) ||
                ( strVariables[iVarCurrent] == ',' ) ||
                ( strVariables[iVarCurrent] == ' ' )
           )
        {
            if ( iVarCurrent == iVarBegin )
            {
                if ( iVarCurrent >= strVariables.length() )
                {
                    break;
                }
                continue;
            }

            vecVariableStrings.push_back (
                strVariables.substr ( iVarBegin, iVarCurrent - iVarBegin ) );

            iVarBegin = iVarCurrent + 1;
        }

        iVarCurrent++;
    }
}


///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOfflineMap::GenerateOfflineMap ( std::string strInputType, std::string strOutputType,
        int nPin, int nPout,
        bool fBubble, int fMonotoneTypeID,
        bool fVolumetric, bool fNoConservation, bool fNoCheck,
        std::string strVariables,
        std::string strInputData, std::string strOutputData,
        std::string strNColName, bool fOutputDouble,
        std::string strPreserveVariables, bool fPreserveAll, double dFillValueOverride,
        bool fInputConcave, bool fOutputConcave )
{
    NcError error ( NcError::silent_nonfatal );

    moab::DebugOutput dbgprint ( std::cout, ( pcomm ? pcomm->rank() : 0 ) );

    try
    {
        // Input / Output types
        enum DiscretizationType
        {
            DiscretizationType_FV,
            DiscretizationType_CGLL,
            DiscretizationType_DGLL
        };

        // Check command line parameters (data arguments)
        if ( ( strInputData != "" ) && ( strOutputData == "" ) )
        {
            _EXCEPTIONT ( "--in_data specified without --out_data" );
        }
        if ( ( strInputData == "" ) && ( strOutputData != "" ) )
        {
            _EXCEPTIONT ( "--out_data specified without --in_data" );
        }

        // Check command line parameters (data type arguments)
        STLStringHelper::ToLower ( strInputType );
        STLStringHelper::ToLower ( strOutputType );

        DiscretizationType eInputType;
        DiscretizationType eOutputType;

        if ( strInputType == "fv" )
        {
            eInputType = DiscretizationType_FV;
        }
        else if ( strInputType == "cgll" )
        {
            eInputType = DiscretizationType_CGLL;
        }
        else if ( strInputType == "dgll" )
        {
            eInputType = DiscretizationType_DGLL;
        }
        else
        {
            _EXCEPTION1 ( "Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
                          strInputType.c_str() );
        }

        if ( strOutputType == "fv" )
        {
            eOutputType = DiscretizationType_FV;
        }
        else if ( strOutputType == "cgll" )
        {
            eOutputType = DiscretizationType_CGLL;
        }
        else if ( strOutputType == "dgll" )
        {
            eOutputType = DiscretizationType_DGLL;
        }
        else
        {
            _EXCEPTION1 ( "Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
                          strOutputType.c_str() );
        }

        // Monotonicity flags
        int nMonotoneType = fMonotoneTypeID;

        // Parse variable list
        std::vector< std::string > vecVariableStrings;
        ParseVariableList ( strVariables, vecVariableStrings );

        // Parse preserve variable list
        std::vector< std::string > vecPreserveVariableStrings;
        ParseVariableList ( strPreserveVariables, vecPreserveVariableStrings );

        if ( fPreserveAll && ( vecPreserveVariableStrings.size() != 0 ) )
        {
            _EXCEPTIONT ( "--preserveall and --preserve cannot both be specified" );
        }

        // Calculate Face areas
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating input mesh Face areas\n" );
        double dTotalAreaInput_loc = m_meshInput->CalculateFaceAreas(fInputConcave);
        Real dTotalAreaInput;
        MPI_Allreduce ( &dTotalAreaInput_loc, &dTotalAreaInput, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput );

        // Input mesh areas
        if ( eInputType == DiscretizationType_FV )
        {
            this->SetSourceAreas ( m_meshInputCov->vecFaceArea );
        }

        // Calculate Face areas
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating output mesh Face areas\n" );
        Real dTotalAreaOutput_loc = m_meshOutput->CalculateFaceAreas(fOutputConcave);
        Real dTotalAreaOutput;
        MPI_Allreduce ( &dTotalAreaOutput_loc, &dTotalAreaOutput, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Output Mesh Geometric Area: %1.15e\n", dTotalAreaOutput );

        // Output mesh areas
        if ( eOutputType == DiscretizationType_FV )
        {
            this->SetTargetAreas ( m_meshOutput->vecFaceArea );
        }

        // Verify that overlap mesh is in the correct order
        int ixSourceFaceMax = ( -1 );
        int ixTargetFaceMax = ( -1 );

        if ( m_meshOverlap->vecSourceFaceIx.size() !=
                m_meshOverlap->vecTargetFaceIx.size()
           )
        {
            _EXCEPTIONT ( "Invalid overlap mesh:\n"
                          "    Possible mesh file corruption?" );
        }

        for ( unsigned i = 0; i < m_meshOverlap->vecSourceFaceIx.size(); i++ )
        {
            if ( m_meshOverlap->vecSourceFaceIx[i] + 1 > ixSourceFaceMax )
            {
                ixSourceFaceMax = m_meshOverlap->vecSourceFaceIx[i] + 1;
            }
            if ( m_meshOverlap->vecTargetFaceIx[i] + 1 > ixTargetFaceMax )
            {
                ixTargetFaceMax = m_meshOverlap->vecTargetFaceIx[i] + 1;
            }
        }

        /*
        // Check for forward correspondence in overlap mesh
        if ( // m_meshInputCov->faces.size() - ixSourceFaceMax == 0 //&&
            ( m_meshOutput->faces.size() - ixTargetFaceMax == 0 )
        )
        {
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Overlap mesh forward correspondence found\n" );
        }
        else if (
            // m_meshOutput->faces.size() - ixSourceFaceMax == 0 //&&
            ( m_meshInputCov->faces.size() - ixTargetFaceMax == 0 )
        )
        {   // Check for reverse correspondence in overlap mesh
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Overlap mesh reverse correspondence found (reversing)\n" );

            // Reorder overlap mesh
            m_meshOverlap->ExchangeFirstAndSecondMesh();
        }
        else
        {   // No correspondence found
            _EXCEPTION4 ( "Invalid overlap mesh:\n"
                          "    No correspondence found with input and output meshes (%i,%i) vs (%i,%i)",
                          m_meshInputCov->faces.size(), m_meshOutput->faces.size(), ixSourceFaceMax, ixTargetFaceMax );
        }
        */

        // Calculate Face areas
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating overlap mesh Face areas\n" );
        Real dTotalAreaOverlap_loc = m_meshOverlap->CalculateFaceAreas(false);
        Real dTotalAreaOverlap;
        MPI_Allreduce ( &dTotalAreaOverlap_loc, &dTotalAreaOverlap, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
        if ( !pcomm->rank() ) dbgprint.printf ( 0, "Overlap Mesh Area: %1.15e\n", dTotalAreaOverlap );

        // Partial cover
        if ( fabs ( dTotalAreaOverlap - dTotalAreaInput ) > 1.0e-10 )
        {
            if ( !fNoCheck )
            {
                if ( !pcomm->rank() ) dbgprint.printf ( 0, "WARNING: Significant mismatch between overlap mesh area "
                                                            "and input mesh area.\n  Automatically enabling --nocheck\n" );
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
        if ( ( eInputType  == DiscretizationType_FV ) &&
                ( eOutputType == DiscretizationType_FV )
           )
        {
            // Generate reverse node array and edge map
            m_meshInputCov->ConstructReverseNodeArray();
            m_meshInputCov->ConstructEdgeMap();

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFV ( *m_meshInputCov );
            this->InitializeTargetCoordinatesFromMeshFV ( *m_meshOutput );

            // Construct OfflineMap
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating offline map\n" );
            LinearRemapFVtoFV_Tempest_MOAB ( nPin );

            // Finite volume input / Finite element output
        }
        else if ( eInputType == DiscretizationType_FV )
        {
            DataMatrix3D<int> dataGLLNodes;
            DataMatrix3D<double> dataGLLJacobian;

            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Generating output mesh meta data\n" );
            double dNumericalArea_loc =
                GenerateMetaData (
                    *m_meshOutput,
                    nPout,
                    fBubble,
                    dataGLLNodes,
                    dataGLLJacobian );

            Real dNumericalArea;
            MPI_Allreduce ( &dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Output Mesh Numerical Area: %1.15e\n", dNumericalArea );

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFV ( *m_meshInputCov );
            this->InitializeTargetCoordinatesFromMeshFE (
                *m_meshOutput, nPout, dataGLLNodes );

            // Generate the continuous Jacobian
            bool fContinuous = ( eOutputType == DiscretizationType_CGLL );

            if ( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodes,
                    dataGLLJacobian,
                    this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobian,
                    this->GetTargetAreas() );
            }

            // Generate reverse node array and edge map
            m_meshInputCov->ConstructReverseNodeArray();
            m_meshInputCov->ConstructEdgeMap();

            // Generate remap weights
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating offline map\n" );

            if ( fVolumetric )
            {
                LinearRemapFVtoGLL_Volumetric_MOAB (
                    dataGLLNodes,
                    dataGLLJacobian,
                    this->GetTargetAreas(),
                    nPin,
                    nMonotoneType,
                    fContinuous,
                    fNoConservation );
            }
            else
            {
                LinearRemapFVtoGLL_MOAB (
                    dataGLLNodes,
                    dataGLLJacobian,
                    this->GetTargetAreas(),
                    nPin,
                    nMonotoneType,
                    fContinuous,
                    fNoConservation );
            }

            // Finite element input / Finite volume output
        }
        else if (
            ( eInputType != DiscretizationType_FV ) &&
            ( eOutputType == DiscretizationType_FV )
        )
        {
            DataMatrix3D<int> dataGLLNodes;
            DataMatrix3D<double> dataGLLJacobian;

            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Generating input mesh meta data\n" );
            double dNumericalArea_loc =
                GenerateMetaData (
                    *m_meshInputCov,
                    nPin,
                    fBubble,
                    dataGLLNodes,
                    dataGLLJacobian );

            Real dNumericalArea;
            MPI_Allreduce ( &dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Input Mesh Numerical Area: %1.15e\n", dNumericalArea );

            if ( fabs ( dNumericalArea - dTotalAreaInput ) > 1.0e-12 )
            {
                dbgprint.printf ( 0, "WARNING: Significant mismatch between input mesh "
                                  "numerical area and geometric area\n" );
            }

            if ( dataGLLNodes.GetSubColumns() != m_meshInputCov->faces.size() )
            {
                _EXCEPTIONT ( "Number of element does not match between metadata and "
                              "input mesh" );
            }

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInputCov, nPin, dataGLLNodes );
            this->InitializeTargetCoordinatesFromMeshFV ( *m_meshOutput );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if ( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodes,
                    dataGLLJacobian,
                    this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobian,
                    this->GetSourceAreas() );
            }

            // Generate offline map
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating offline map\n" );

            if ( fVolumetric )
            {
                _EXCEPTIONT ( "Unimplemented: Volumetric currently unavailable for"
                              "GLL input mesh" );
            }

            LinearRemapSE4_Tempest_MOAB (
                dataGLLNodes,
                dataGLLJacobian,
                nMonotoneType,
                fContinuousIn,
                fNoConservation
            );

            // Finite element input / Finite element output
        }
        else if (
            ( eInputType  != DiscretizationType_FV ) &&
            ( eOutputType != DiscretizationType_FV )
        )
        {
            DataMatrix3D<int> dataGLLNodesIn;
            DataMatrix3D<double> dataGLLJacobianIn;

            DataMatrix3D<int> dataGLLNodesOut;
            DataMatrix3D<double> dataGLLJacobianOut;

            // Input metadata
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Generating input mesh meta data" );
            double dNumericalAreaIn_loc =
                GenerateMetaData (
                    *m_meshInputCov,
                    nPin,
                    fBubble,
                    dataGLLNodesIn,
                    dataGLLJacobianIn );

            Real dNumericalAreaIn;
            MPI_Allreduce ( &dNumericalAreaIn_loc, &dNumericalAreaIn, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Input Mesh Numerical Area: %1.15e", dNumericalAreaIn );

            if ( fabs ( dNumericalAreaIn - dTotalAreaInput ) > 1.0e-12 )
            {
                dbgprint.printf ( 0, "WARNING: Significant mismatch between input mesh "
                                  "numerical area and geometric area" );
            }

            // Output metadata
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Generating output mesh meta data" );
            double dNumericalAreaOut_loc =
                GenerateMetaData (
                    *m_meshOutput,
                    nPout,
                    fBubble,
                    dataGLLNodesOut,
                    dataGLLJacobianOut );

            Real dNumericalAreaOut;
            MPI_Allreduce ( &dNumericalAreaOut_loc, &dNumericalAreaOut, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );

            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Output Mesh Numerical Area: %1.15e", dNumericalAreaOut );

            if ( fabs ( dNumericalAreaOut - dTotalAreaOutput ) > 1.0e-12 )
            {
                if ( !pcomm->rank() ) dbgprint.printf ( 0, "WARNING: Significant mismatch between output mesh "
                                                            "numerical area and geometric area" );
            }

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInputCov, nPin, dataGLLNodesIn );
            this->InitializeTargetCoordinatesFromMeshFE (
                *m_meshOutput, nPout, dataGLLNodesOut );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if ( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesIn,
                    dataGLLJacobianIn,
                    this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobianIn,
                    this->GetSourceAreas() );
            }

            // Generate the continuous Jacobian for output mesh
            bool fContinuousOut = ( eOutputType == DiscretizationType_CGLL );

            if ( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesOut,
                    dataGLLJacobianOut,
                    this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobianOut,
                    this->GetTargetAreas() );
            }

            // Generate offline map
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Calculating offline map" );

            LinearRemapGLLtoGLL2_MOAB (
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
                fNoConservation
            );

        }
        else
        {
            _EXCEPTIONT ( "Not implemented" );
        }

//#pragma warning "NOTE: VERIFICATION DISABLED"
        // Verify consistency, conservation and monotonicity
        if ( !fNoCheck )
        {
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Verifying map" );
            this->IsConsistent ( 1.0e-8 );
            if ( !fNoConservation ) this->IsConservative ( 1.0e-8 );

            if ( nMonotoneType != 0 )
            {
                this->IsMonotone ( 1.0e-12 );
            }
        }

        // Apply Offline Map to data
        if ( strInputData != "" )
        {
            if ( !pcomm->rank() ) dbgprint.printf ( 0, "Applying offline map to data\n" );

            this->SetFillValueOverride ( static_cast<float> ( dFillValueOverride ) );
            this->Apply (
                strInputData,
                strOutputData,
                vecVariableStrings,
                strNColName,
                fOutputDouble,
                false );
        }

        // Copy variables from input file to output file
        if ( ( strInputData != "" ) && ( strOutputData != "" ) )
        {
            if ( fPreserveAll )
            {
                if ( !pcomm->rank() ) dbgprint.printf ( 0, "Preserving variables" );
                this->PreserveAllVariables ( strInputData, strOutputData );

            }
            else if ( vecPreserveVariableStrings.size() != 0 )
            {
                if ( !pcomm->rank() ) dbgprint.printf ( 0, "Preserving variables" );
                this->PreserveVariables (
                    strInputData,
                    strOutputData,
                    vecPreserveVariableStrings );
            }
        }

    }
    catch ( Exception & e )
    {
        dbgprint.printf ( 0, "%s", e.ToString().c_str() );
        return ( moab::MB_FAILURE );

    }
    catch ( ... )
    {
        return ( moab::MB_FAILURE );
    }
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsConsistent (
    double dTolerance
)
{
    // Get map entries
    DataVector<int> dataRows;
    DataVector<int> dataCols;
    DataVector<double> dataEntries;

    // Calculate row sums
    DataVector<double> dRowSums;
    if ( pcomm->size() > 1 )
    {
        if ( !m_globalMapAvailable )
            this->GatherAllToRoot();
        if ( pcomm->size() > 1 && pcomm->rank() ) return true;
        SparseMatrix<double>& m_mapRemapGlobal = m_weightMapGlobal->GetSparseMatrix();
        m_mapRemapGlobal.GetEntries ( dataRows, dataCols, dataEntries );
        dRowSums.Initialize ( m_mapRemapGlobal.GetRows() );
    }
    else
    {
        m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );
        dRowSums.Initialize ( m_mapRemap.GetRows() );
    }

    for ( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        dRowSums[dataRows[i]] += dataEntries[i];
    }

    // Verify all row sums are equal to 1
    bool fConsistent = true;
    for ( unsigned i = 0; i < dRowSums.GetRows(); i++ )
    {
        if ( fabs ( dRowSums[i] - 1.0 ) > dTolerance )
        {
            fConsistent = false;
            Announce ( "TempestOfflineMap is not consistent in row %i (%1.15e)",
                       i, dRowSums[i] );
        }
    }

    return fConsistent;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsConservative (
    double dTolerance
)
{
    // Get map entries
    DataVector<int> dataRows;
    DataVector<int> dataCols;
    DataVector<double> dataEntries;
    const DataVector<double>& dTargetAreas = this->GetGlobalTargetAreas();
    const DataVector<double>& dSourceAreas = this->GetGlobalSourceAreas();

    // Calculate column sums
    DataVector<double> dColumnSums;
    if ( pcomm->size() > 1 )
    {
        if ( !m_globalMapAvailable )
            this->GatherAllToRoot();
        if ( pcomm->size() > 1 && pcomm->rank() ) return true;
        SparseMatrix<double>& m_mapRemapGlobal = m_weightMapGlobal->GetSparseMatrix();
        m_mapRemapGlobal.GetEntries ( dataRows, dataCols, dataEntries );
        dColumnSums.Initialize ( m_mapRemapGlobal.GetColumns() );
    }
    else
    {
        m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );
        dColumnSums.Initialize ( m_mapRemap.GetColumns() );
    }

    for ( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        dColumnSums[dataCols[i]] +=
            dataEntries[i] * dTargetAreas[dataRows[i]];
    }

    // Verify all column sums equal the input Jacobian
    bool fConservative = true;
    for ( unsigned i = 0; i < dColumnSums.GetRows(); i++ )
    {
        if ( fabs ( dColumnSums[i] - dSourceAreas[i] ) > dTolerance )
        {
            fConservative = false;
            Announce ( "TempestOfflineMap is not conservative in column "
                       "%i (%1.15e / %1.15e)",
                       i, dColumnSums[i], dSourceAreas[i] );
        }
    }

    return fConservative;
}

///////////////////////////////////////////////////////////////////////////////

bool moab::TempestOfflineMap::IsMonotone (
    double dTolerance
)
{
    // Get map entries
    DataVector<int> dataRows;
    DataVector<int> dataCols;
    DataVector<double> dataEntries;

    if ( pcomm->size() > 1 )
    {
        if ( !m_globalMapAvailable )
            this->GatherAllToRoot();
        if ( pcomm->size() > 1 && pcomm->rank() ) return true;
        SparseMatrix<double>& m_mapRemapGlobal = m_weightMapGlobal->GetSparseMatrix();
        m_mapRemapGlobal.GetEntries ( dataRows, dataCols, dataEntries );
    }
    else
        m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );

    // Verify all entries are in the range [0,1]
    bool fMonotone = true;
    for ( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        if ( ( dataEntries[i] < -dTolerance ) ||
                ( dataEntries[i] > 1.0 + dTolerance )
           )
        {
            fMonotone = false;

            Announce ( "TempestOfflineMap is not monotone in entry (%i): %1.15e",
                       i, dataEntries[i] );
        }
    }

    return fMonotone;
}


///////////////////////////////////////////////////////////////////////////////
moab::ErrorCode moab::TempestOfflineMap::GatherAllToRoot()   // Collective
{
    Mesh globalMesh;
    int ierr, rootProc = 0;
    moab::ErrorCode rval;

    // Write SparseMatrix entries
    DataVector<int> vecRow;
    DataVector<int> vecCol;
    DataVector<double> vecS;

    moab::DebugOutput dbgprint ( std::cout, ( pcomm ? pcomm->rank() : 0 ) );

    m_mapRemap.GetEntries ( vecRow, vecCol, vecS );
    const DataVector<double>& dSourceAreas = m_meshInputCov->vecFaceArea;
    const DataVector<double>& dTargetAreas = m_meshOutput->vecFaceArea;

    // Translate the index in Row and Col to global_id and dump it out

    // Communicate the necessary data
    int grows = 0, gcols = 0, gsrc = 0, gtar = 0;
    std::vector<int> rowcolss, rowcolsv;
    DataVector<int> rows, cols, srcelmindx, tgtelmindx;
    {
        // First, accumulate the sizes of rows and columns of the matrix
        if ( !pcomm->rank() ) rowcolss.resize ( pcomm->size() * 4 );

        int sendarray[4];
        sendarray[0] = vecRow.GetRows();
        sendarray[1] = vecCol.GetRows();
        sendarray[2] = dSourceAreas.GetRows();
        sendarray[3] = dTargetAreas.GetRows();

        ierr = MPI_Gather ( sendarray, 4, MPI_INTEGER, rowcolss.data(), 4, MPI_INTEGER, rootProc, pcomm->comm() );
        if ( ierr != MPI_SUCCESS ) return moab::MB_FAILURE;

        if ( !pcomm->rank() )
        {
            for ( unsigned i = 0; i < pcomm->size(); ++i )
            {
                grows += rowcolss[4 * i];
                gcols += rowcolss[4 * i + 1];
                gsrc += rowcolss[4 * i + 2];
                gtar += rowcolss[4 * i + 3];
            }
            rowcolsv.resize ( grows + gcols + gsrc + gtar );
            rows.Initialize ( grows ); // we are assuming rows = cols
            cols.Initialize ( gcols ); // we are assuming rows = cols

            // Let us allocate our global offline map object
            m_weightMapGlobal = new OfflineMap();
            // m_weightMapGlobal->InitializeSourceDimensionsFromMesh(*m_meshInput);
            m_weightMapGlobal->InitializeSourceDimensionsFromMesh ( gsrc, 0 );
            // m_weightMapGlobal->InitializeTargetDimensionsFromMesh(*m_meshOutput);
            m_weightMapGlobal->InitializeTargetDimensionsFromMesh ( gtar, 0 );

            m_weightMapGlobal->InitializeSourceCoordinatesFromMeshFV ( *m_meshInput );
            m_weightMapGlobal->InitializeTargetCoordinatesFromMeshFV ( *m_meshOutput );

            m_weightMapGlobal->GetSourceAreas().Initialize ( gsrc ); srcelmindx.Initialize ( gsrc );
            m_weightMapGlobal->GetTargetAreas().Initialize ( gtar ); tgtelmindx.Initialize ( gtar );

#ifdef VERBOSE
            dbgprint.printf ( 0, "Received global dimensions: %d, %d\n", vecRow.GetRows(), rows.GetRows() );
            dbgprint.printf ( 0, "Global: n(source) = %d, and n(target) = %d\n", gsrc, gtar );
            dbgprint.printf ( 0, "Operator size = %d\n", grows + gcols );
#endif
        }
    }

#ifdef VERBOSE
    {
        std::stringstream sstr;
        sstr << "rowscols_" << pcomm->rank() << ".txt";
        std::ofstream output_file ( sstr.str().c_str() );
        output_file << "VALUES\n";
        for ( unsigned ip = 0; ip < vecRow.GetRows(); ++ip )
        {
            output_file << ip << " (" << vecRow[ip] << ", " << vecCol[ip] << ") = " << vecS[ip] << "\n";
        }
        output_file.flush(); // required here
        output_file.close();
    }
#endif

    {
        // Next, accumulate the row and column values for the matrices (in local indexing)
        const int nR = vecRow.GetRows() + vecCol.GetRows() + dSourceAreas.GetRows() + dTargetAreas.GetRows();
        std::vector<int> sendarray ( nR );
        for ( unsigned ix = 0; ix < vecRow.GetRows(); ++ix )
        {
            sendarray[ix] = m_remapper->GetGlobalID ( moab::Remapper::TargetMesh, vecRow[ix] );
        }
        for ( unsigned ix = 0, offset = vecRow.GetRows(); ix < vecCol.GetRows(); ++ix )
        {
            sendarray[offset + ix] = m_remapper->GetGlobalID ( moab::Remapper::CoveringMesh, vecCol[ix] );
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
            rval = mbCore->tag_get_handle ( "GLOBAL_ID", gidtag ); MB_CHK_ERR ( rval );
            rval = mbCore->tag_get_data ( gidtag, m_remapper->m_covering_source_entities, &sendarray[vecRow.GetRows() + vecCol.GetRows()] ); MB_CHK_ERR ( rval );
            rval = mbCore->tag_get_data ( gidtag, m_remapper->m_target_entities, &sendarray[vecRow.GetRows() + vecCol.GetRows() + dSourceAreas.GetRows()] ); MB_CHK_ERR ( rval );
        }

        std::vector<int> displs, rcount;
        if ( !pcomm->rank() )
        {
            displs.resize ( pcomm->size(), 0 );
            rcount.resize ( pcomm->size(), 0 );
            int gsum = 0;
            for ( unsigned i = 0; i < pcomm->size(); ++i )
            {
                displs[i] = gsum;
                rcount[i] = rowcolss[4 * i] + rowcolss[4 * i + 1] + rowcolss[4 * i + 2] + rowcolss[4 * i + 3];
                gsum += rcount[i];
            }
            assert ( rowcolsv.size() - gsum == 0 );
        }

        // Both rows and columns have a size of "rowsize"
        ierr = MPI_Gatherv ( &sendarray[0], nR, MPI_INTEGER, &rowcolsv[0], &rcount[0], &displs[0], MPI_INTEGER, rootProc, pcomm->comm() );

        if ( !pcomm->rank() )
        {
#ifdef VERBOSE
            std::ofstream output_file ( "rows-cols.txt", std::ios::out );
            output_file << "ROWS\n";
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip], iend = istart + rowcolss[4 * ip];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
#ifdef VERBOSE
                    output_file << offset << " " << rowcolsv[i] << "\n";
#endif
                    rows[offset] = rowcolsv[i];
                }
            }
#ifdef VERBOSE
            output_file << "COLS\n";
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip] + rowcolss[4 * ip], iend = istart + rowcolss[4 * ip + 1];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
#ifdef VERBOSE
                    output_file << offset << " " << rowcolsv[i] << "\n";
#endif
                    cols[offset] = rowcolsv[i];
                }
            }
#ifdef VERBOSE
            output_file.flush(); // required here
            output_file.close();
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip] + rowcolss[4 * ip] + rowcolss[4 * ip + 1], iend = istart + rowcolss[4 * ip + 2];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
                    srcelmindx[offset] = rowcolsv[i];
                }
            }
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip] + rowcolss[4 * ip] + rowcolss[4 * ip + 1] + rowcolss[4 * ip + 2], iend = istart + rowcolss[4 * ip + 3];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
                    tgtelmindx[offset] = rowcolsv[i];
                }
            }
        } /* (!pcomm->rank()) */
        rowcolsv.clear();
    }

    {
        // Next, accumulate the row and column values for the matrices (in local indexing)
        const int nR = vecS.GetRows() + dSourceAreas.GetRows() + dTargetAreas.GetRows();

        std::vector<double> sendarray ( nR );
        std::copy ( ( double* ) vecS, ( double* ) vecS + vecS.GetRows(), sendarray.begin() );
        std::copy ( ( const double* ) dSourceAreas, ( const double* ) dSourceAreas + dSourceAreas.GetRows(), sendarray.begin() + vecS.GetRows() );
        std::copy ( ( const double* ) dTargetAreas, ( const double* ) dTargetAreas + dTargetAreas.GetRows(), sendarray.begin() + vecS.GetRows() + dSourceAreas.GetRows() );

        std::vector<int> displs, rcount;
        int gsum = 0;
        if ( !pcomm->rank() )
        {
            displs.resize ( pcomm->size(), 0 );
            rcount.resize ( pcomm->size(), 0 );
            for ( unsigned i = 0; i < pcomm->size(); ++i )
            {
                displs[i] = gsum;
                rcount[i] = rowcolss[4 * i] + rowcolss[4 * i + 2] + rowcolss[4 * i + 3];
                gsum += rcount[i];
            }
        }

        std::vector<double> rowcolsvals ( grows + gsrc + gtar );
        // Both rows and columns have a size of "rowsize"
        ierr = MPI_Gatherv ( sendarray.data(), sendarray.size(), MPI_DOUBLE, rowcolsvals.data(), &rcount[0], &displs[0], MPI_DOUBLE, rootProc, pcomm->comm() );

        if ( !pcomm->rank() )
        {
            DataVector<double> globvalues ( grows );
#ifdef VERBOSE
            std::ofstream output_file ( "rows-cols.txt", std::ios::app );
            output_file << "VALUES\n";
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                for ( int i = displs[ip]; i < displs[ip] + rowcolss[4 * ip]; ++i, ++offset )
                {
                    globvalues[offset] = rowcolsvals[i];
#ifdef VERBOSE
                    output_file << offset << " (" << rows[offset] << ", " << cols[offset] << ") = " << rowcolsvals[i] << "\n";
#endif
                }
            }
#ifdef VERBOSE
            output_file.flush(); // required here
            output_file.close();
#endif
            // m_mapRemapGlobal.SetEntries(rows, cols, globvalues);
            m_weightMapGlobal->GetSparseMatrix().SetEntries ( rows, cols, globvalues );
        }

        if ( !pcomm->rank() )
        {
            DataVector<double>& m_areasSrcGlobal = m_weightMapGlobal->GetSourceAreas();
            DataVector<double>& m_areasTgtGlobal = m_weightMapGlobal->GetTargetAreas();
            // Store the global source and target elements areas
#ifdef VERBOSE
            std::ofstream output_file ( "source-target-areas.txt" );
            output_file << "Source areas\n";
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip] + rowcolss[4 * ip], iend = istart + rowcolss[4 * ip + 2];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
#ifdef VERBOSE
                    output_file << srcelmindx[offset] << " " << rowcolsvals[i] << "\n";
#endif
                    m_areasSrcGlobal[srcelmindx[offset]] = rowcolsvals[i];
                }
            }
#ifdef VERBOSE
            output_file << "Target areas\n";
#endif
            for ( unsigned ip = 0, offset = 0; ip < pcomm->size(); ++ip )
            {
                int istart = displs[ip] + rowcolss[4 * ip] + rowcolss[4 * ip + 2], iend = istart + rowcolss[4 * ip + 3];
                for ( int i = istart; i < iend; ++i, ++offset )
                {
#ifdef VERBOSE
                    output_file << tgtelmindx[offset] << " " << rowcolsvals[i] << "\n";
#endif
                    m_areasTgtGlobal[tgtelmindx[offset]] = rowcolsvals[i];
                }
            }
#ifdef VERBOSE
            output_file.flush(); // required here
            output_file.close();
#endif
        }

    }
    rowcolss.clear();

    // Update all processes that we have a global view of the map available
    // on the root process
    m_globalMapAvailable = true;

    if ( !pcomm->rank() && false )
    {
        dbgprint.printf ( 0, "Writing out file outGlobalView.nc\n" );
        // m_dSourceCenterLon
        // m_dSourceCenterLat
        // m_dTargetCenterLon
        // m_dTargetCenterLat
        // m_dSourceVertexLon
        // m_dSourceVertexLat
        m_weightMapGlobal->Write ( "outGlobalView.nc" );
    }

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
