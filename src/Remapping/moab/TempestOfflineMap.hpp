///////////////////////////////////////////////////////////////////////////////
///
///	\file    OfflineMap.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#ifndef _TEMPESTOFFLINEMAP_H_
#define _TEMPESTOFFLINEMAP_H_

#include "SparseMatrix.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "OfflineMap.h"
#include <string>
#include <vector>

#include "moab/TempestRemapper.hpp"

class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class TempestOfflineMap : public OfflineMap {

public:

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	TempestOfflineMap(moab::Interface* mbCore_in, moab::ParallelComm* pcomm_in) : OfflineMap(), mbCore(mbCore_in), pcomm(pcomm_in)
	{ }

	///	<summary>
	///		Define a virtual destructor.
	///	</summary>
	virtual ~TempestOfflineMap()
	{
		mbCore = NULL;
		pcomm = NULL;
	}

	///	<summary>
	///		Initialize the necessary data so that OfflineMap can be generated in parallel.
	///	</summary>
	virtual void Initialize();

public:
	///	<summary>
	///		Generate the offline map, given the source and target mesh and discretization details.
	///     This method generates the mapping between the two meshes based on the overlap and stores 
	///     the result in the SparseMatrix.
	///	</summary>
	moab::ErrorCode GenerateOfflineMap( Mesh& meshInput, Mesh& meshOutput, Mesh& meshOverlap,
                                        std::string strInputType, std::string strOutputType,
                                        int nPin=4, int nPout=4,
                                        bool fBubble=false, int fMonotoneTypeID=0,
                                        bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                        std::string strVariables="", std::string strOutputMap="",
                                        std::string strInputData="", std::string strOutputData="",
                                        std::string strNColName="", bool fOutputDouble=false,
                                        std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0 );

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	// moab::ErrorCode GenerateMetaData();

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	virtual bool IsConsistent(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	virtual bool IsConservative(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	virtual bool IsMonotone(
		double dTolerance
	);

protected:
	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	moab::Interface* mbCore;

	///	<summary>
	///		Vector of areas associated with input degrees of freedom.
	///	</summary>
	moab::ParallelComm* pcomm;

};

///////////////////////////////////////////////////////////////////////////////

#endif

