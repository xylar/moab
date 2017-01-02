/*
 * =====================================================================================
 *
 *       Filename:  Remapper.hpp
 *
 *    Description:  Interface to the a general remapping capability on arbitrary topology
 *                  that performs both mesh intersection between a source and target grid,
 *                  with arbitrary decompositions. The intersections can then be used to
 *                  either evaluate interpolation weights or to perform high-order
 *                  conservative remapping of solutions defined on the source grid.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#ifndef MB_REMAPPER_HPP
#define MB_REMAPPER_HPP

#include "moab/Interface.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

namespace moab
{
enum IntersectionContext { SourceMesh = 0, TargetMesh = 1, IntersectedMesh = 2 };

class Remapper
{
public:
	Remapper(moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL) : m_interface(mbInt), m_pcomm(pcomm)
	{ }

	virtual void initialize() = 0;

protected:

	// member data
	Interface* m_interface;
	ParallelComm* m_pcomm;
};

}

#endif /* MB_REMAPPER_HPP */
