/*
 * =====================================================================================
 *
 *       Filename:  TempestRemapper.hpp
 *
 *    Description:  Interface to the TempestRemap library to enable intersection and
 *                  high-order conservative remapping of climate solution from
 *                  arbitrary resolution of source and target grids on the sphere.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#ifndef MB_TEMPESTREMAPPER_HPP
#define MB_TEMPESTREMAPPER_HPP

#include "moab/Remapper.hpp"
#include "moab/Intx2MeshOnSphere.hpp"
#include "moab/IntxUtils.hpp"

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

namespace moab
{
enum TempestMeshType { CS = 0, RLL = 1, ICO = 2, OVERLAP = 3, OVERLAP_V2 = 4, OVERLAP_MOAB = 5 };

class TempestRemapper : public Remapper
{
public:
    TempestRemapper(moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL) :
        Remapper(mbInt, pcomm), meshValidate(false), constructEdgeMap(false)
    {
    }

    virtual void initialize();

    moab::ErrorCode GenerateMesh(IntersectionContext ctx, TempestMeshType type);

    moab::ErrorCode LoadMesh(IntersectionContext ctx, std::string inputFilename, TempestMeshType type);

    moab::ErrorCode ComputeOverlapMesh();

    moab::ErrorCode ConvertTempestMesh(IntersectionContext ctx, moab::EntityHandle& meshset);

    moab::ErrorCode ConvertMeshToTempest(moab::EntityHandle meshset, Mesh* mesh);

    // public members
    bool meshValidate;  // Validate the mesh after loading from file

    bool constructEdgeMap;  //  Construct the edge map within the TempestRemap datastructures

    Mesh* GetSourceMesh();

    Mesh* GetTargetMesh();

    Mesh* GetOverlapMesh();

    static moab::ErrorCode LoadTempestMesh(std::string inputFilename, Mesh** tempest_mesh, bool meshValidate = false, bool constructEdgeMap = false);

    static moab::ErrorCode ConvertTempestMeshToMOAB(TempestMeshType type, moab::Interface* mb, Mesh* mesh, moab::EntityHandle& meshset);

    static moab::ErrorCode ConvertMOABMeshToTempest(moab::Interface * mb, ParallelComm* pcomm, Mesh * mesh, moab::EntityHandle meshset);

    static moab::ErrorCode AssociateSrcTargetInOverlap(Interface* mb, Mesh* mesh, EntityHandle* meshsets);

    static moab::ErrorCode ExchangeGhostWeights(Interface* mb, OfflineMap* weightMap);

    static const bool verbose = true;

private:

    // Source and Target meshes
    Mesh* m_source;
    TempestMeshType m_source_type;
    moab::EntityHandle m_source_set;

    Mesh* m_target;
    TempestMeshType m_target_type;
    moab::EntityHandle m_target_set;

    // Overlap meshes
    Mesh* m_overlap;
    moab::EntityHandle m_overlap_set;

};

// Inline functions
inline
Mesh* TempestRemapper::GetSourceMesh()
{ return m_source; }

inline
Mesh* TempestRemapper::GetTargetMesh()
{ return m_target; }

inline
Mesh* TempestRemapper::GetOverlapMesh()
{ return m_overlap; }

}

#endif // MB_TEMPESTREMAPPER_HPP
