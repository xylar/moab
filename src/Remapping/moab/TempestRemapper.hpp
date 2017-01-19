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

class TempestRemapper : public Remapper
{
public:

    TempestRemapper(moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL) :
        Remapper(mbInt, pcomm), meshValidate(false), constructEdgeMap(false), m_source_type(DEFAULT), m_target_type(DEFAULT)
    {
    }

    virtual ErrorCode initialize();

    // Mesh type with a correspondence to Tempest/Climate formats
    enum TempestMeshType { 
        DEFAULT = -1,
        CS = 0,
        RLL = 1,
        ICO = 2,
        OVERLAP_FILES = 3,
        OVERLAP_MEMORY = 4,
        OVERLAP_MOAB = 5
    };

    moab::ErrorCode GenerateMesh(Remapper::IntersectionContext ctx, TempestMeshType type);

    moab::ErrorCode LoadMesh(Remapper::IntersectionContext ctx, std::string inputFilename, TempestMeshType type);

    moab::ErrorCode ComputeOverlapMesh(double tolerance=1e-8, bool use_tempest=false);

    // Converters between MOAB and Tempest representations
    moab::ErrorCode ConvertTempestMesh(Remapper::IntersectionContext ctx);

    moab::ErrorCode ConvertMeshToTempest(Remapper::IntersectionContext ctx);

    // Parallel utilities
    moab::ErrorCode ExchangeGhostWeights(OfflineMap* weightMap);

    Mesh* GetMesh(Remapper::IntersectionContext ctx);

    void SetMesh(Remapper::IntersectionContext ctx, Mesh* mesh, bool overwrite=true);

    moab::EntityHandle GetMeshSet(Remapper::IntersectionContext ctx) const;

    void SetMeshType(Remapper::IntersectionContext ctx, TempestMeshType type);

    TempestMeshType GetMeshType(Remapper::IntersectionContext ctx) const;

    // public members
    bool meshValidate;  // Validate the mesh after loading from file

    bool constructEdgeMap;  //  Construct the edge map within the TempestRemap datastructures

    static const bool verbose = true;

private:

    // private methods
    moab::ErrorCode LoadTempestMesh_Private(std::string inputFilename, Mesh** tempest_mesh);

    moab::ErrorCode ConvertMOABMeshToTempest_Private(Mesh* mesh, moab::EntityHandle meshset);

    moab::ErrorCode ConvertTempestMeshToMOAB_Private(TempestMeshType type, Mesh* mesh, moab::EntityHandle& meshset);

    moab::ErrorCode AssociateSrcTargetInOverlap(Mesh* mesh, EntityHandle* meshsets);

    // Source and Target meshes
    Mesh* m_source;
    TempestMeshType m_source_type;
    moab::EntityHandle m_source_set;

    Mesh* m_target;
    TempestMeshType m_target_type;
    moab::EntityHandle m_target_set;

    // Overlap meshes
    Mesh* m_overlap;
    TempestMeshType m_overlap_type;
    moab::EntityHandle m_overlap_set;

};

// Inline functions
inline
Mesh* TempestRemapper::GetMesh(Remapper::IntersectionContext ctx)
{
    switch(ctx)
    {
        case Remapper::SourceMesh:
            return m_source;
        case Remapper::TargetMesh:
            return m_target;
        case Remapper::IntersectedMesh:
            return m_overlap;
        default:
            // Nothing to return.
            return NULL;
    }
}

inline
void TempestRemapper::SetMesh(Remapper::IntersectionContext ctx, Mesh* mesh, bool overwrite)
{
    switch(ctx)
    {
        case Remapper::SourceMesh:
            if (!overwrite && m_source) return;
            if (overwrite && m_source) delete m_source;
            m_source = mesh;
            break;
        case Remapper::TargetMesh:
            if (!overwrite && m_target) return;
            if (overwrite && m_target) delete m_source;
            m_target = mesh;
            break;
        case Remapper::IntersectedMesh:
            if (!overwrite && m_overlap) return;
            if (overwrite && m_overlap) delete m_source;
            m_overlap = mesh;
            break;
        case Remapper::DEFAULT:
            // Nothing to do.
            break;
    }
}

inline
moab::EntityHandle TempestRemapper::GetMeshSet(Remapper::IntersectionContext ctx) const
{
    switch(ctx)
    {
        case Remapper::SourceMesh:
            return m_source_set;
        case Remapper::TargetMesh:
            return m_target_set;
        case Remapper::IntersectedMesh:
            return m_overlap_set;
        default:
            // Nothing to return, so return the rootset
            return 0;
    }
}

inline
void TempestRemapper::SetMeshType(Remapper::IntersectionContext ctx, TempestRemapper::TempestMeshType type)
{
    switch(ctx)
    {
        case Remapper::SourceMesh:
            m_source_type = type;
            break;
        case Remapper::TargetMesh:
            m_target_type = type;
            break;
        case Remapper::IntersectedMesh:
            m_overlap_type = type;
            break;
        case Remapper::DEFAULT:
            // Nothing to do.
            break;
    }
}

inline
TempestRemapper::TempestMeshType TempestRemapper::GetMeshType(Remapper::IntersectionContext ctx) const
{
    switch(ctx)
    {
        case Remapper::SourceMesh:
            return m_source_type;
        case Remapper::TargetMesh:
            return m_target_type;
        case Remapper::IntersectedMesh:
            return m_overlap_type;
        case Remapper::DEFAULT:
            // Nothing to do.
            return TempestRemapper::DEFAULT;
    }
}

}

#endif // MB_TEMPESTREMAPPER_HPP
