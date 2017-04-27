/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */


#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include <iostream>
#include "moab/Core.hpp"
#include "moab/Remapping/TempestRemapper.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include "Internals.hpp"
#include "TestRunner.hpp"

using namespace moab;

const static double radius = 1.0;
const double MOAB_PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
const static double surface_area = 4.0 * MOAB_PI * radius * radius;
const static std::string outFilenames[5] = {"outTempestCS.g", "outTempestRLL.g", "outTempestICO.g", "outTempestICOD.g", "outTempestOV.g"};

void test_tempest_cs_create();
void test_tempest_rll_create();
void test_tempest_ico_create();
void test_tempest_mpas_create();
void test_tempest_overlap_combinations();
void test_tempest_to_moab_convert();

int main(int argc, char** argv)
{
  REGISTER_TEST( test_tempest_cs_create );
  REGISTER_TEST( test_tempest_rll_create );
  REGISTER_TEST( test_tempest_ico_create );
  REGISTER_TEST( test_tempest_mpas_create );
  REGISTER_TEST( test_tempest_overlap_combinations );
  REGISTER_TEST( test_tempest_to_moab_convert );

#ifdef MOAB_HAVE_MPI
  MPI_Init( &argc, &argv );
#endif
  int result = RUN_TESTS( argc, argv );
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return result;
}


void test_tempest_cs_create()
{
  NcError error(NcError::verbose_nonfatal);
  const int blockSize = 5;
  const std::string outFilename = outFilenames[0];

  std::cout << "Creating TempestRemap Cubed-Sphere Mesh ...\n";
  Mesh *tempest_mesh = GenerateCSMesh(blockSize, false, outFilename);

  // Compute the surface area of CS mesh
  const double sphere_area = tempest_mesh->CalculateFaceAreas();
  CHECK_REAL_EQUAL(sphere_area, surface_area, 1e-10);

  delete tempest_mesh;
}


void test_tempest_rll_create()
{
  NcError error(NcError::verbose_nonfatal);
  const int blockSize = 5;
  const std::string outFilename = outFilenames[1];

  std::cout << "Creating TempestRemap Latitude-Longitude Mesh ...\n";
  Mesh *tempest_mesh = GenerateRLLMesh(blockSize * 2, blockSize, 0.0, 360.0, -90.0, 90.0, false, outFilename);

  // Compute the surface area of RLL mesh
  const double sphere_area = tempest_mesh->CalculateFaceAreas();
  CHECK_REAL_EQUAL(sphere_area, surface_area, 1e-10);

  delete tempest_mesh;
}


void test_tempest_ico_create()
{
  NcError error(NcError::verbose_nonfatal);
  const int blockSize = 5;
  const bool computeDual = false;
  const std::string outFilename = outFilenames[2];

  std::cout << "Creating TempestRemap Icosahedral Mesh ...\n";
  Mesh *tempest_mesh = GenerateICOMesh(blockSize, computeDual, outFilename);

  // Compute the surface area of ICO mesh
  const double sphere_area = tempest_mesh->CalculateFaceAreas();
  CHECK_REAL_EQUAL(sphere_area, surface_area, 1e-10);

  delete tempest_mesh;
}


void test_tempest_mpas_create()
{
  NcError error(NcError::verbose_nonfatal);
  const int blockSize = 5;
  const bool computeDual = true;
  const std::string outFilename = outFilenames[3];

  std::cout << "Creating TempestRemap MPAS Mesh (dual of the Icosahedral) ...\n";
  Mesh *tempest_mesh = GenerateICOMesh(blockSize, computeDual, outFilename);

  // Compute the surface area of MPAS mesh
  const double sphere_area = tempest_mesh->CalculateFaceAreas();
  CHECK_REAL_EQUAL(sphere_area, surface_area, 1e-10);

  delete tempest_mesh;
}


void test_tempest_overlap_combinations()
{
  NcError error(NcError::verbose_nonfatal);
  const std::string outFilename = outFilenames[4];

  Mesh* inpMesh = new Mesh(outFilenames[0]);
  // verify input mesh area first
  const double inpArea = inpMesh->CalculateFaceAreas();
  CHECK_REAL_EQUAL(inpArea, surface_area, 1e-10);

  for (int isrc=0; isrc < 4; ++isrc) {
    for (int jsrc=0; jsrc < 4; ++jsrc) {
      std::cout << "Computing Overlap between " << outFilenames[isrc] << " and " << outFilenames[jsrc] << " ...\n";
      Mesh *tempest_mesh = GenerateOverlapMesh(outFilenames[isrc], outFilenames[jsrc], outFilename, "exact", false, false);
      // verify overlap mesh area
      const double ovArea = tempest_mesh->CalculateFaceAreas();
      CHECK_REAL_EQUAL(ovArea, surface_area, 1e-10);
      delete tempest_mesh;
    }
  }
  delete inpMesh;
}


void test_tempest_to_moab_convert()
{
  NcError error(NcError::verbose_nonfatal);

  // Allocate and create MOAB Remapper object
  moab::ErrorCode rval;
  moab::Interface* mbCore = new (std::nothrow) moab::Core;
  CHECK (NULL != mbCore);

  moab::ParallelComm* pcomm = new moab::ParallelComm(mbCore, MPI_COMM_WORLD, 0);

  moab::TempestRemapper *remapper = new moab::TempestRemapper(mbCore, pcomm);
  remapper->meshValidate = true;
  remapper->constructEdgeMap = true;
  remapper->initialize();

  rval = pcomm->check_all_shared_handles();CHECK_ERR(rval);

  rval = remapper->LoadMesh(moab::Remapper::SourceMesh, outFilenames[0], moab::TempestRemapper::CS);CHECK_ERR(rval);

  // Load the meshes and validate
  rval = remapper->ConvertTempestMesh(moab::Remapper::SourceMesh);CHECK_ERR(rval);

  Mesh* srcTempest = remapper->GetMesh(moab::Remapper::SourceMesh);

  moab::EntityHandle srcset = remapper->GetMeshSet(moab::Remapper::SourceMesh);

  moab::EntityHandle& tgtset = remapper->GetMeshSet(moab::Remapper::TargetMesh);

  tgtset = srcset;

  // Load the meshes and validate
  rval = remapper->ConvertMeshToTempest(moab::Remapper::TargetMesh);CHECK_ERR(rval);

  Mesh* tgtTempest = remapper->GetMesh(moab::Remapper::TargetMesh);

  const size_t tempest_nodes_src = srcTempest->nodes.size(), tempest_elems_src = srcTempest->faces.size();
  const size_t tempest_nodes_tgt = tgtTempest->nodes.size(), tempest_elems_tgt = tgtTempest->faces.size();
  CHECK_EQUAL(tempest_nodes_src, tempest_nodes_tgt);
  CHECK_EQUAL(tempest_elems_src, tempest_elems_tgt);

  delete remapper;
  delete pcomm;
  delete mbCore;
}

