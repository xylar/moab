/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

// copied and modified for use with MOAB, from chaco's code/main/main.c file

#include <iostream>
#include <stdio.h>

#include "ComputePartition.hpp"
#include "MBInterface.hpp"

int main(int argc, char* argv[])
{

  ComputePartition cp;

  if (argc < 5) {
    std::cout << "Usage: " << argv[0] << " <mesh_file> <write_out(y/n)> <output_mesh_file> <nprocs>" 
              << std::endl;
    return 1;
  }
  
  std::cout << "                    Chaco 2.0" << std::endl;
  std::cout << "                      (MOAB-based version)" << std::endl;
  std::cout << "          Sandia National Laboratories" << std::endl;

  bool write_output = false;
  if (argv[2][0] == 'y') write_output = true;
  else if (argv[2][0] == 'n') write_output = false;
  else {
    std::cout << "Argument 2 must be 'y' or 'n', not '" << argv[2][0] << "'" << std::endl;
    return 1;
  }
  
  int nprocs = -1;
  std::
  sscanf(argv[4], "%d", &nprocs);

  MBErrorCode result = cp.compute_partition(nprocs, argv[1], write_output, argv[3]);

  if (MB_SUCCESS != result) {
    std::cout << "compute_partition returned non-success error code = " << result << std::endl;
    return 1;
  }

  return 0;
}
