# MOAB: Mesh-Oriented datABase

MOAB is a component for representing and evaluating mesh data. MOAB can store structured and unstructured mesh, consisting of elements in the finite element "zoo". The functional interface to MOAB is simple yet powerful, allowing the representation of many types of metadata commonly found on the mesh. MOAB is optimized for efficiency in space and time, based on access to mesh in chunks rather than through individual entities, while also versatile enough to support individual entity access.
MOAB can be used in several ways: as the underlying mesh data representation for applications (MOAB is used in this way in the VERDE mesh verification code, as a mesh input mechanism (using mesh readers included with MOAB), or as a translator between mesh formats (using readers and writers included with MOAB).

MOAB was developed as part of the CUBIT project at Sandia National Laboratories, and was partially funded by the DOE SciDAC program as part of the Terascale Tools and Technologies (TSTT) center.

# Dependencies
MOAB depends on the NetCDF libraries (C and C++) to compile the ExodusII reader/writer. Support for C++ was added to netcdf in version 3.5.1, and took a bit of time to get compiling ironed out, so make sure you have version 3.6 or later. To get netcdf, search the web or try [NetCDF].
The only MOAB file format that can represent the entire MOAB data model is MOAB's native HDF5-based file format. Support for this file format requires version 5 of the HDF library. It may be obtained at [HDF5].

# Compiling

  - Unpack the source code in a local directory.
  - Run ```autoreconf -fi``` to generate the configure script
  - Run the ```configure --help``` script in the top source directory to see a list of available configure options.
      - Use ```--prefix=INSTALL_DIR``` to specify installation directory
      - Override default compilers with environment or user options: ```CC, CXX, FC, F77```
      - If you have MPI installed, use ```--with-mpi=$MPI_DIR```
      - If you have HDF5 and NetCDF installed, use ```--with-netcdf=$NETCDF_DIR``` and ```--with-hdf5=$HDF5_DIR``` to specify dependencies.
  - Now run the ```configure``` script with desired configuration options either in-source or out-of-source (build) directory.
  - In the build directory, run the following:
      - Compile MOAB library and supported tools: ```make -j4```.
      - Verify configuration and build setup: ```make check```.
  - To install the compiled libraries, headers and tools, run ```make install```.
  - You can now use the ```makefile``` generated under the ```build/examples``` folder and modify it to compile user code dependent on MOAB libraries

# Continuous Integration

There are several hooks to online continuous integration systems, nightly and commit-based Buildbot builds that are constantly run during a development day to check the integrity and robustness of the MOAB library.

# Bugs, Correspondence, Contributing
MOAB is LGPL code, and we encourage users to submit bug reports (and, if desired, fixes) to moab-dev@mcs.anl.gov. Users are encouraged to check [SIGMA-MOAB] documentation pages often for news and updates. Please submit pull requests (PR) with a Bitbucket fork or send us patches that you would like merged upstream.

[NetCDF]: http://www.unidata.ucar.edu/software/netcdf/
[HDF5]: https://www.hdfgroup.org/HDF5/
[SIGMA-MOAB]: http://sigma.mcs.anl.gov/moab-library
[CodeShip]: [![CodeShip Build Status](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master)](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master)
[Drone.io]: [![Drone.io Build Status](https://drone.io/bitbucket.org/fathomteam/moab/status.png)](https://drone.io/bitbucket.org/fathomteam/moab/latest)