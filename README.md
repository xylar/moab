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

# Current overall build status
Buildbot: [![Buildbot Status](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master)](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master) \
CodeShip: [![CodeShip Build Status](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master)](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master) \
Drone.io: [![Drone.io Build Status](https://drone.io/bitbucket.org/fathomteam/moab/status.png)](https://drone.io/bitbucket.org/fathomteam/moab/latest)

Where's the configure script?

  After installing the GNU autotools suite, execute the following
  command to generate the configure script (and other necessary
  generated files):
  
    autoreconf -fi
    
  If for some reason, the autoreconf command is not available,
  the following sequence of commands should have the same result:
    
    autoheader
    aclocal -I m4
    libtoolize -f
    autoconf 
    automake -a
    

Why aren't the configure script and other generated files in the
CVS repository?
  
  Because they are generated files.  Why save a version history for
  them?  Further, some of the above commands get re-run automatically
  when Makefile.am's or other files are changed.  This could lead to 
  compatibility problems if some of the generated files in the CVS 
  repository are from a different version of the GNU autotools.
  

Aren't we requiring users to have GNU autotools installed in order 
to configure MOAB?

  No.  Developers (or anyone else using source directly from the CVS
  repository) must have the autotools installed.  When creating a
  tarball for distribution of MOAB, the commands below should be run.
  The resulting tarball will contain all necessary generated files,
  including the configure script.


What needs to be done to prepare MOAB for distribution as a tarball?

  Check out a clean copy of MOAB.
  Execute the following commands in the top-most directory:
    autoreconf -fi
    ./configure
    make dist
To create a distributable tarball from a working source directory, do
  make dist
  
To create a distributable tarball from a clean check-out of the MOAB
source, do:
  cd MOAB
  aclocal
  libtoolize
  automake -a
  autoconf
  ./configure --without-netcdf --without-hdf5
  make dist

MOAB iMesh Interface Implementation, iMesh v1.2
===============================================

1. Non-compliant iMesh functionality:

1.1 Iterators for list-type entity sets: 

The iMesh 1.2 specification requires that iterators over list-type entity sets be updated
in response to membership changes in the set.  Specifically, if entities are added to or
removed from the set, the spec requires that added entities be iterated without needing to
reset the iterator, and that removed entities not be iterated.  MOAB will not support this
capability in the iMesh 1.2 release.  Future support will depend on whether this can be
implemented efficiently, without degrading performance for static mesh applications.

1.2 No support for septahedron entities:

MOAB does not support septahedron entities at this time (though such entities could be
represented as general polyhedra).

2. MOAB capabilities not accessible through iMesh:

2.1 Dense tags: MOAB supports two kinds of tag storage: dense tags, where tag values are
stored in sequence for sequences of contiguous entity handles; and sparse tags, which are
stored in (entity handle, tag value) tuples.  iMesh does not support creation of a tag
with a default value, nor does it have a mechanism for passing general options to the tag
creation function.  Therefore, MOAB's iMesh implementation creates sparse tags by
default.  Alternatives for specifying the tag creation type will be explored for future
iMesh releases.

2.2 Variable-length tags: MOAB supports a variable-length tag, where a tag value can have
a different length for each entity to which it is assigned.  This functionality is not
supported in iMesh.

2.3 Direct access to tag data (tag_iterate functionality): MOAB 4.0.1 introduced the
ability for applications to get direct access to tag storage for dense-type tags (see the
tag_iterate function in src/moab/Interface.hpp).  This functionality is not supported in
iMesh.

2.4 Corner vs. all vertices in connectivity list: MOAB represents so-called "higher-order
entities", e.g. quadratic tetrahedra, by allowing the connectivity list to be an
application-specified size.  The connectivity array returned by MOAB's iMesh
implementation will always be the total number of vertices, including any high-order
vertices.  MOAB's interface allows applications to specify whether all or just corner
vertices are requested.

2.5 Retrieval of entity, set handles in order from list-type sets: MOAB uses the same
handle type for both entities and entity sets.  The order of addition of entities and
entity sets to list-type sets is therefore preserved, since these handles are stored in
the same list.  Since iMesh uses a different handle type for entities than for sets, and
returns those handles in different API functions (iMesh_getEntities versus
iMesh_getEntSets), the original order of addition in the case where entities and entity
sets are interspersed cannot be recovered.

2.6 No support for knife-type elements: MOAB supports a seven-sided element referred to as
a "Knife" element; this element results from collapsing a quadrilateral bounding a
hexahedron by merging opposite nodes.  This element type is not supported in iMesh.

Other errata:

None

This document describes file formats supported by MOAB and options
controlling file input and output.


Supported File Formats
----------------------

Some of the file formats listed below may not be supported by a particular
build of MOAB depending on the availability of external libraries.  An 
up-to-date list of file formats supported in a particular build of MOAB
can be obtained programatically using the MBReaderWriterSet API or as 
a simple list using the '-l' option of the mbconvert utility.

Format              Name    Read         Write   File name suffixes
------------------  ------  ------------ ------  --------------------
MOAB native         MOAB    yes          yes     h5m mhdf
Exodus II           EXODUS  yes          yes     exo exoII exo2 g gen
Kitware VTK         VTK     up to v3.0   v3.0    vtk
Cubit               CUBIT   yes          no      cub
SLAC                SLAC    no           yes     slac
GMV                 GMV     no           yes     gmv
Ansys               ANSYS   no           yes     ans
Gmsh                GMSH    v1.0, v2.0   v2.0    msh gmsh
Stereo Lithography  STL     yes          yes     stl

Any of the values from the 'Name' column may be passed as an option
to the file IO methods to request a particular file.  If no file 
format is specified, the default is to choose the write format using
the file extension and to try all file readers until one succeeds.


---------------
File IO Options
---------------

An options list as passed to MOAB file IO routines is a single C-style
string containing the concatenation of a list of string options, were
individual options are separated by a designated separator character.
The default separator character is a semicolon (;).  To specify an alternate
separator character, begin the options string with a semicolon followed
by the desired separator.  Options are not case sensitive.


---------------
Common Options
---------------

  PRECISION=<N>

Specify the precision to use when writing float and double values
(such as node coordinates) to text-based file formats.


  CREATE
  
Do not overwrite existing file.

  FACET_DISTANCE_TOLERANCE=<D>

Max. distance deviation between facet and geometry, default:0.001.

  FACET_NORMAL_TOLERANCE=<N>

Max. normal angle deviation (degrees) between facet and geometry, default:5.

  MAX_FACET_EDGE_LENGTH=<D>

Max. facet edge length, default:0.

  CGM_ATTRIBS={yes|no}

Actuation of all CGM attributes, default:no.

  DEBUG_IO=n
  
Set threashold for debug messages from low-level (format-specific) 
reader/writer.  Default is 0 (none).

-------------------
Parallel IO Options
-------------------

MOAB must be built with parallel support enabled before these options
can be used.  

Parallel Read Options:
----------------------

  PARALLEL={NONE|BCAST|BCAST_DELETE|READ_DELETE|READ_PART}

Set parallel read mode.  Options are:
  - NONE  - force serial read/write on each processor (default)
  - BCAST - read on one processor and broadcast a copy to all others
  - BCAST_DELETE - read file on one processor, broadcasting partitioned
                   data to other processors.
  - READ_DELETE - read file on all processors and partition by deleting
                  mesh from non-local partitions.
  - READ_PART - read only relevant subset of file on each processor,
                utilizing format-specific parallel I/O if available
                (this option is not supported for all file formats.)
  - WRITE_PART - write only relevant subset of file on each processor,
                 creating a single file utilizing format-specific parallel 
                 I/O if available (this option is not supported for all 
                 file formats.)
  - FORMAT - depricated (use WRITE_PART)


  PARALLEL_RESOLVE_SHARED_ENTS

Resolve which entities are shared between which processes, such
that propogation of data across processes can be done.  This should
probably be the defautl behavior, as you almost certainly want this
unless, perhaps, PARALLEL=BCAST.

  
  PARTITION
  PARTITION=<tag_name>

Specify that mesh should be partitioned using partition IDs stored
in a tag.  If the tag name is not specified, the default ("PARTITION")
is used.


  PARTITION_VAL=<int_list>

If a tag name is specified to the 'PARTITION' option, then treat as
partitions only those sets for which the tag *value* is a single integer
for which the value is one of integers in thie specified list.


  PARTITION_DISTRIBUTE

Depricated.  Implied by "PARTITION" option.


  PARTITION_BY_RANK
  
Assume 1-1 mapping between MPI rank and part ID.  Assing parts
to processors for which rank == part id.


  MPI_IO_RANK=<RANK>

For IO modes in which a single processor handles all disk access, the
MPI rank of the processor to use.  Default is 0.

  PARALLEM_COMM=<id>

Specify moab::ParallelComm to use as communicator.

-----------------------
Format-specific options
-----------------------

Stereo Lithography (STL) files
------------------------------
  
   BINARY|ASCII
   
Write binary or text STL file.  Default is text.

   BIG_ENDIAN|LITTLE_ENDIAN
   
Force byte ordering of binary data.  Default is BIG_ENDIAN for writing
and autodetection for reading (BIG_ENDIAN if autodetect fails).


MOAB Native (HDF5-based MHDF) format
------------------------------------

   ELEMENTS={EXPLICIT|NODES|SIDES}

If reading only part of a file, specify which elements to read
in addition to those contained in the specified set.  The options are:
   EXPLICIT    - read only explicitly designated elements
   NODES       - read any element for which all the nodes are being read.
   SIDES       - read explicilty specified elements and any elements that
                 are sides of those elements.
Default is SIDES unless only nodes are explicitly specified, in which
case NODES will be used.

   CHILDREN={NONE|SETS|CONTENTS}
   SETS={NONE|SETS|CONTETS}

If reading only part of a file, specify whether or not child or
contained sets (CHILDREN an SETS, respectively) of input sets are to be read.  
The options are:
   NONE     - Do not read sets because they are children of designated sets.
   SETS     - Read all child sets of designated input sets.
   CONTENTS - (Default).  Read all child sets and any entities contained
              in those sets.


   BUFFER_SIZE=<BYTES>
   
Reduce buffer size for debugging purposes.

   KEEP

For debugging purposes, retain partially written file if a failure
occurs during the write process.

   BLOCKED_COORDINATE_IO={yes|no}
   
During read of HDF5 file, read vertex coordinates in blocked format 
(read all X coordinates, followed by all Y coordinates, etc.)  Default
is 'no'.

   BCAST_SUMMARY={yes|no}
   
During parallel read of HDF5 file, read file summary data on root process
as serial IO and broadcast summary to other processes.  All processes then
re-open file for parallel IO.  If 'no', file is opened only once by all
processes for parallel IO and all processes read summary data.  Default is
'yes'.

   BCAST_DUPLICATE_READS={yes|no}

Work around MPI IO limitations when all processes read an indentical region of
a file by reading the data only on the root process and using MPI_Bcast to
send the data to other processes.  Default is 'yes'.

   HYPERSLAB_OR 

During partial or parallel read, use documented (default) method for
building HDF5 hyperslabs.  This option is deprecated.  Reader should
automatically detect if HYPERSLAB_APPEND works and if not will default
to HYPERSLAB_OR.  This option is retaind for debugging purposes only.

   HYPERSLAB_APPEND
  
During partial or parallel read using a modified (i.e. hacked) HDF5
library, utilize hack for more efficient hyperslab selection construction.
This option is deprecated.  Reader should automatically detect if 
HYPERSLAB_APPEND works and if not will default to HYPERSLAB_OR.

   HYPERSLAB_SELECT_LIMIT=n
   
Set upper bound on the number of HDF5 hyperslabs that can be combined for
a single read from an dataset.  The defualt is defined by 
DEFUALT_HYPERSLAB_SELECT_LIMIT in ReadHDF5Dataset.cpp (200 at the time
this was written).  If HYPERSLAB_APPEND is specified and this option is
not, then the default is no limit.  This limit should be removed for 
future HDF5 releases (> 1.8.x) as said future releases will no longer 
require the HYPERSLAB_APPEND hack in order to avoid O(n^2) hyperslab 
selection behavior.

