# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I config

# Subdirectories to build
SUBDIRS = src MeshFiles $(itaps_dir) tools test 

if ENABLE_fbigeom
  itaps_dir_fbigeom = itaps
else
  itaps_dir_fbigeom =
endif

if ENABLE_imesh
  itaps_dir = itaps
else
  itaps_dir = $(itaps_dir_fbigeom)
endif

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

doc_DATA = ANNOUNCE KNOWN_ISSUES LICENSE README.md RELEASE_NOTES

examples_DATA = examples/CrystalRouterExample.cpp \
                examples/DirectAccessNoHoles.cpp \
                examples/DirectAccessNoHolesF90.F90 \
                examples/DirectAccessWithHoles.cpp \
                examples/ErrorHandlingSimulation.cpp \
                examples/GenLargeMesh.cpp \
                examples/GetEntities.cpp \
                examples/HelloMOAB.cpp \
                examples/HelloParMOAB.cpp \
                examples/LloydRelaxation.cpp \
                examples/LoadPartial.cpp \
                examples/PointInElementSearch.cpp \
                examples/PushParMeshIntoMoabF90.F90 \
                examples/ReadWriteTest.cpp \
                examples/ReduceExchangeTags.cpp \
                examples/StructuredMesh.cpp \
                examples/StructuredMeshSimple.cpp \
                examples/SetsNTags.cpp \
                examples/TestErrorHandling.cpp \
                examples/TestErrorHandlingPar.cpp \
                examples/TestExodusII.cpp \
                examples/UniformRefinement.cpp \
                examples/VisTags.cpp \
                examples/ExtrudePoly.cpp \
                examples/makefile.in \
                examples/makefile

CMAKE_CONFIGURATION_FILES = CMakeLists.txt \
														config/AutoconfHeader.cmake      config/FindCGM.cmake        config/FindMetis.cmake     config/FindZoltan.cmake    config/ForceAddFlags.cmake \
														config/ResolveCompilerPaths.cmake  config/CMakeDistCheck.cmake      config/FindEigen3.cmake     config/FindNetCDF.cmake    config/Findgraphviz.cmake \
														config/GetAcInitVersion.cmake  config/TestBigEndian.cmake    config/CheckCompilerFlags.cmake  config/FindHDF5_MOAB.cmake  config/FindParMetis.cmake  \
														config/Findverdict.cmake   config/MOABConfig.cmake.cmake    config/CMakeReplicateConfig.sh    config/cmake_uninstall.cmake.in \
														config/TemplateSpecialization.cxx   config/MOABConfig.cmake.h.in

EXTRA_DIST = $(doc_DATA) $(examples_DATA) $(CMAKE_CONFIGURATION_FILES)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = moab.make moab.config
cmakedir = $(libdir)/cmake/MOAB/
cmake_DATA = MOABConfig.cmake

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	$(AM_V_at)echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/moab.make
	$(AM_V_at)echo "MOAB_INCLUDES=-I${includedir}" >> $(DESTDIR)$(cfgdir)/moab.make

# Generate a file to be installed in $libdir containing the configuration
# options used for this MOAB build.
# Note: If you modify this, verify that it works correctly for option
#       values containing spaces.
# First line: get unformatted data from config.status
# Second line: extract just the list of options
# Third line: put individual options on separate lines
# Fourth line: remove any empty lines and trim spaces and then write to config file.
moab.config: config.status Makefile.am
	$(AM_V_GEN)./config.status --version | \
	  sed -e 's/.*options "\(.*\)"/\1/p' -e 'd' | \
	  tr "'" "\n" | \
	  sed -e '/^ *$$/d' > $@

CLEANFILES = moab.config a.out configs.sed

