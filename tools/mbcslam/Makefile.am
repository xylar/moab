MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/tools/mbcoupler \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/test \
               -I$(top_srcdir)/src/IntxMesh

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
LDADD =  $(top_builddir)/src/libMOAB.la

TESTS = case1_test
if PARALLEL
if HAVE_HDF5
   TESTS += linear_remap
endif
endif

noinst_PROGRAMS =  proj1 create_dp
if PARALLEL
endif

check_PROGRAMS = $(TESTS)
DEFS += -DBUILDING_MB \
        -DSRCDIR=$(srcdir) \
        -DMESHDIR=$(MESHDIR)

case1_test_SOURCES = case1_test.cpp
proj1_SOURCES = proj1.cpp
create_dp_SOURCES = create_dp.cpp
# this test is for netcdf only ; maybe convert to vtk for the test
linear_remap_SOURCES = linear_remap.cpp

if PARALLEL
  linear_advection_SOURCES = linear_advection.cpp
  noinst_PROGRAMS += linear_advection
if HAVE_HDF5
  diffusion_SOURCES = diffusion.cpp
  intx_imesh_SOURCES = intx_imesh.cpp update_tracer.cpp
  intx_imesh_LDADD =  $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la
  cslam_par_test_SOURCES = cslam_par_test.cpp
  intx_mpas_SOURCES = intx_mpas.cpp
  noinst_PROGRAMS += cslam_par_test intx_imesh intx_mpas diffusion
if ENABLE_FORTRAN
  advection_SOURCES = advection.F90
  advection_LDADD =  update_tracer.o $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la $(FCLIBS)
  noinst_PROGRAMS += advection
endif
endif
endif

FORTDEFS = $(FCDEFS) #$(DEFS:-D=$(FC_DEFINE))
PPFCCOMPILE = $(FC) $(FORTDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
    $(AM_CPPFLAGS) $(FCPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
    $(LIBTOOLFLAGS) --mode=compile $(FC) $(FORTDEFS) \
    $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(FCPPFLAGS) \
    $(AM_FCFLAGS) $(FCFLAGS)

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = intersect1.h5m \
                   lagr.h5m \
                   polyWithEdges.vtk \
                   lagr00.h5m \
                   intersect00.h5m \
                   lagrIni.h5m \
                   intx.vtk \
                   spectral.vtk \
                   intx1.vtk \
                   SpanEdges0.h5m \
                   SpanQuads0.h5m
