AM_CPPFLAGS = -I$(top_srcdir)/src \
              -I$(top_builddir)/src \
	      -I$(top_srcdir)/tools

AM_CXXFLAGS =
if HAVE_CGM
  AM_CPPFLAGS += $(CGM_INCLUDES) $(CGM_CPPFLAGS)
  AM_CXXFLAGS += $(CGM_CXXFLAGS)

if HAVE_NETCDF
  netcdf_progs = cub2h5m
else
  netcdf_progs =
endif

endif

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

bin_PROGRAMS = cgm2moab cub2h5m
noinst_PROGRAMS = quads_to_tris $(netcdf_progs)

CGM_CPPFLAGS = @CGM_CPPFLAGS@
CGM_LDFLAGS = @CGM_LDFLAGS@
CGM_LTFLAGS = @CGM_LTFLAGS@
CGM_LIBS = @CGM_LIBS@

LDADD =  $(top_builddir)/src/libMOAB.la
if HAVE_CGM
if USING_DEV_CGM
 LDADD += $(CGM_DIR)/src/libcgm.la
else
 LDADD += $(CGM_DIR)/lib/libcgm.la
endif
endif

cgm2moab_SOURCES = cgm2moab.cpp obb_analysis.cpp cgm2moab.hpp
cgm2moab_LDADD = $(top_builddir)/src/libMOAB.la

quads_to_tris_SOURCES = quads_to_tris_driver.cpp quads_to_tris.cpp
quads_to_tris_LDADD = $(top_builddir)/src/libMOAB.la

cub2h5m_SOURCES = cub2h5m.cpp quads_to_tris.cpp quads_to_tris.hpp
cub2h5m_LDADD = $(top_builddir)/src/libMOAB.la

EXTRA_DIST = CMakeLists.txt

