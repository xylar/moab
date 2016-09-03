# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libIntxMesh.la         
AM_CPPFLAGS =

# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB -DSRCDIR=$(srcdir) -DLOCDIR=src/IntxMesh/ \
               -I$(srcdir)/..  -I.. -I$(srcdir)/../moab -I$(srcdir)/moab \
	       -I$(srcdir)/../parallel -I$(top_srcdir)/itaps -I$(top_srcdir)/itaps/imesh

# The directory in which to install headers
libIntxMesh_la_includedir = $(includedir)

# Conditional sources
MOAB_LOCALDISCR_SRCS = \
     Intx2Mesh.cpp  \
     Intx2MeshInPlane.cpp  \
     Intx2MeshOnSphere.cpp \
     IntxRllCssphere.cpp \
     IntxUtils.cpp \
     IntxUtils.hpp 

if PARALLEL
    MOAB_LOCALDISCR_SRCS += wrap_intx.cpp
endif

MOAB_LOCALDISCR_HDRS = \
	moab/Intx2Mesh.hpp \
        moab/Intx2MeshInPlane.hpp \
        moab/Intx2MeshOnSphere.hpp \
        moab/IntxRllCssphere.hpp   

# The list of source files, and any header files that do not need to be installed
libIntxMesh_la_SOURCES = \
   $(MOAB_LOCALDISCR_SRCS)

# The list of header files which are to be installed
nobase_libIntxMesh_la_include_HEADERS = \
  $(MOAB_LOCALDISCR_HDRS)
