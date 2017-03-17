# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libRemapping.la
AM_CPPFLAGS =

# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB -DSRCDIR=$(srcdir) -DLOCDIR=src/Remapping/ \
               -I$(srcdir)/..  -I.. -I$(srcdir)/../moab -I$(srcdir)/moab -I$(top_srcdir)/src/IntxMesh \
	       -I$(srcdir)/../parallel -I$(top_srcdir)/itaps -I$(top_srcdir)/itaps/imesh

# The directory in which to install headers
libRemapping_la_includedir = $(includedir)

# Conditional sources
MOAB_REMAPPING_SRCS = \
     Remapper.cpp  \
     RemapperUtils.cpp

MOAB_REMAPPING_HDRS = \
	moab/Remapper.hpp \
  moab/RemapperUtils.hpp

if ENABLE_mbtempest
  AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
  AM_LDFLAGS += $(TEMPESTREMAP_LIBS)

  MOAB_REMAPPING_SRCS += \
     TempestRemapper.cpp  \
     TempestLinearRemap.cpp \
     TempestOfflineMap.cpp

  MOAB_REMAPPING_HDRS += \
     moab/TempestRemapper.hpp \
     moab/TempestOfflineMap.hpp

endif

# The list of source files, and any header files that do not need to be installed
libRemapping_la_SOURCES = \
   $(MOAB_REMAPPING_SRCS)

# The list of header files which are to be installed
nobase_libRemapping_la_include_HEADERS = \
  $(MOAB_REMAPPING_HDRS)
