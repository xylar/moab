# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
MOAB_BUILD_DIR := /homes/iulian/source/MOAB14
MOAB_DIR := /usr/software/moab/dev/feature_check

ifneq ($(wildcard ${MOAB_DIR}/lib/moab.make),)
include ${MOAB_DIR}/lib/moab.make
include ${MOAB_DIR}/lib/iMesh-Defs.inc
else
include ${MOAB_BUILD_DIR}/moab.make
include ${MOAB_BUILD_DIR}/itaps/imesh/iMesh-Defs.inc
endif

.SUFFIXES: .o .cpp

all: libiMOAB

libiMOAB: iMOAB.o
	${AR} rcs libiMOAB.a $< 

clean:
	rm -rf *.o *.h5m 

.cpp.o:
	${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

