# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
MOAB_BUILD_DIR := /homes/iulian/source/MOAB14
MOAB_DIR := /homes/iulian/lib/moab14

ifneq ($(wildcard ${MOAB_DIR}/lib/moab.make),)
include ${MOAB_DIR}/lib/moab.make
include ${MOAB_DIR}/lib/iMesh-Defs.inc
else
include ${MOAB_BUILD_DIR}/moab.make
include ${MOAB_BUILD_DIR}/itaps/imesh/iMesh-Defs.inc
endif

.SUFFIXES: .o .cpp .F90


driver: driver.o ../src/imoab.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ../src/imoab.o ${MOAB_LIBS_LINK}

clean:
	rm -rf *.o *.mod *.h5m ${ALLEXAMPLES}

.cpp.o:
	${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.c.o:
	${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o:
	${IMESH_FC} ${FCFLAGS} ${IMESH_FCFLAGS} ${MOAB_CPPFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

