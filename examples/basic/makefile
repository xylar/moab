# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES = HelloMOAB GetEntities SetsNTags StructuredMeshSimple UniformRefinement
PAREXAMPLES = HelloParMOAB ReadWriteTest
EXOIIEXAMPLES = TestExodusII
ifeq ("$(MOAB_DEV)","yes")
ERROREXAMPLES = ErrorHandlingSimulation TestErrorHandling TestErrorHandlingPar
endif
ALLEXAMPLES = ${EXAMPLES}

ifeq ("$(MOAB_NETCDF_ENABLED)","yes")
ALLEXAMPLES += $(EXOIIEXAMPLES)
endif

ifeq ("$(MOAB_MPI_ENABLED)","yes")
parallel: ${PAREXAMPLES}
ALLEXAMPLES += ${PAREXAMPLES}
else
parallel:
endif

all: $(ALLEXAMPLES)

HelloMOAB: HelloMOAB.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloMOAB: HelloMOAB
	${RUNSERIAL} ./HelloMOAB

GetEntities: GetEntities.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-GetEntities: GetEntities
	${RUNSERIAL} ./GetEntities

SetsNTags: SetsNTags.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-SetsNTags: SetsNTags
	${RUNSERIAL} ./SetsNTags

StructuredMeshSimple: StructuredMeshSimple.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-StructuredMeshSimple: StructuredMeshSimple
	${RUNSERIAL} ./StructuredMeshSimple -d 3 -n 5

HelloParMOAB: HelloParMOAB.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloParMOAB: HelloParMOAB
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./HelloParMOAB
endif

TestExodusII: TestExodusII.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-TestExodusII: TestExodusII
	${RUNSERIAL} ./TestExodusII

ReadWriteTest: ReadWriteTest.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReadWriteTest: ReadWriteTest
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)-$(MOAB_NETCDF_ENABLED)","yes-yes-yes")
	${RUNPARALLEL} ./ReadWriteTest
endif

ErrorHandlingSimulation: ErrorHandlingSimulation.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandling: TestErrorHandling.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandlingPar: TestErrorHandlingPar.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

UniformRefinement: UniformRefinement.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-UniformRefinement: UniformRefinement GenLargeMesh
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
ifeq ("$(wildcard GenLargeMesh.h5m)","")
	${RUNPARALLEL} ./GenLargeMesh -M @NP@ -N 1 -K 1 -A 1 -B 2 -C 2 -b 5
endif
	${RUNSERIAL}   ./UniformRefinement GenLargeMesh.h5m
	${RUNPARALLEL} ./UniformRefinement GenLargeMesh.h5m
else
ifeq ("$(wildcard GenLargeMesh.vtk)","")
	${RUNSERIAL} ./GenLargeMesh -A 1 -B 2 -C 2 -b 5
endif
	${RUNSERIAL} ./UniformRefinement GenLargeMesh.vtk
endif

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

