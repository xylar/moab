# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include makefile.config

default: all

EXAMPLES = HelloMOAB GetEntities SetsNTags LoadPartial StructuredMesh StructuredMeshSimple DirectAccessWithHoles DirectAccessNoHoles PointInElementSearch GenLargeMesh UniformRefinement ExtrudePoly
PAREXAMPLES = HelloParMOAB ReduceExchangeTags LloydRelaxation CrystalRouterExample ReadWriteTest
EXOIIEXAMPLES = TestExodusII VisTags
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
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloMOAB: HelloMOAB
	${RUNSERIAL} ./HelloMOAB

ExtrudePoly: ExtrudePoly.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ExtrudePoly: ExtrudePoly
	${RUNSERIAL} ./ExtrudePoly

GetEntities: GetEntities.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-GetEntities: GetEntities
	${RUNSERIAL} ./GetEntities

SetsNTags: SetsNTags.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-SetsNTags: SetsNTags
	${RUNSERIAL} ./SetsNTags

LloydRelaxation: LloydRelaxation.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-LloydRelaxation: LloydRelaxation
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./LloydRelaxation
	${RUNPARALLEL} ./LloydRelaxation
endif

LoadPartial: LoadPartial.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} ${MOAB_LIBS_LINK}

run-LoadPartial: LoadPartial
	${RUNSERIAL} ./LoadPartial

StructuredMesh: StructuredMesh.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-StructuredMesh: StructuredMesh
	${RUNSERIAL} ./StructuredMesh

StructuredMeshSimple: StructuredMeshSimple.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-StructuredMeshSimple: StructuredMeshSimple
	${RUNSERIAL} ./StructuredMeshSimple -d 3 -n 5

DirectAccessWithHoles: DirectAccessWithHoles.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-DirectAccessWithHoles: DirectAccessWithHoles
	${RUNSERIAL} ./DirectAccessWithHoles -n 1000 -H 5

DirectAccessNoHoles: DirectAccessNoHoles.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-DirectAccessNoHoles: DirectAccessNoHoles
	${RUNSERIAL} ./DirectAccessNoHoles -n 1000

ReduceExchangeTags: ReduceExchangeTags.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReduceExchangeTags: ReduceExchangeTags
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./ReduceExchangeTags
	${RUNPARALLEL} ./ReduceExchangeTags
endif

HelloParMOAB: HelloParMOAB.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloParMOAB: HelloParMOAB
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./HelloParMOAB
endif

CrystalRouterExample: CrystalRouterExample.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-CrystalRouterExample: CrystalRouterExample
	${RUNSERIAL} ./CrystalRouterExample
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./CrystalRouterExample -n 2 -t 10 -r 0
endif

TestExodusII: TestExodusII.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-TestExodusII: TestExodusII
	${RUNSERIAL} ./TestExodusII

PointInElementSearch: PointInElementSearch.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-PointInElementSearch: PointInElementSearch
	${RUNSERIAL} ./PointInElementSearch
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./PointInElementSearch
endif

VisTags: VisTags.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-VisTags: VisTags
	${RUNSERIAL} ./VisTags

ReadWriteTest: ReadWriteTest.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReadWriteTest: ReadWriteTest
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)-$(MOAB_NETCDF_ENABLED)","yes-yes-yes")
	${RUNPARALLEL} ./ReadWriteTest
endif

GenLargeMesh: GenLargeMesh.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-GenLargeMesh: GenLargeMesh
	${RUNSERIAL} ./GenLargeMesh
	${RUNSERIAL} ./GenLargeMesh -A 2 -B 2 -C 5
	${RUNSERIAL} ./GenLargeMesh -A 2 -B 2 -C 5 -f
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./GenLargeMesh -M @NP@ -N 1 -K 1 -A 1 -B 2 -C 2 -r
	${RUNPARALLEL} ./GenLargeMesh -M @NP@ -N 1 -K 1 -A 1 -B 2 -C 2 -d 10.0 -b 2 -r
	${RUNPARALLEL} ./GenLargeMesh -M @NP@ -N 1 -K 1 -A 1 -B 2 -C 2 -d 10.0 -b 2 -r -f
	${RUNPARALLEL} ./GenLargeMesh -M @NP@  -b 10  -A 1 -C 1 -B 1 -r -G
endif

ErrorHandlingSimulation: ErrorHandlingSimulation.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandling: TestErrorHandling.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandlingPar: TestErrorHandlingPar.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

UniformRefinement: UniformRefinement.o ${MOAB_LIBDIR}/libMOAB.la
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

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

