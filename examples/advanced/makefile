# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES =  ExtrudePoly test_trivial_distribution
PAREXAMPLES = LloydRelaxation CrystalRouterExample GenLargeMesh CoupleMGen
EXOIIEXAMPLES = VisTags
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

ExtrudePoly: ExtrudePoly.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ExtrudePoly: ExtrudePoly
	${RUNSERIAL} ./ExtrudePoly

LloydRelaxation: LloydRelaxation.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

test_trivial_distribution: test_trivial_distribution.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-LloydRelaxation: LloydRelaxation
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./LloydRelaxation
	${RUNPARALLEL} ./LloydRelaxation
endif

CrystalRouterExample: CrystalRouterExample.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-CrystalRouterExample: CrystalRouterExample
	${RUNSERIAL} ./CrystalRouterExample
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./CrystalRouterExample -n 2 -t 10 -r 0
endif

VisTags: VisTags.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-VisTags: VisTags
	${RUNSERIAL} ./VisTags

GenLargeMesh: GenLargeMesh.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

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

CoupleMGen: CoupleMGen.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< -lmbcoupler ${MOAB_LIBS_LINK}

run-CoupleMGen: CoupleMGen
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./CoupleMGen -M @NP@
endif


run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

