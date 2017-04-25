# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES = LoadPartial DirectAccessWithHoles DirectAccessNoHoles PointInElementSearch
PAREXAMPLES = ReduceExchangeTags
ALLEXAMPLES = ${EXAMPLES}

ifeq ("$(MOAB_MPI_ENABLED)","yes")
parallel: ${PAREXAMPLES}
ALLEXAMPLES += ${PAREXAMPLES}
else
parallel:
endif

all: $(ALLEXAMPLES)

LoadPartial: LoadPartial.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} ${MOAB_LIBS_LINK}

run-LoadPartial: LoadPartial
	${RUNSERIAL} ./LoadPartial

DirectAccessWithHoles: DirectAccessWithHoles.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-DirectAccessWithHoles: DirectAccessWithHoles
	${RUNSERIAL} ./DirectAccessWithHoles -n 1000 -H 5

DirectAccessNoHoles: DirectAccessNoHoles.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-DirectAccessNoHoles: DirectAccessNoHoles
	${RUNSERIAL} ./DirectAccessNoHoles -n 1000

ReduceExchangeTags: ReduceExchangeTags.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReduceExchangeTags: ReduceExchangeTags
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./ReduceExchangeTags
	${RUNPARALLEL} ./ReduceExchangeTags
endif

PointInElementSearch: PointInElementSearch.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-PointInElementSearch: PointInElementSearch
	${RUNSERIAL} ./PointInElementSearch
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./PointInElementSearch
endif

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

