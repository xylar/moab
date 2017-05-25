# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../../makefile.config

default: all

ifeq ("$(MOAB_HDF5_ENABLED)","yes")
ALLEXAMPLES += ContinentsOnGlobe BoundaryContinents
endif

all: $(ALLEXAMPLES)

ContinentsOnGlobe: ContinentsOnGlobe.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ContinentsOnGlobe: ContinentsOnGlobe
	${RUNSERIAL} ./ContinentsOnGlobe

BoundaryContinents: BoundaryContinents.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-BoundaryContinents: BoundaryContinents
	${RUNSERIAL} ./BoundaryContinents

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

