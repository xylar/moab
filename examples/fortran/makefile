# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

F90EXAMPLES = DirectAccessNoHolesF90
F90PAREXAMPLES = PushParMeshIntoMoabF90 MigrateMesh

ALLEXAMPLES =
ifeq ("$(MOAB_FORTRAN_ENABLED)-$(MOAB_IMESH_ENABLED)","yes-yes")
ALLEXAMPLES += ${F90EXAMPLES} 
endif

ifeq ("$(MOAB_MPI_ENABLED)","yes")
ifeq ("$(MOAB_FORTRAN_ENABLED)-$(MOAB_IMESH_ENABLED)","yes-yes")
parallel: ${F90PAREXAMPLES}
ALLEXAMPLES += ${F90PAREXAMPLES}
endif
else
parallel:
endif

all: $(ALLEXAMPLES)

DirectAccessNoHolesF90: DirectAccessNoHolesF90.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${IMESH_LIBS}

run-DirectAccessNoHolesF90: DirectAccessNoHolesF90
	${RUNSERIAL} ./DirectAccessNoHolesF90 -n 1000

PushParMeshIntoMoabF90: PushParMeshIntoMoabF90.o
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${IMESH_LIBS}

run-PushParMeshIntoMoabF90: PushParMeshIntoMoabF90
	${RUNSERIAL} ./PushParMeshIntoMoabF90
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./PushParMeshIntoMoabF90
endif

MigrateMesh: MigrateMesh.o
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${IMESH_LIBS}

run-MigrateMesh: MigrateMesh
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./MigrateMesh
endif

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

