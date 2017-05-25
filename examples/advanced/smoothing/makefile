# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../../makefile.config

default: all

ALLEXAMPLES = LaplacianSmoother OptimizeMeshMesquite

all: $(ALLEXAMPLES)

LaplacianSmoother: LaplacianSmoother.o
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-LaplacianSmoother: LaplacianSmoother
	${RUNSERIAL} ./LaplacianSmoother

OptimizeMeshMesquite: OptimizeMeshMesquite.o
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-OptimizeMeshMesquite: OptimizeMeshMesquite
	${RUNSERIAL} ./OptimizeMeshMesquite

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

