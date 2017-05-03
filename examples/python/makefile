# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../makefile.config

default: all

ALLEXAMPLES = yt2moab laplaciansmoother

all: $(ALLEXAMPLES)

yt2moab: yt2moab.py
	@echo "[Python]  $@"
	${PYTHON} yt2moab.py

laplaciansmoother: laplaciansmoother.py
	@echo "[Python]  $@"
	${PYTHON} laplaciansmoother.py ${MESH_DIR}/surfrandomtris-4part.h5m 25

run: all 

clean: clobber
	rm -rf ${ALLEXAMPLES}

