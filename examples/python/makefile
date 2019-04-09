# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../makefile.config

default: all

ALLEXAMPLES = yt2moab laplaciansmoother H5MToScrip ScripToH5M

all: $(ALLEXAMPLES)

laplaciansmoother: laplaciansmoother.py
	@echo "[Python]  $@"
	${PYTHON} laplaciansmoother.py ${MESH_DIR}/surfrandomtris-4part.h5m 25

H5MToScrip: H5MToScrip.py
	@echo "[Python]  $@"
	${PYTHON} H5MToScrip.py -m ${MESH_DIR}/surfrandomtris-4part.h5m -s surfrandomtris-4part.nc

ScripToH5M: ScripToH5M.py
	@echo "[Python]  $@"
	${PYTHON} ScripToH5M.py -m ${MESH_DIR}/surfrandomtris-4part-converted.h5m -s surfrandomtris-4part.nc

run: all 

clean: clobber
	rm -rf ${ALLEXAMPLES}

