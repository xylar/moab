# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../makefile.config

default: all

ALLEXAMPLES = yt2moab laplaciansmoother

all: $(ALLEXAMPLES)

yt2moab: yt2moab.py
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

laplaciansmoother: laplaciansmoother.py
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

