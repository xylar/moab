# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../makefile.config

default: all

ALLEXAMPLES = mbex1 mbex2 mbex3 mbex4

all: $(ALLEXAMPLES)

mbex1: mbex1.o
	$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

mbex2: mbex2.o
	$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

mbex3: mbex3.o
	$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

mbex4: mbex4.o
	$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

