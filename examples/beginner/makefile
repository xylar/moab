# Include the MOAB configuration information so that  
# all required flags and libs are populated correctly
include ../makefile.config

default: all

ALLEXAMPLES = mbex1 mbex2 mbex3 mbex4

all: $(ALLEXAMPLES)

mbex1: mbex1.o
mbex2: mbex2.o
mbex3: mbex3.o
mbex4: mbex4.o

$(ALLEXAMPLES):
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-mbex1: mbex1
#The first example must be run before the second.
run-mbex2: mbex2 run-mbex1
run-mbex3: mbex3
run-mbex4: mbex4

$(addprefix run-,$(ALLEXAMPLES)):
	./$<

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}
