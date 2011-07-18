#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'
ifndef VERBOSE
  QUIET := @
endif

CC = g++
CMSROOT = ./
ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

CORESOURCES=./CORE/CMS2.cc \
./CORE/conversionTools.cc \
./CORE/electronSelections.cc \
./CORE/electronSelectionsParameters.cc \
./CORE/eventSelections.cc \
./CORE/jetSelections.cc \
./CORE/mcSelections.cc \
./CORE/metSelections.cc \
./CORE/MITConversionUtilities.cc \
./CORE/muonSelections.cc \
./CORE/trackSelections.cc \
./CORE/triggerUtils.cc \
./CORE/utilities.cc 
COREOBJECTS=$(CORESOURCES:.cc=.o)
CORELIB=libCMS2NtupleMacrosCORE.so

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = liblooper.so

LIBS = $(CORELIB) $(LIB)

#
# how to make it
#

libs:	$(LIBS)

$(CORELIB): $(COREOBJECTS)
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(INCLUDE) LinkDef.h; \
	cat LinkDef.h LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx looper.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

#
# target to build
#
build: $(LIBS) # compile code
	@#root -b -q ".L MyScanChain.C++"
b: build

clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.cxx \
	rm -f *.so \
	rm -f LinkDef_out* \
	rm -f ./CORE/*.o \
	rm -f ./CORE/*.d \
	rm -f ./CORE/*.cxx \
	rm -f ./CORE/*.so

endif
