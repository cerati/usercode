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

#ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lEG -lGenVector

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h
#include "Math/QuantFuncMathCore.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

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
OBJECTS = $(SOURCES:.cc=.o)
LIB = liblooper.so

DICT = LinkDef_out.o

LIBS = $(CORELIB) $(LIB)

EXE = main.exe

#
# how to make it
#

$(CORELIB): $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

$(LIB):	$(DICT) $(OBJECTS) 
	$(QUIET) echo "Linking $@"; \
	$(LINKER) $(CFLAGS) $(LINKERFLAGS) -shared $(OBJECTS) $(DICT) -o $@

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE)  LinkDef.h; \
	cat LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

%.exe:  $(LIBS) 
	$(QUIET) echo "Building $@"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) ${@:.exe=.cc} 

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

loopclean:
	rm -f \
	*_out.*	 \
	*.o \
	*.*~ \
	$(LIB) \
	main.exe

clean: loopclean
	rm -f \
	$(CORELIB) \
	./CORE/*.o

endif
