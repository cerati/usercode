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

LIBS = $(CORELIB)

DICT = LinkDef_out.cc

EXE = main.exe

#
# how to make it
#

$(CORELIB): $(COREOBJECTS)
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

$(DICT): LinkDef.h
	rootcint -f $@ -c -p $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h LinkDef.h

%.exe: $(LIBS) $(DICT) looper.cc
	$(CC) -o $@ `root-config --ldflags --cflags --libs` -lEG -lGenVector -L. -lCMS2NtupleMacrosCORE LinkDef_out.cc looper.cc ${@:.exe=.cc} 

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

clean:
	rm -f *.o \
	rm -f *.cxx \
	rm -f *_out.* \
	rm -f *.*~ \
	rm -f *.so \
	rm -f ./CORE/*.o \
	rm -f ./CORE/*.so \
	rm main.exe

endif
