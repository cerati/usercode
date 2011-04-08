include Makefile.arch

#
# stuff to make
#

CORESOURCES=../CORE/CMS2.cc \
	../CORE/eventSelections.cc \
	../CORE/muonSelections.cc \
	../CORE/electronSelections.cc \
	../CORE/electronSelectionsParameters.cc \
	../CORE/MITConversionUtilities.cc \
	../CORE/trackSelections.cc 
COREOBJECTS=$(CORESOURCES:.cc=.o)
CORELIB=libCMS2NtupleMacrosCORE.so

#
# how to make it
#

$(CORELIB): $(COREOBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) $(COREOBJECTS) -o $@


# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.o: 	%.cc 
	$(QUIET) echo "Compiling $<"; \
	$(CXX) -shared $(CXXFLAGS) $(ROOTLIBS) -I$(INCLUDE) -c $< -o $@

%.o: 	%.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

#
# target to build
#

clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.cxx \
	rm -f *.so
	rm -f ../CORE/*.o \
	rm -f ../CORE/*.d \
	rm -f ../CORE/*.cxx \
	rm -f ../CORE/*.so
