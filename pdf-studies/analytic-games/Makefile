
# Makefile generated automatically by /home/soyez/local/bin/mkcxx.pl '-s' '-c'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

INCLUDE += -I../../helpers
LIBRARIES  += -L../../helpers -lCmdLine
LIBRARIES += -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

INCLUDE += $(LCLINCLUDE)

COMMONSRC = qcd-constants.cc resum-blocks.cc
F77SRC = 
COMMONOBJ = qcd-constants.o resum-blocks.o

PROGSRC = test-lh-multiplicity.cc
PROGOBJ = test-lh-multiplicity.o

INCLUDE += 
LIBRARIES += 


all:  test-lh-multiplicity 


test-lh-multiplicity: test-lh-multiplicity.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	/home/soyez/local/bin/mkcxx.pl '-s' '-c'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  test-lh-multiplicity 

.cc.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.C.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.f.o:         $<
	$(F77) $(FFLAGS) -c $< -o $@


depend:
	makedepend  $(LCLINCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE

qcd-constants.o: qcd-constants.hh
resum-blocks.o: resum-blocks.hh qcd-constants.hh
