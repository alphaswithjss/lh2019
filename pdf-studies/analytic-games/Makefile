
# Makefile generated automatically by /home/soyez/local/bin/mkcxx.pl '-s' '-c' '-l' '-lpython2.7'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

INCLUDE += -I../../helpers
LIBRARIES  += -L../../helpers -lCmdLine
LIBRARIES += -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

INCLUDE += $(LCLINCLUDE)

COMMONSRC = ecfs.cc hypergeometric.cc qcd-constants.cc resum-blocks-vincent.cc resum-blocks.cc
F77SRC = 
COMMONOBJ = ecfs.o hypergeometric.o qcd-constants.o resum-blocks-vincent.o resum-blocks.o

PROGSRC = test-ecfs.cc test-lh-multiplicity.cc
PROGOBJ = test-ecfs.o test-lh-multiplicity.o

INCLUDE += 
LIBRARIES += -lpython2.7


all:  test-ecfs test-lh-multiplicity 


test-ecfs: test-ecfs.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

test-lh-multiplicity: test-lh-multiplicity.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	/home/soyez/local/bin/mkcxx.pl '-s' '-c' '-l' '-lpython2.7'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  test-ecfs test-lh-multiplicity 

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

ecfs.o: qcd-constants.hh ecfs.hh resum-blocks-vincent.hh hypergeometric.hh
hypergeometric.o: hypergeometric.hh
qcd-constants.o: qcd-constants.hh
resum-blocks-vincent.o: qcd-constants.hh resum-blocks-vincent.hh
resum-blocks.o: resum-blocks.hh qcd-constants.hh
test-ecfs.o: ecfs.hh
test-lh-multiplicity.o: lh-multiplicity.hh resum-blocks.hh qcd-constants.hh
