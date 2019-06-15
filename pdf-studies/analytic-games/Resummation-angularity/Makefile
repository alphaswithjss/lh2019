# 05.08.2014, Vincent Theeuwes

FC = g++
LDLIBS=-lgsl -lgslcblas -lm
FFLAGS=-Wall -W
CXXFLAGS=-Wall -W -g -w
LIBS = -I/usr/include/python2.7 -lpython2.7 

SD-Angularity: main.o Integrand.o ResumFunc.o HypGeo.o
	$(FC) $^ -o $@ $(LDLIBS) $(LIBS)
	
main.o: main.cc Definitions.h Integrand.h HypGeo.h
	$(FC) -c $< -o $@  $(LIBS)

Integrand.o: Integrand.cc Definitions.h Integrand.h ResumFunc.h HypGeo.h

ResumFunc.o: ResumFunc.cc Definitions.h

HypGeo.o: HypGeo.cc Definitions.h

clean:
	rm -f SD-Angularity *.o 
