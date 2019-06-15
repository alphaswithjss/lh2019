# dummy Makefile which just invokes scons (a more modern build tool)

all: libCmdLine.a

CXXFLAGS=-g -ansi -pedantic -Wall -O3 -fPIC -DPIC -std=c++11

libCmdLine.a: CmdLine.o
	ar rc libCmdLine.a CmdLine.o
	ranlib libCmdLine.a

clean:
	rm -f *.o

CmdLine.o: CmdLine.cc CmdLine.hh
