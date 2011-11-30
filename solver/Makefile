include ../common.mk

all: 
	${CXX} ${CXXFLAGS} ${IFLAGS} -c *.cpp
	${AR} -rc ../lib/libsmpc_solver.a *.o


clean:
	rm -f *.o


# dummy targets
.PHONY: clean
