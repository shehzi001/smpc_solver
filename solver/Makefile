include ../common.mk

all: 
	echo "#define HAVE_FEENABLEEXCEPT" > solver_config.h
	${CXX} ${CXXFLAGS} ${IFLAGS} -c *.cpp
	${AR} -rc ../lib/libsmpc_solver.a *.o


clean:
	rm -f *.o
	rm solver_config.h


# dummy targets
.PHONY: clean
