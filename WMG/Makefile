include ../common.mk

all: 
	${CXX} ${CXXFLAGS_EIGEN} ${IFLAGS_EIGEN} -c *.cpp
	${AR} -rc ../lib/libwmg.a *.o


clean:
	rm -f *.o


# dummy targets
.PHONY: clean
