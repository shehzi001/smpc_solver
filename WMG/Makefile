include ../common.mk

CXXFLAGS+=-I../solver/

all: 
	${CXX} ${CXXFLAGS} ${IFLAGS} -c *.cpp
	${AR} -rc ../lib/libwmg.a *.o


clean:
	rm -f *.o


# dummy targets
.PHONY: clean
