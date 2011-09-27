include ../common.mk

CXXFLAGS+=-I../solver/

all: 
	${CXX} ${CXXFLAGS} ${IFLAGS} -c *.cpp
	${AR} -rc libwmg.a *.o


clean:
	rm -f *.o *.a


# dummy targets
.PHONY: clean
