include ../common.mk

all: 
	${CXX} ${CXXFLAGS} ${IFLAGS} -c *.cpp
	${AR} -rc libwmg.a *.o


clean:
	rm -f *.o *.a


# dummy targets
.PHONY: clean
