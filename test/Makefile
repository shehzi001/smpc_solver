include ../common.mk

IFLAGS+=-I../src/
LDFLAGS+=-L../src/ -lsmpc_solver


TESTS=\
	  test_01 \
	  test_02 \
	  test_03 \
	  test_04


all: ${TESTS} demo

${TESTS}:
	${CXX} ${CXXFLAGS} ${IFLAGS} -c $@.cpp WMG.cpp 
	${CXX} -o $@.a $@.o WMG.o ${LDFLAGS}

demo:
	${CXX} ${CXXFLAGS} ${IFLAGS} -c $@.cpp WMG.cpp 
	${CXX} -o $@.a $@.o WMG.o ${LDFLAGS}

clean:
	rm -f *.a *.o

# dummy targets
.PHONY: clean
