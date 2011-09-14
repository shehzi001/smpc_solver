include ../common.mk

IFLAGS+=-I../src/
LDFLAGS+=-L../src/ -lsmpc_solver


TESTS=\
	  test_01 \
	  test_02 \
	  test_03 \
	  test_04 \
	  test_05 \
	  test_06 \
	  test_07


all: ${TESTS}

${TESTS}:
	${CXX} ${CXXFLAGS} ${IFLAGS} -c $@.cpp WMG.cpp init_WMG.cpp
	${CXX} -o $@.a $@.o WMG.o init_WMG.o ${LDFLAGS}

clean:
	rm -f *.a *.o

# dummy targets
.PHONY: clean
