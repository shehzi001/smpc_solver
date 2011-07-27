#CXXFLAGS=-g -Wall
CXXFLAGS+=-O3 -Wall
IFLAGS+=-I../src/
OBJFILES=../src/chol_solve.o ../src/qp_as.o
COMPILER=${CXX} ${CXXFLAGS} ${IFLAGS} ${OBJFILES}


all: test_01 test_02 test_03

test_01: 
	${COMPILER} $@.cpp -o $@.a

test_02: 
	${COMPILER} $@.cpp -o $@.a

test_03: 
	${COMPILER} $@.cpp -o $@.a

clean:
	rm -f ./*.a

# dummy targets
.PHONY: clean
