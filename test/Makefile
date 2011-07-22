CXXFLAGS+=-O3 -Wall
IFLAGS+=-I../src/
OBJFILES=../src/chol_solve.o ../src/qp_as.o


all: test_01 test_02

test_01: 
	${CXX} ${CXXFLAGS} ${IFLAGS} ${OBJFILES} $@.cpp -o $@.a

test_02: 
	${CXX} ${CXXFLAGS} ${IFLAGS} ${OBJFILES} $@.cpp -o $@.a

clean:
	rm -f ./*.a

# dummy targets
.PHONY: clean
