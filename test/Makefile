CXXFLAGS+=-O3 -Wall
IFLAGS+=-I../src/
OBJFILES=../src/chol_solve.o ../src/qp_as.o


all: test_01 test_02

test_01: 
	${CXX} ${CXXFLAGS} ${IFLAGS} ${OBJFILES} ${.TARGET}.cpp -o ${.TARGET}.a

test_02: 
	${CXX} ${CXXFLAGS} ${IFLAGS} ${OBJFILES} ${.TARGET}.cpp -o ${.TARGET}.a

clean:
	rm -f ./*.a

# dummy targets
.PHONY: clean
