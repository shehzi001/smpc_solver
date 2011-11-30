include ../common.mk


TESTS=\
	  test_01 \
	  test_02 \
	  test_03 \
	  test_04 \
	  test_05 \
	  test_06


all: ${TESTS} demo

${TESTS}:
	${CXX} ${CXXFLAGS} -c $@.cpp
	${CXX} -o $@.a $@.o ${LDFLAGS}

demo:
	${CXX} ${CXXFLAGS} -c $@.cpp
	${CXX} -o $@.a $@.o ${LDFLAGS}

clean:
	rm -f *.a *.o test_*.m

# dummy targets
.PHONY: clean
