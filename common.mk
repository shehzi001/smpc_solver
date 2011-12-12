CXX_WARN_FLAGS=-Wall -Wfloat-equal -Wshadow -pedantic -std=c++98
IFLAGS+=-I../include
LDFLAGS+=-L../lib/ -lsmpc_solver -lwmg

ifdef DEBUG
CXXFLAGS=-g -DSMPCS_DEBUG=1 ${CXX_WARN_FLAGS} ${IFLAGS}
CMAKEFLAGS=-DCMAKE_BUILD_TYPE=DEBUG
else
CXXFLAGS+=-O3 ${CXX_WARN_FLAGS} ${IFLAGS}
CMAKEFLAGS=-DCMAKE_BUILD_TYPE=Release
endif
