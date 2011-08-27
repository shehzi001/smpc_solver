CXX_WARN_FLAGS=-Wall -Wfloat-equal -Wshadow -pedantic

ifdef DEBUG
CXXFLAGS=-g -DSMPCS_DEBUG=1 ${CXX_WARN_FLAGS}
else
CXXFLAGS+=-O3 ${CXX_WARN_FLAGS}
endif
