include ./common.mk

libs: smpc_solver wmg

smpc_solver:
	cd solver; ${MAKE}

wmg:
	cd WMG; ${MAKE}

test: smpc_solver wmg
	cd test; ${MAKE}

cmake: 
	-mkdir build;
ifdef TOOLCHAIN
	cd build; cmake -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN} ${CMAKEFLAGS} ..;
else
	cd build; cmake ${CMAKEFLAGS} ..;
endif
	cd build; ${MAKE}

clean:
	cd test; ${MAKE} clean
	cd solver; ${MAKE} clean
	cd WMG; ${MAKE} clean
	rm -f docs/*.html
	rm -f docs/*.png
	rm -f docs/*.css
	rm -f docs/*.js
	rm -f docs/formula.repository
	rm -Rf build/*
	rm -f lib/*.a

copydoc: clean
	doxygen
	rm ../smpc_solver_docs/*
	cp docs/*.html docs/*.png docs/*.css ../smpc_solver_docs/

# dummy targets
.PHONY: clean
