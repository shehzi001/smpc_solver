smpc_solver:
	cd solver; ${MAKE}

wmg:
	cd WMG; ${MAKE}

test: smpc_solver wmg
	cd test; ${MAKE}

cmake: 
	-mkdir build;
ifdef TOOLCHAIN
	cd build; cmake -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN} ..;
else
	cd build; cmake ..;
endif
	cd build; ${MAKE}

clean:
	cd test; ${MAKE} clean
	cd solver; ${MAKE} clean
	cd WMG; ${MAKE} clean
	rm -f docs/*.html
	rm -f docs/*.png
	rm -f docs/*.css
	rm -f docs/formula.repository
	rm -Rf build/*
	rm -f lib/*.a

copydoc: clean
	doxygen
	rm ../smpc_solver_docs/*
	cp docs/*.html docs/*.png docs/*.css ../smpc_solver_docs/

# dummy targets
.PHONY: clean
