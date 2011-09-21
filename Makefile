smpc_solver:
	cd src; ${MAKE}

test: smpc_solver
	cd test; ${MAKE}

clean:
	cd test; ${MAKE} clean
	cd src; ${MAKE} clean
	rm -Rf docs

# dummy targets
.PHONY: clean
