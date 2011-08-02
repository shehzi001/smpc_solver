test: smpc_solver
	cd test; ${MAKE}

smpc_solver:
	cd src; ${MAKE}

clean:
	cd test; ${MAKE} clean
	cd src; ${MAKE} clean
	rm -Rf docs

# dummy targets
.PHONY: clean
