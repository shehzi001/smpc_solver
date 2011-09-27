smpc_solver:
	cd solver; ${MAKE}

wmg:
	cd WMG; ${MAKE}

test: smpc_solver wmg
	cd test; ${MAKE}

clean:
	cd test; ${MAKE} clean
	cd solver; ${MAKE} clean
	cd WMG; ${MAKE} clean
	rm -Rf docs

# dummy targets
.PHONY: clean
