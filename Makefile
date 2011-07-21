test: smpc_solver
	cd test; make

smpc_solver:
	cd src; make

clean:
	cd test; make clean
	cd src; make clean

# dummy targets
.PHONY: clean
