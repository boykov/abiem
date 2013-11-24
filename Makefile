gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

dirichlet-helmholtz.out: dirichlet-helmholtz.max
	maxima -b dirichlet-helmholtz.max > /dev/null

jacobian.out: jacobian.max
	maxima -b jacobian.max > /dev/null

singular.out: singular.max
	maxima -b singular.max > /dev/null

beta.f90: singular.out jacobian.out dirichlet-helmholtz.out
	python dbsym.py

toms_mod.o: toms_mod.f90
	$(gf) -c toms_mod.f90

fast_dbsym.o: fast_dbsym.f90
	$(gf) -c fast_dbsym.f90

dbsym.o: dbsym.f90 beta.f90 toms_mod.o
	$(gf) -c dbsym.f90

dbsym.so: dbsym.f90 beta.f90
	$(f2) -m dbsym -c toms_mod.f90 dbsym.f90

test: beta.f90 dbsym.so dbsym.o
	cd .. && make test

clear:
	rm -f `git ls-files -o`
