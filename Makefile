gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

specfun.out: specfun.max
	maxima -b specfun.max > /dev/null

dirichlet-helmholtz.out: dirichlet-helmholtz.max
	maxima -b dirichlet-helmholtz.max > /dev/null

jacobian.out: jacobian.max
	maxima -b jacobian.max > /dev/null

singular.out: singular.max
	maxima -b singular.max > /dev/null

spherical_bessel_j.f90: specfun.out
	python dbsym.py spec

xx_rho1.f90: jacobian.out
	python dbsym.py xx_rhopy

beta.f90: jacobian.out
	python dbsym.py jacobianpy

beta2.f90: jacobian.out
	python dbsym.py jacobian2py

singular.f90: singular.out
	python dbsym.py singularpy

Amn.f90: dirichlet-helmholtz.out
	python dbsym.py dirichlet_helmholtzpy

toms_mod.o: toms_mod.f90
	$(gf) -c toms_mod.f90

fast_dbsym.o: fast_dbsym.f90 spherical_bessel_j.f90
	$(gf) -c fast_dbsym.f90

dbsym.o: dbsym.f90 beta.f90 beta2.f90 singular.f90 Amn.f90 xx_rho1.f90 toms_mod.o
	$(gf) -c dbsym.f90

dbsym.so: dbsym.f90 beta.f90
	$(f2) -m dbsym -c toms_mod.f90 dbsym.f90

test: beta.f90 dbsym.so dbsym.o
	cd .. && make test

clear:
	rm -f `git ls-files -o`
