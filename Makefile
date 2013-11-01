export PYTHONPATH := /home/eab/git/difwave/gsie/:$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

dirichlet-helmholtz.out: dirichlet-helmholtz.max
	maxima -b dirichlet-helmholtz.max > /dev/null

jacobian.out: jacobian.max
	maxima -b jacobian.max > /dev/null

singular.out: singular.max
	maxima -b singular.max > /dev/null

jacobian.gen:
	python dbsym.py

dbsym.o: dbsym.f90
	$(gf) -c toms_mod.f90 dbsym.f90

dbsym.so: dbsym.f90
	$(f2) -m dbsym -c toms_mod.f90 dbsym.f90

test: dbsym.so dbsym.o
	echo ok

clear:
	rm -f *.o *.so *.mod beta.f90 x.f90 jacobian.f90 jacobian.out
