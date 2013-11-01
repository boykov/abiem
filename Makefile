export PYTHONPATH := /home/eab/git/difwave/gsie/:$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

jacobian.out:
	maxima -b /home/eab/git/difwave/gsie/jacobian.max > /dev/null
	python dbsym.py

dbsym.o: dbsym.f90 jacobian.out
	$(gf) -c dbsym.f90

dbsym.so: dbsym.f90 jacobian.out
	$(f2) -m dbsym -c dbsym.f90

test: dbsym.so dbsym.o
	echo ok

clear:
	rm -f *.o *.so *.mod beta.f90 x.f90 jacobian.f90 jacobian.out
