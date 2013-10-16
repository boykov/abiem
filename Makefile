export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(shell python defaults.py gsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py dotgsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py petsc4py_dir):$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

params.o: params.f90
	$(gf) -c params.f90
	ln -f params.mod dbsym/params.mod

phi.o: phi.f90 params.o
	$(gf) -c phi.f90

dbsym/dbsym.o: dbsym/dbsym.f90 dbsym/jacobian.out params.o
	cd dbsym && $(gf) -c dbsym.f90
	ln -f dbsym/dbsym.mod dbsym.mod

dbsym.so: dbsym/dbsym.f90 dbsym/jacobian.out params.o
	$(f2) -m dbsym -c params.o dbsym/dbsym.f90

dbsym/jacobian.out:
	cd dbsym && maxima -b /home/eab/git/difwave/gsie/jacobian.max > /dev/null
	cd dbsym && python dbsym.py

intrate.o: dbsym/dbsym.o intrate.f90 dbsym/jacobian.out
	$(gf) -c intrate.f90

libintphi.so: intphi.f90
	$(gf) -shared -c intphi.f90 -o libintphi.so

phi.so: phi.f90 params.o
	$(f2) -m phi -c params.o phi.f90

integ.so: integ.f90 intrate.o params.o phi.o libintphi.so
	test -s integ.so || $(f2) -m integ -lgomp -L. -lintphi -c params.o dbsym/dbsym.o intrate.o phi.o integ.f90 

test: testcase.py phi.so dbsym.so integ.so
	python -m unittest testcase$(tn)

clear:
	rm -f *.o *.so *.pyc *.tem *.mod dbsym/beta.f90 dbsym/x.f90 dbsym/jacobian.f90 dbsym/jacobian.out dbsym/*.mod dbsym/*.o
