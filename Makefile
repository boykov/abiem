export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(shell python defaults.py dbsym_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py gsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py dotgsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py petsc4py_dir):$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

params.o: params.f90
	$(gf) -c params.f90

phi.o: phi.f90 params.o
	$(gf) -c phi.f90

dbsym/dbsym.o:
	cd dbsym && make dbsym.o

dbsym/dbsym.so:
	cd dbsym && make dbsym.so

intrate.o: dbsym/dbsym.o intrate.f90
	$(gf) -c -I$(shell python defaults.py dbsym_dir) intrate.f90

libintphi.so: intphi.f90
	$(gf) -shared -c intphi.f90 -o libintphi.so

phi.so: phi.f90 params.o
	$(f2) -m phi -c params.o phi.f90

integ.so: integ.f90 intrate.o params.o phi.o libintphi.so
	test -s integ.so || $(f2) -m integ -lgomp -L. -lintphi -c params.o dbsym/dbsym.o intrate.o phi.o integ.f90 

test: testcase.py phi.so dbsym/dbsym.so integ.so
	python -m unittest testcase$(tn)

clear:
	rm -f *.o *.so *.pyc *.tem *.mod
