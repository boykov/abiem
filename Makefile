export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(shell python defaults.py dbsym_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py gsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py dotgsie_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py petsc4py_dir):$(PYTHONPATH)
export PYTHONPATH := $(shell python defaults.py ahmed_dir):$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

export LD_LIBRARY_PATH=/home/eab/git/difwave/bie/

params.o: params.f90
	$(gf) -c params.f90

phi.o: phi.f90 params.o
	$(gf) -c phi.f90

dbsym/dbsym.o: dbsym/dbsym.f90
	cd dbsym && make dbsym.o

dbsym/fast_dbsym.o: dbsym/fast_dbsym.f90
	cd dbsym && make fast_dbsym.o

dbsym/dbsym.so: dbsym/dbsym.f90
	cd dbsym && make dbsym.so

libphi.so: phi.f90 params.o set_params.f90
	$(gf) -shared params.o phi.f90 -o libphi.so

phi.so: libphi.so
	test -s phi.so || f2py -m phi --overwrite-signature -h phi.pyf phi.f90
	test -s phi.so || $(f2) -m phi -L. -lphi -c phi.pyf phi.f90

libinteg.so: dbsym/dbsym.o dbsym/fast_dbsym.o integ.f90 params.o set_params.f90 libphi.so
	$(gf) -shared -I$(shell python defaults.py dbsym_dir) dbsym/dbsym.o dbsym/toms_mod.o dbsym/fast_dbsym.o params.o phi.o integ.f90 -o libinteg.so

integ.so: integ.f90 params.o phi.o libinteg.so
	test -s integ.so || f2py -m integ --overwrite-signature -h integ.pyf integ.f90
	test -s integ.so || $(f2) -m integ -lgomp -I$(shell python defaults.py dbsym_dir) -L. -linteg -c integ.pyf cover.f90

test:
	rm -f libinteg.so
	make testcase tn=.testBIEsmall

test2:
	rm -f libinteg.so
	@make testcase tn=.testBIEmedium

testcase: testcase.py phi.so integ.so
	python -m unittest testcase$(tn)

clear:
	rm -f *.o *.so *.pyc *.tem *.mod *.pyf
