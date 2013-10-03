export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py gsie_dir)
export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py dotgsie_dir)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

params.o: params.f90
	$(gf) -c params.f90

phi.o: phi.f90 params.o
	$(gf) -c params.o phi.f90

intrate.o: intrate.f90 params.o
	$(gf) -c params.o intrate.f90

phi.so: phi.f90 params.o
	$(f2) -m phi -c params.o phi.f90

integ.so: integ.f90 intrate.o params.o phi.o
	$(f2) -m integ -lgomp -c params.o intrate.o phi.o integ.f90 

test: testcase.py phi.so integ.so
	python -m unittest testcase$(tn)

clear:
	rm -f intphi.f90 *.o *.so *.pyc *.tem *.mod

