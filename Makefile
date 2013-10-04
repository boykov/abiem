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

phi.so: phi.f90 params.o
	$(f2) -m phi -c params.o phi.f90

test: testcase.py phi.so
	python -m unittest testcase$(tn)

clear:
	rm -f *.o *.so *.pyc *.tem *.mod

