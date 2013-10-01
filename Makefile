export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py gsie_dir)
export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py dotgsie_dir)

phi.so: phi.f90
	f2py -m phi --f90flags="-ffree-line-length-none -fopenmp" -c phi.f90

intrate.o: intrate.f90
	gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops -c params.f90 intrate.f90

integ.so: integ.f90 intrate.o
	f2py -m integ --f90flags="-ffree-line-length-none -fopenmp" -c params.f90 intrate.o integ.f90 

test: testcase.py phi.so
	python -m unittest testcase$(tn)

clear:
	rm -f intphi.f90 *.o *.so *.pyc *.tem *.mod

