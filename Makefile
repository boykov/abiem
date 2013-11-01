export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py gsie_dir)
export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py dotgsie_dir)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

export LD_LIBRARY_PATH=/home/eab/git/difwave/bie/

params.o: params.f90
	$(gf) -c params.f90

libphi.so: phi.f90 params.o set_params.f90
	$(gf) -shared params.o phi.f90 -o libphi.so

phi.so: libphi.so
	test -s phi.so || f2py -m phi --overwrite-signature -h phi.pyf phi.f90
	test -s phi.so || $(f2) -m phi -L. -lphi -c phi.pyf phi.f90

test: testcase.py libphi.so phi.so
	python -m unittest testcase$(tn)

clear:
	rm -f *.o *.so *.pyc *.tem *.mod *.pyf

