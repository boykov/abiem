export PETSC_DIR = $(shell python defaults.py petsc_dir)

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py gsie_dir)
export PYTHONPATH := $(PYTHONPATH):$(shell python defaults.py dotgsie_dir)

test: testcase.py
	python -m unittest testcase.testEllipsoid$(tn)

clear:
	rm -f intphi.f90 *.o *.so *.pyc *.tem *.mod

