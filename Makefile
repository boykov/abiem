export PETSC_DIR=/home/eab/data/src/petsc-3.3-p4

FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

clear:
	rm -f intphi.f90 *.o *.so *.pyc *.tem *.mod

