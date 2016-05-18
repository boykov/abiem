dbsym_dir = $(shell python defaults.py dbsym_dir)
ahmed_dir = $(shell python defaults.py ahmed_dir)
sql_dir   = $(shell python defaults.py sql_dir)
petsc4py_dir = $(shell python defaults.py petsc4py_dir)
bie_dir = $(shell python defaults.py bie_dir)

export PYTHONPATH := ${petsc4py_dir}:$(sql_dir):$(dbsym_dir):$(PYTHONPATH)

gf = gfortran -fopenmp -ffree-line-length-none -fPIC -O3 -funroll-loops
f2 = f2py --f90flags="-ffree-line-length-none -fopenmp"

export LD_LIBRARY_PATH=${bie_dir}/src

include ahmed.mk

slaeahmed.cpp.o: slaeahmed.cpp bemmatrix.h
	/usr/bin/g++ -I. ${CXX_FLAGS} ${CXX_DIST_FLAGS} -o slaeahmed.cpp.o -c slaeahmed.cpp

slaeahmed.so: slaeahmed.cpp.o
	/usr/bin/g++   -fopenmp -DMETIS_VERSION=5  -std=c++0x -I/usr/local/include -pthread  -O3 -DNDEBUG slaeahmed.cpp.o  -o slaeahmed.so ${LINK_FLAGS} ${LINK_DIST_FLAGS}

solve:
	mpiexec.mpich -n 4 python mpislae.py -ksp_monitor

params.f90: common.py
	python -c "from common import common; c = common(); c.write_module()"

set_params.f90: common.py
	python -c "from common import common; c = common(); c.write_set_module()"

params.o: params.f90
	$(gf) -c params.f90

phi.o: phi.f90 params.o
	$(gf) -c phi.f90

dbsym.o: ${dbsym_dir}/dbsym.f90
	cd ${dbsym_dir} && make dbsym.o

fast_dbsym.o: ${dbsym_dir}/fast_dbsym.f90
	cd ${dbsym_dir} && make fast_dbsym.o

libphi.so: phi.f90 params.o set_params.f90
	$(gf) -shared params.o phi.f90 -o libphi.so

phi.so: libphi.so
	test -s phi.so || f2py -m phi --overwrite-signature -h phi.pyf phi.f90
	test -s phi.so || $(f2) -m phi -L. -lphi -c phi.pyf phi.f90

libinteg.so: dbsym.o fast_dbsym.o integ.f90 kernels.f90 params.o set_params.f90 libphi.so
	$(gf) -shared -I$(dbsym_dir) ${dbsym_dir}/dbsym.o ${dbsym_dir}/toms_mod.o ${dbsym_dir}/fast_dbsym.o params.o phi.o integ.f90 -o libinteg.so

integ.so: params.o phi.o libinteg.so
	test -s integ.so || f2py -m integ --overwrite-signature -h integ.pyf integ.f90
	test -s integ.so || $(f2) -m integ -lgomp -I$(dbsym_dir) -L. -linteg -c integ.pyf cover.f90

test:
	rm -f libinteg.so
	make testcase tn=.testBIEsmall

test3:
	rm -f libinteg.so
	@make testcase tn=.testBIEmedium

testcase: testcase.py phi.so integ.so slaeahmed.so
	python -m unittest testcase$(tn)

clear:
	rm -f `git ls-files -o`
