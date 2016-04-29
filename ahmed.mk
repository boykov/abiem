CXX_FLAGS = -I/usr/include/python2.7 -fopenmp -DMETIS_VERSION=5  -std=c++0x -I/usr/local/include -pthread  -O3 -DNDEBUG -I/home/eab/data/gitno/github/AHMED/examples/Include -I/home/eab/data/gitno/github/AHMED/libsH/Include -I/home/eab/data/gitno/github/AHMED/basmod/Include -I/home/eab/data/gitno/github/AHMED/matrix/Include -I/usr/local/include

LINK_DIST_FLAGS = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro

CXX_DIST_FLAGS = -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7

LINK_FLAGS = -rdynamic -L/usr/include/python2.7/ -lpython2.7 /home/eab/data/gitno/github/AHMED/examples/build/libsH/libAHMED.so /home/eab/data/gitno/github/AHMED/examples/build/Basmod/libBASMOD.so /home/eab/data/gitno/github/AHMED/examples/build/Matrix/libMATRIX.so -lblas -llapack /home/eab/data/gitno/github/AHMED/examples/build/Basmod/libBASMOD.so -lblas -llapack -pthread -L/usr/local/lib -lmpi_cxx -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl -Wl,-rpath,/home/eab/data/gitno/github/AHMED/examples/build/libsH:/home/eab/data/gitno/github/AHMED/examples/build/Basmod:/home/eab/data/gitno/github/AHMED/examples/build/Matrix
