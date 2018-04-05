CXX_FLAGS = -I/usr/include/python2.7 \
            -fopenmp -DMETIS_VERSION=5  -std=c++0x \
            -I/usr/local/include \
            -pthread  -O3 -DNDEBUG \
            -I${ahmed_dir}/libsH/Include \
            -I${ahmed_dir}/basmod/Include \
            -I${ahmed_dir}/matrix/Include

LINK_DIST_FLAGS = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro

CXX_DIST_FLAGS = -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7

LINK_FLAGS = -rdynamic \
             -L/usr/include/python2.7/ -lpython2.7 \
             ${ahmed_dir}/libsH/build/libAHMED.so \
             ${ahmed_dir}/libsH/build/Basmod/libBASMOD.so \
             ${ahmed_dir}/libsH/build/Matrix/libMATRIX.so \
             -lblas -llapack -pthread \
             -L/usr/local/lib -lmpi_cxx -lmpi \
             -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl \
             -Wl,-rpath,/usr/local/lib:${ahmed_dir}/libsH/build:${ahmed_dir}/libsH/build/Basmod:${ahmed_dir}/libsH/build/Matrix
