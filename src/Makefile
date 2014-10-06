OMP_NUM_THREADS=12

LIBS=-lz -ltbb
EIGEN_PATH=-I../lib/Eigen
ARCH=native

OPTFLAGS=-g2 -O3 -ffast-math -Wall -mfpmath=sse -msse2 -mssse3 \
	 -funroll-loops -fopenmp -march=$(ARCH)

INCLUDE_PATH=-I../lib -I../lib/imagelib -I../lib/matrix -I../lib/flann

LIB_PATH=-L../lib -L../lib/ann_1.1_char/lib

CC=g++

all: point_reduce

point_reduce.o: point_reduce.cpp
	$(CC) -c $(OPTFLAGS) $(EIGEN_PATH) $(INCLUDE_PATH) $(LIB_PATH) point_reduce.cpp 

keys.o: keys.cpp
	$(CC) -c $(OPTFLAGS) -lANN_char -lz keys.cpp

point_reduce: point_reduce.o keys.o
	$(CC) $(INCLUDE_PATH) $(LIB_PATH) -o $(OPTFLAGS) -o point_reduce point_reduce.o keys.o $(LIBS)

clean:
	rm -rf *.o point_reduce
