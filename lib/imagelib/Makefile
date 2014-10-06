# Makefile for imagelib

MACHTYPE=$(shell uname -m)

ifeq ($(MACHTYPE), x86_64)
P=x64
else
P=ia32
endif

GCC			= g++

CC=gcc
# OPTFLAGS=-g2
OPTFLAGS=-O3 -g2 -ffast-math -Wall -mfpmath=sse -msse2 -funroll-loops -march=core2

ifdef PROFILE
PROFILE_FLAGS=-pg
else
PROFILE_FLAGS=
endif

OTHERFLAGS=-Wall -Wno-unused-result

IMAGELIB_OBJS= affine.o bmp.o canny.o color.o	\
	fileio.o filter.o fit.o fmatrix.o		\
	homography.o horn.o image.o lerp.o		\
	morphology.o pgm.o poly.o qsort.o		\
	resample.o tps.o transform.o triangulate.o	\
	util.o kruppa.o vanish.o

INCLUDE_PATH=-I../matrix

CFLAGS=$(OTHERFLAGS) $(OPTFLAGS) $(PROFILE_FLAGS) $(INCLUDE_PATH)

IMAGELIB=libimage.a

%.o : %.c
	$(CC) -c -o $@ $(CFLAGS) $<

all: $(IMAGELIB)

$(IMAGELIB): $(IMAGELIB_OBJS)
	ar r $(IMAGELIB) $(IMAGELIB_OBJS)
	cp $(IMAGELIB) ..

clean:
	rm -f *.o *~ $(IMAGELIB)
