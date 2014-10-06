/* 
 *  Copyright (c) 2008  Noah Snavely (snavely (at) cs.washington.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* keys.h */
/* Class for SIFT keypoints */

#include <vector>
#include <stdio.h>
#include <zlib.h>

//#include "ANN/ANN.h"
//using namespace ann_1_1_char;

class Keypoint {
public:    
    Keypoint(float x, float y, float scale, float ori, short int *d) :
	m_x(x), m_y(y), m_scale(scale), m_ori(ori), m_d(d)
    { }

    float m_x, m_y;        /* Subpixel location of keypoint. */
    float m_scale, m_ori;  /* Scale and orientation (range [-PI,PI]) */
    short int *m_d;    /* Vector of descriptor values */
};


/* Data struct for matches */
class KeypointMatch {
public:
    KeypointMatch() 
    { }

    KeypointMatch(int idx1, int idx2) :
	m_idx1(idx1), m_idx2(idx2)
    { }

    int m_idx1, m_idx2;
    // float m_x1, m_y1;
    // float m_x2, m_y2;
};

typedef struct {
    float x, y;
    float scale;
    float orient;
} keypt_t;

/* Returns the number of keys in a file */
int GetNumberOfKeys(const char *filename);

/* This reads a keypoint file from a given filename and returns the list
 * of keypoints. */
int ReadKeyFile(const char *filename, unsigned char **keys, 
                keypt_t **info = NULL);

int ReadKeyPositions(const char *filename, keypt_t **info);

/* Read keypoints from the given file pointer and return the list of
 * keypoints.  The file format starts with 2 integers giving the total
 * number of keypoints and the size of descriptor vector for each
 * keypoint (currently assumed to be 128). Then each keypoint is
 * specified by 4 floating point numbers giving subpixel row and
 * column location, scale, and orientation (in radians from -PI to
 * PI).  Then the descriptor vector for each keypoint is given as a
 * list of integers in range [0,255]. */
int ReadKeys(FILE *fp, unsigned char **keys, keypt_t **info = NULL);
int ReadKeysGzip(gzFile fp, unsigned char **keys, keypt_t **info = NULL);

