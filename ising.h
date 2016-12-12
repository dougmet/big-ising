// MIT License
//
// Copyright (c) 2016 Douglas Ashton
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Ising model designed to work around Tc

#ifndef ISING_H // header protection
#define ISING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include <stdint.h>
#include "MersenneTwister.h"

#ifdef PNG_DUMP
#include "lodepng.h"
#endif

// Moved these to main.cpp
//#define WIDTH 2048 // Giant 131072 // Large 32768 // Med 4096 // Sm 1024
//#define N (WIDTH*WIDTH) // Giant 17179869184 // Large 1073741824 // Med 16777216 // Sm 1048576
//#define MAX_CLUSTER (64*WIDTH) // Giant 8388608 // Large 2097152 // Med 262144 // Sm 65536
//#define T 2.28
//#define T 2.269158
//uint64_t tint;

class ising_class
{
	public:

	unsigned char *spin;
	unsigned char *incluster;
	uint64_t *cluster;
	int64_t Nc;
	int64_t mag;
	double cprob;
	MTRand mt;

	ising_class();
	inline bool getspin(int64_t index);
	inline bool getincluster(int64_t index);
	inline void flipspin(int64_t index);
	inline void flipincluster(int64_t index);

	void wolff(double sweeps);
	double energy();
	int64_t magnetisation();
/*	void draw_lattice(long Wframes, long block_length);
	void draw_whole_lattice(); */
	void clusters();
	void save_config(const char * filename);
	void load_config(const char * filename);

#ifdef PNG_DUMP
	void draw_xy_L(long x, long y, long Wrn, double Lfac, int findex);
#endif
};

#endif // end header protection
