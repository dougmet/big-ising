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

#define WIDTH 2048 // Giant 131072 // Large 32768 // Med 4096 // Sm 1024
#define N (WIDTH*WIDTH) // Giant 17179869184 // Large 1073741824 // Med 16777216 // Sm 1048576
#define MAX_CLUSTER (64*WIDTH) // Giant 8388608 // Large 2097152 // Med 262144 // Sm 65536
//#define T 2.28
#define T 2.269158

//uint64_t tint;

int main(int argc, char *argv[])
{
	int i, iframe=0;
	bool x;
	ising_class ising;
	ofstream clearfile, datafile;


	cout << "Everyone's alive" << endl;

	clearfile.open("data");
	clearfile << "N " << N << " T " << T  << endl;
	clearfile.close();

	//ising.mt.seed(10);
//	ising.wolff(25);

	cout << "Loading..." << flush;
	ising.load_config("lattice.pos");
	cout << " Done." << endl;

	datafile.open("data",fstream::app | fstream::ate);
	i=0;
	while (i>-1)// || fabs(ising.mag)/((double) N) > 0.1)
	{
	cout << i << endl;
		ising.wolff(0.1);
#ifdef PNG_DUMP
		ising.draw_xy_L(0,0, 720, 1.0, iframe++);
#endif
		datafile << ising.mag << " " << ising.energy() << endl;

		if ((fabs(ising.mag)/((double) N) < 0.03) && (i>50))
		{
			ising.save_config("lattice.pos");
			cout << "Close to zero magnetisation. Stopping." << endl;
			exit(0);
		}

		if (i%20==0)
			ising.save_config("lattice.pos");

		i++;
	}
	datafile.close();




//	ising.load_config("lattice1.pos");
//	ising.clusters();

/*	ising.draw_lattice(4, 2);
	ising.draw_lattice(1, 4);
	ising.draw_whole_lattice();
*/


	return 0;
}
