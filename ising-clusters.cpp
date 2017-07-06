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

// This is nearly a copy of ising64.cpp and should be refactored out later when
// the classes are fixed.

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include <stdint.h>
using namespace std;
#include "MersenneTwister.h"

#ifdef PNG_DUMP
#include "lodepng.h"
#endif

#define WIDTH 2048 // Giant 131072 // Large 32768 // Med 4096 // Sm 1024
#define N (WIDTH*WIDTH) // Giant 17179869184 // Large 1073741824 // Med 16777216 // Sm 1048576
#define MAX_CLUSTER (64*WIDTH) // Giant 8388608 // Large 2097152 // Med 262144 // Sm 65536
//#define T 2.28
 #define T 2.269158

uint64_t tint;


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



int main(int argc, char *argv[])
{
	int i;
#ifdef PNG_DUMP
	int iframe=0;
#endif
	ising_class ising;
	ofstream metadatafile, datafile;

	cout << "Everyone's alive" << endl;

#ifdef SEED
    unsigned int seed = SEED;
#else
    unsigned int seed = ising.mt.randInt();
#endif
    ising.mt.seed(seed);

	// metadatafile.open("metadata.yaml");
	// metadatafile << "\"N\": " << N << endl << // N is false in yaml...
	// "\"T\": " << T  << endl <<
	// "Commit: " << COMMIT << endl << 
	// "Branch: " << BRANCH << endl <<
	// "Seed: " << seed << endl;
	// metadatafile.close();


//	ising.wolff(25);

	cout << "Loading..." << flush;
	ising.load_config("lattice.pos");
	cout << " Done." << endl;
	ising.clusters();
	
	ising.draw_xy_L(0,0, 1024, 1.0, 0);



	return 0;
}

ising_class::ising_class()
{
	int64_t i;

	spin = new unsigned char [N/8];
	incluster = new unsigned char [N/8];
	cluster = new uint64_t [MAX_CLUSTER];


	mag = -N;
	for (i=0;i<N/8;i++)
	{
		spin[i] = 0;
		incluster[i] = 0;
	}
	Nc = 0;
	cprob = 1 - exp(-2.0/T);

}

void ising_class::wolff(double sweeps)
{
	int64_t i,j,x,modi,current, nc, row, col;
	double t=0, teq = (sweeps)*((double) N);
	bool spin0;

//cprob=1.0;

	for (t=0; t<teq; t++)
	{
		// CLUSTER FLIP
		// choose start of cluster
		// if N>2^32 we need to randInts
#if N<4294967296
		x = mt.randInt(N-1);
#else
		x = mt.randInt(4294967295) + (mt.randInt(N/4294967296 - 1) * 4294967296);
#endif
		spin0 = getspin(x);

		// Add it to the cluster
		nc = 1;
		cluster[0] = x;
		flipincluster(x);

		i=0;
		while (i < nc)
		{
			row = cluster[(i%MAX_CLUSTER)]/WIDTH;
			col = cluster[(i%MAX_CLUSTER)] - row*WIDTH;

//cout << row << " " << col << endl;

			for (j=0;j<4;j++)
			{
				if (j < 2)
				{
					if (j==0)
					current = (col + 1)%WIDTH + row*WIDTH;
					else
					current = (col-1+WIDTH)%WIDTH + row*WIDTH;
				}
				else
				{
					if (j==2)
					current = col + ((row+1)%WIDTH)*WIDTH;
					else
					current = col + ((row-1+WIDTH)%WIDTH)*WIDTH;
				}

				// try to add all his neighbours to the stack


				if (!getincluster(current))
				{
					if (getspin(current) == spin0)
					{
						if (mt() < cprob)
						{
							cluster[(nc%MAX_CLUSTER)] = current;
							nc++;
							flipincluster(current);
						}
					}
				}

				if (nc - i > MAX_CLUSTER) {
					cout << "Cluster got too big" << endl;
					exit(1);
				}

			}


			i++;

		}
		/* Flip the cluster. We don't have it all stored so
		we'll have to form it again using the incluster array */

/*		modi=0;
		for (i=0;i<N;i++){
			modi += getincluster(i);}
*/

		// Start the cluster again
		nc = 1;
		cluster[0] = x;
		flipincluster(x);

		i=0;
		while (i < nc)
		{
			modi = i%MAX_CLUSTER;

			flipspin(cluster[modi]);

			row = cluster[modi]/WIDTH;
			col = cluster[modi] - row*WIDTH;

			for (j=0;j<4;j++)
			{
				if (j < 2)
				{
					if (j==0)
					current = (col + 1)%WIDTH + row*WIDTH;
					else
					current = (col-1+WIDTH)%WIDTH + row*WIDTH;
				}
				else
				{
					if (j==2)
					current = col + ((row+1)%WIDTH)*WIDTH;
					else
					current = col + ((row-1+WIDTH)%WIDTH)*WIDTH;
				}

				// try to add all his neighbours to the stack

				if (getincluster(current))
				{
					cluster[(nc%MAX_CLUSTER)] = current;
					flipincluster(current);
					nc++;
				}

			}

			i++;

			if (nc - i > MAX_CLUSTER) {
				cout << "Cluster got too big" << endl;
				exit(1);
			}

		}


		t+= nc;
		mag -= 2*nc*(2*spin0 - 1);

	}
}

void ising_class::save_config(const char * filename)
{
	ofstream lattice;

	lattice.open(filename,ofstream::binary);
	lattice.write((char *) spin,N/8);
	lattice.close();
}

void ising_class::load_config(const char * filename)
{
	ifstream lattice;

	lattice.open(filename,ifstream::binary);
	if (lattice.good())
		lattice.read((char *) spin,N/8);
	else
		cout << "Input file not found, starting with M=-N" << endl;
	lattice.close();

	magnetisation();
}

/*
void ising_class::draw_whole_lattice()
{
	long col, row;
	double spin;
	std::ostringstream filename;


	filename.str("");
	filename << "ising-whole.png";
	pngwriter png(WIDTH,WIDTH,0.0,filename.str().c_str());

	for (col=0;col<WIDTH;col++)
	{
		for (row=0;row<WIDTH;row++)
		{
			spin = (double) (getspin(row*WIDTH + col) + 1e-10);
			png.plot(col, row, spin, spin, spin);
		}
	}

	png.close();

}

void ising_class::draw_lattice(long Wframes, long block_length)
{
	long i, j, k, Nrn, Wrn, row, col;
	long i_fr, j_fr, kb, jb;
	long corner, bigcorner, Nblock;
	long Nframes = Wframes*Wframes;
	double av_spin;
	std::ostringstream filename;


	// Number of renormalised sites
	Wrn = (WIDTH/Wframes) / block_length;
	Nrn = Wrn * Wrn;
	Nblock = block_length * block_length;


	for (i_fr=0;i_fr<Wframes;i_fr++)
	{
	   for(j_fr=0;j_fr<Wframes;j_fr++)
	   {

		bigcorner = (j_fr*WIDTH*WIDTH)/Wframes + (i_fr*WIDTH)/Wframes;

		filename.str("");
		filename << "ising-" << i_fr << "x" << j_fr << "-";
		if (block_length < 10) filename << "0";
		if (block_length < 100) filename << "0";
		if (block_length < 1000) filename << "0";
		filename << block_length << ".png";
		pngwriter png(Wrn,Wrn,0.0,filename.str().c_str());

cout << filename.str().c_str() << endl;
		for (col=0;col<Wrn;col++)
		{
		    for (row=0;row<Wrn;row++)
		    {
		    	corner = bigcorner + col*block_length +
				row*block_length*WIDTH;

			// Block_length must divide evenly (no wrap around)
			av_spin = 0;
			for (jb=0;jb<block_length;jb++)
			   for (kb=0;kb<block_length;kb++)
			      av_spin += getspin(corner + WIDTH*jb + kb);

			av_spin /= Nblock;

			png.plot(col, row, av_spin, av_spin, av_spin);

		    }
		}

		png.close();
	   }
	}
}
*/
int64_t ising_class::magnetisation()
{
	int64_t i;

	for (mag=0, i=0;i<N;i++)
		mag += 2*getspin(i) - 1;

	return mag;
}

double ising_class::energy()
{
	int64_t i, row,col;
	double energy = 0;
	bool si;

	for (i=0;i<N;i++)
	{
		row = i/WIDTH;
		col = i - row*WIDTH;
		si = getspin(i);
		energy -= (2*si - 1) *
			(2*getspin(row*WIDTH + ((col+1)%WIDTH)) - 1);

		energy -= (2*si - 1) *
			(2*getspin(((row+1)%WIDTH)*WIDTH + col) - 1);
	}

	return energy;
}

inline bool ising_class::getspin(int64_t index)
{
	int shift = (index & 7);
	return ((spin[(index >> 3)] & (1 << shift)) >> shift);
}

inline bool ising_class::getincluster(int64_t index)
{
	int shift = (index & 7);
	return ((incluster[(index >> 3)] & (1 << shift)) >> shift);
}

inline void ising_class::flipspin(int64_t index)
{
	int shift = (index & 7);
	spin[(index >> 3)] ^= (1 << shift);
}

inline void ising_class::flipincluster(int64_t index)
{
	int shift = (index & 7);
	incluster[(index >> 3)] ^= (1 << shift);
}

void ising_class::clusters()
{
	int64_t i,j,x,current, nc, row, col;
	bool spin0;
	ofstream datafile;
//cprob=1.0;

	datafile.open("clusters.csv");

	for (x=0; x<N; x++)
	{
	if (!getincluster(x))
	{
		spin0 = getspin(x);

		// Add it to the cluster
		nc = 1;
		cluster[0] = x;
		flipincluster(x);

		i=0;
		while (i < nc)
		{
			row = cluster[(i%MAX_CLUSTER)]/WIDTH;
			col = cluster[(i%MAX_CLUSTER)] - row*WIDTH;

			for (j=0;j<4;j++)
			{
				if (j < 2)
				{
					if (j==0)
					current = (col + 1)%WIDTH + row*WIDTH;
					else
					current = (col-1+WIDTH)%WIDTH + row*WIDTH;
				}
				else
				{
					if (j==2)
					current = col + ((row+1)%WIDTH)*WIDTH;
					else
					current = col + ((row-1+WIDTH)%WIDTH)*WIDTH;
				}

				// try to add all his neighbours to the stack


				if (!getincluster(current))
				{
					if (getspin(current) == spin0)
					{
						cluster[(nc%MAX_CLUSTER)] = current;
						nc++;
						flipincluster(current);
					}
				}

				if (nc - i > MAX_CLUSTER) {
					cout << "Cluster got too big" << endl;
					exit(1);
				}
			}

			i++;

		}

		datafile << nc << endl;
	}
	}

	datafile.close();
}


#ifdef PNG_DUMP
void ising_class::draw_xy_L(long x, long y, long Wrn, double Lfac, int findex)
{
	// In this one we input L as a fraction of total length, from here we make our best guess for the
	// block size (rounding down where necessary).

	long i, j, k, Nrn, row, col;
	long kb, jb;
	long corner, bigcorner, Nblock;
	double av_spin, steep;
	char filename[20];
	int block_length;


	// Apply some contrast
	steep=6.0;

	// Rounding down
	block_length = (int) (WIDTH * Lfac / Wrn);
	// We might want to round Wrn up
	Wrn = lround(WIDTH * Lfac / block_length);

	// Number of renormalised sites
	Nrn = Wrn * Wrn;
	Nblock = block_length * block_length;

	cout << "Wrn=" << Wrn << ", block_length=" << block_length << endl;

	if ((x + Wrn*block_length <= WIDTH) && (y + Wrn*block_length <= WIDTH))
	{
		bigcorner = x + y*WIDTH;
        std::vector<unsigned char> pngbuf(4*Nrn);

		sprintf(filename, "swolff%4.4d.png", findex);

		cout << filename << endl;
		for (col=0;col<Wrn;col++)
		{
			for (row=0;row<Wrn;row++)
			{
				corner = bigcorner + col*block_length +
				row*block_length*WIDTH;

				// Block_length must divide evenly (no wrap around)
				av_spin = 0;
				for (jb=0;jb<block_length;jb++)
					for (kb=0;kb<block_length;kb++)
						av_spin += getspin(corner + WIDTH*jb + kb);

				av_spin /= Nblock;


				// av_spin = -pow(av_spin,2)*(2*av_spin-3);									// cubic (dy/dx=0 at x=0,1)
				av_spin = tanh(steep * (av_spin - 0.5)) / tanh(0.5 * steep) / 2.0 + 0.5;	// tanh

                pngbuf[4*(col + Wrn*row) + 0] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 1] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 2] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 3] = 255;


			}
		}

        //Encode the image
        unsigned error = lodepng::encode(filename, pngbuf, Wrn, Wrn);

        //if there's an error, display it
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
	}
	else
	{
		cout << "Haven't built PBC into this yet, choose different x,y" << endl;
		cout << "x=" << x << " y=" << y << " block_length=" << block_length << " Wrn=" << Wrn << endl;
	}

}

#endif
