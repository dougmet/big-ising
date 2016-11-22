// Ising model designed to work around Tc

#include <iostream>
#include <fstream>

#include <sstream>
#include <math.h>
#include <string>
#include <stdint.h>
using namespace std;
#include "MersenneTwister.h"

#include "lodepng.h"

#define WIDTH 2048 // Giant 131072 // Large 32768 // Med 4096 // Sm 1024
#define N (WIDTH*WIDTH) // Giant 17179869184 // Large 1073741824 // Med 16777216 // Sm 1048576
// For MAX_CLUSTER YOU MUST NOT PUT SOMETHING LIKE 64*WIDTH, PUT THE NUMBER
#define MAX_CLUSTER 65536 // Large 2097152 // Med 262144 // Sm 65536
#define T 2.2694
// #define T 2.269158

long tint;

class ising_class
{
	public:

	unsigned char *spin;
	unsigned char *incluster;
	unsigned long *cluster;
	long Nc;
	long mag;
	double cprob;
	MTRand mt;

	ising_class();
	inline bool getspin(long index);
/*	inline bool getincluster(long index);
	inline void flipspin(long index);
	inline void flipincluster(long index);

	void wolff(double sweeps);
 */
	double energy();
	long magnetisation();
	void draw_lattice(long Wframes, long block_length);
	void draw_xy(long x, long y, long Wrn, long block_length);
	void draw_xy_L(long x, long y, long Wrn, double Lfac, int findex);
	void draw_whole_lattice();
	void save_config(const char * filename);
	void load_config(const char * filename);
};



int main(int argc, char *argv[])
{
	long i, Nframes, x, y;
	double Lfac, scale, cx, cy;
	long inmag, inx, iny;
	char contin = 'y';
	ising_class ising;
	ofstream clearfile, datafile;

	// Check inputs are okay
/*	if (argc <4)
	{
		cout << "Not enough inputs, need x, y (0..1) and magnification" << endl;
		exit(1);
	}

	inx = argv[1];
	iny = argv[2];
	inmag = argv[3];
*/

	cout << "Everyone's alive" << endl;

	cout << "Loading..."; cout.flush();
	ising.load_config("lattice.pos");
	cout << " done." << endl;

	while (contin == 'y' || contin == 'Y')
	{
		cout << "Enter x: "; cout.flush();
		cin >> inx;

		cout << "Enter y: "; cout.flush();
		cin >> iny;

		cout << "Enter length (max=" << WIDTH << "): "; cout.flush();
		cin >> inmag;

		ising.draw_xy_L(inx, iny, 1024, ((double) inmag)/WIDTH, inmag);

		cout << "Done. Continue? (y/n) "; cout.flush();
		cin >> contin; cout << endl;
	}



//	ising.draw_lattice(4, 2);
//	ising.draw_lattice(1, 1);

	cx=0.5;
	cy=0.66;
	Lfac = 1.0;
	Nframes=1000;
	scale = pow(1.0/128.0,1.0/Nframes);
	for (i=0;i<Nframes;i++)
	{
		// First x
		if (cx - Lfac/2 < 0)
			x=0;
		else if (cx + Lfac/2 > 1)
			x = (int) (WIDTH * (1 - Lfac));
		else
			x = (int) (WIDTH * (cx - Lfac/2));

		// Now y
		if (cy - Lfac/2 < 0)
			y=0;
		else if (cy + Lfac/2 > 1)
			y = (int) (WIDTH * (1 - Lfac));
		else
			y = (int) (WIDTH * (cy - Lfac/2));


		ising.draw_xy_L(x,y, 720, Lfac, i);
		Lfac *= scale;
	}

//	ising.draw_whole_lattice();



	return 0;
}

ising_class::ising_class()
{

	spin = new unsigned char [N/8];

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
	lattice.read((char *) spin,N/8);
	lattice.close();
}


void ising_class::draw_whole_lattice()
{
	long col, row;
	double spin;
	std::ostringstream filename;


	filename.str("");
	filename << "ising-whole.png";
    std::vector<unsigned char> pngbuf(4*N);

	for (col=0;col<WIDTH;col++)
	{
		for (row=0;row<WIDTH;row++)
		{
			spin = (double) (getspin(row*WIDTH + col) + 1e-10);
			pngbuf[4*(col + WIDTH*row) + 0] = 255*spin;
			pngbuf[4*(col + WIDTH*row) + 1] = 255*spin;
			pngbuf[4*(col + WIDTH*row) + 2] = 255*spin;
			pngbuf[4*(col + WIDTH*row) + 3] = 255;
		}
	}

	unsigned error = lodepng::encode(filename.str().c_str(), pngbuf, WIDTH, WIDTH);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

}

void ising_class::draw_xy(long x, long y, long Wrn, long block_length)
{
	long i, j, k, Nrn, row, col;
	long kb, jb;
	long corner, bigcorner, Nblock;
	double av_spin, steep;
	std::ostringstream filename;


	// Number of renormalised sites
	Nrn = Wrn * Wrn;
	Nblock = block_length * block_length;


	if ((x + Wrn*block_length <= WIDTH) && (y + Wrn*block_length <= WIDTH))
	{
		bigcorner = x + y*WIDTH;

		filename.str("");
		filename << "ising-" << x << "-" << y << "-";
		if (block_length < 10) filename << "0";
		if (block_length < 100) filename << "0";
		if (block_length < 1000) filename << "0";
		filename << block_length << ".png";

		std::vector<unsigned char> pngbuf(4*Nrn);

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


				// av_spin = -pow(av_spin,2)*(2*av_spin-3);									// cubic (dy/dx=0 at x=0,1)
				av_spin = tanh(steep * (av_spin - 0.5)) / tanh(0.5 * steep) / 2.0 + 0.5;	// tanh

                pngbuf[4*(col + Wrn*row) + 0] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 1] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 2] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 3] = 255;
			}
		}

		//Encode the image
        unsigned error = lodepng::encode(filename.str().c_str(), pngbuf, Wrn, Wrn);

        //if there's an error, display it
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
	}
	else
	{
		cout << "Haven't built PBC into this yet, choose different x,y" << endl;
		cout << "x=" << x << " y=" << y << " block_length=" << block_length << " Wrn=" << Wrn << endl;
	}

}


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
	cout << Wrn <<endl;
	while ((x + Wrn*block_length > WIDTH) || (y + Wrn*block_length > WIDTH))
		Wrn --;


	// Number of renormalised sites
	Nrn = Wrn * Wrn;
	Nblock = block_length * block_length;


	if ((x + Wrn*block_length <= WIDTH) && (y + Wrn*block_length <= WIDTH))
	{
		cout << "x=" << x << " y=" << y << " block_length=" << block_length << " Wrn=" << Wrn << endl;
		bigcorner = x + y*WIDTH;

		sprintf(filename, "frame%4.4d.png", findex);
		std::vector<unsigned char> pngbuf(4*Nrn);

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


void ising_class::draw_lattice(long Wframes, long block_length)
{
	long i, j, k, Nrn, Wrn, row, col;
	long i_fr, j_fr, kb, jb;
	long corner, bigcorner, Nblock;
	long Nframes = Wframes*Wframes;
	double av_spin, steep;
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
		std::vector<unsigned char> pngbuf(4*Nrn);

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

				// Apply some contrast
				steep=6.0;
				// av_spin = -pow(av_spin,2)*(2*av_spin-3);									// cubic (dy/dx=0 at x=0,1)
				av_spin = tanh(steep * (av_spin - 0.5)) / tanh(0.5 * steep) / 2.0 + 0.5;	// tanh

				pngbuf[4*(col + Wrn*row) + 0] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 1] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 2] = lround(255*av_spin);
                pngbuf[4*(col + Wrn*row) + 3] = 255;

		    }
		}

		//Encode the image
        unsigned error = lodepng::encode(filename.str().c_str(), pngbuf, Wrn, Wrn);

        //if there's an error, display it
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
	   }
	}
}

long ising_class::magnetisation()
{
	long i, mag=0;

	for (i=0;i<N;i++)
		mag += 2*getspin(i) - 1;

	return mag;
}

double ising_class::energy()
{
	long i, row,col;
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

inline bool ising_class::getspin(long index)
{
	int shift = (index & 7);
	return ((spin[(index >> 3)] & (1 << shift)) >> shift);
}
