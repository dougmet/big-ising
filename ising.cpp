// Ising model designed to work around Tc

#include <iostream>
#include <fstream> 
#include <stdlib.h>
#include <sstream>      
#include <math.h>  
#include <string> 
using namespace std; 
#include "/home/da246/MersenneTwister.h"

/*#include "/home/douglas/cprogs/pngwriter.h"
#include "/home/douglas/cprogs/pngwriter.cc"
*/
#define N 17179869184 // Large 1073741824 // Med 16777216 // Sm 1048576
#define WIDTH  131072 // Large 32768 // Med 4096 // Sm 1024
// For MAX_CLUSTER YOU MUST NOT PUT SOMETHING LIKE 64*WIDTH, PUT THE NUMBER
#define MAX_CLUSTER 8388608 // Large 2097152 // Med 262144 // Sm 65536
//#define T 2.28
 #define T 2.269158

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
	inline bool getincluster(long index);
	inline void flipspin(long index);
	inline void flipincluster(long index);
	
	void wolff(double sweeps);
	double energy();
	long magnetisation();
/*	void draw_lattice(long Wframes, long block_length);
	void draw_whole_lattice(); */
	void clusters();
	void save_config(const char * filename);
	void load_config(const char * filename);
};



int main(int argc, char *argv[])
{
	long i;
	bool x;
	ising_class ising;
	ofstream clearfile, datafile;


	cout << "Everyone's alive" << endl;

	clearfile.open("data");
	clearfile << "N " << N << " T " << T  << endl;
	clearfile.close();
	
	//ising.mt.seed(10);
//	ising.wolff(25);
	
//	ising.load_config("lattice.pos");

	datafile.open("data",fstream::app | fstream::ate);
	i=0;
	while (i<200)// || fabs(ising.mag)/((double) N) > 0.1)
	{
	cout << i << endl;
		ising.wolff(0.1);
		datafile << ising.mag << " " << ising.energy() << endl;
		
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

ising_class::ising_class()
{
	long i;
	
	spin = new unsigned char [N/8];
	incluster = new unsigned char [N/8];
	cluster = new unsigned long [MAX_CLUSTER];
	
	
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
	long i,j,x,modi,current, nc, row, col;
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
		
		modi=0;
		for (i=0;i<N;i++){
			modi += getincluster(i);}
			
		
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
	lattice.read((char *) spin,N/8);
	lattice.close();
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

inline bool ising_class::getincluster(long index)
{
	int shift = (index & 7);
	return ((incluster[(index >> 3)] & (1 << shift)) >> shift);
}

inline void ising_class::flipspin(long index)
{
	int shift = (index & 7);
	spin[(index >> 3)] ^= (1 << shift);
}

inline void ising_class::flipincluster(long index)
{
	int shift = (index & 7);
	incluster[(index >> 3)] ^= (1 << shift);
}

void ising_class::clusters()
{
	long i,j,x,modi,current, nc, row, col;
	bool spin0;
	ofstream datafile;
//cprob=1.0;
	
	datafile.open("clusters.dat");

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
