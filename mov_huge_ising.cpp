// Given a configuration and a starting point, zoom out and generate frames

#include <iostream>
#include <fstream> 

#include <sstream>      
#include <math.h>  
#include <string>
#include <float.h>

#ifdef __APPLE__
 #include "/Users/da246/Cprogs/pngwriter.h"
 #include "/Users/da246/Cprogs/pngwriter.cc"
#else
 #include "/home/da246/pngwriter.h"
 #include "/home/da246/pngwriter.cc"
#endif


using namespace std; 
#include "../MersenneTwister.h"



#define DIMENSION 2
#define WIDTH 131072
#define N 17179869184
#define MAXNP 1073741824 // 2^30
#define L 2048.0
#define T 0.4814
#define MU 1.66

#define EPS 1e-14
#define LEPS (L*(1-FLT_EPSILON))


// For the visualisation this is no longer the real cutoff. It's the size of coarse graining.
#define RCUT 2.5

class ising_class
	{
	public:
		
		unsigned char *spin;
		
		ising_class();
		void draw_xy_L(long x, long y, long Wrn, double Lfac, int findex);
		inline bool getspin(long index);
		void load_config(const char * filename);
	};

class particle_class
	{
	public:
		float r[DIMENSION];
		int next;
	};



class box_class
{
public:
	
	int Np;						// Number of active particles and queued particles in slabs
	int Wcells, Ncells;			// Linear number of cells and total number
	int *cell;					// The cells linked list
	float vol;
	double cell_length;
	
	particle_class *particle;	// The particles
	
	// Methods
	box_class();		// Constructor
	void start_job();	// Where we'll do everything
	void add_particle(int in_particle, int ic);
	int assign_cell(int inpart);
	void load_positions(ising_class *ising, double cx, double cy);
	
	
	// Just for this programme
	double startx, starty, Lwin;
	int resolution;
	
};

						

int main(int argc, char *argv[])
{
	int i, ix, iy;
	box_class *box;
	ising_class *ising;
	
	double startx, starty, Lwin;
	int Nf, t, resolution;
	double cx,cy;
	char convert[50];
	
	ising = new ising_class;
	cout << "Everyone's alive. Loading..." << flush;
	ising->load_config("lattice.pos");
	cout << " done." << endl;
	
	box = new box_class();
	
	
	if (argc<5)
	{
		cout << "Need x,y,L,resolution." << endl;
		exit(1);
	}
	startx = atof(argv[1]);
	starty = atof(argv[2]);
	Lwin = atof(argv[3]);
	resolution = atoi(argv[4]);
	
	cout << "Loading particles..." << flush;
	box->load_positions(ising, startx, starty);
	
	
	Nf=1000;
	
	for (t=0;t<Nf;t++)
	{
		// Decide on the parameters		
		box->Lwin = Lwin * pow(WIDTH/Lwin,((double) t)/((double) Nf-1));
		if (startx > WIDTH/2.0)
			cx = startx+1 - pow(1+startx-WIDTH/2.0,((double) t)/((double) Nf-1));
		else
			cx = startx-1 + pow(1+WIDTH/2.0-startx,((double) t)/((double) Nf-1));
		
		if (starty > WIDTH/2.0)
			cy = starty+1 - pow(1+starty-WIDTH/2.0,((double) t)/((double) Nf-1));
		else
			cy = starty-1 + pow(1+WIDTH/2.0-starty,((double) t)/((double) Nf-1));
		
		cout << "Global: cx" << cx << " cy" << cy << endl;
		
		// Decide if we'll use the blocking (large Lwin) or squares (small Lwin)
		if ((box->Lwin < L*0.9) && (cx - startx + L/2 - box->Lwin/2 > 0) && (cy - starty + L/2 - box->Lwin/2 > 0))
		{
			// coordinates are reversed
			cx = cx - startx + L/2;
			cy = cy - starty + L/2;
			box->startx = cx - box->Lwin/2;
			box->starty = cy - box->Lwin/2;
			box->resolution = resolution;
				
			cout << "Local sq: " << box->startx << ", " << box->starty << ": " << box->Lwin << endl;
		
			box->start_job();
		
			sprintf(convert, "convert -resize 1024x1024 squares.png squares%4.4d.jpg", t);
			printf("%s\n\n", convert);
			system(convert);
		}
		else
		{
			ix = (int) (cx - box->Lwin/2);
			iy = (int) (cy - box->Lwin/2);
			
			cout << "Local fr: " << ix << ", " << iy << ": " << box->Lwin << endl;
			
			ising->draw_xy_L(ix,iy,resolution, box->Lwin / WIDTH, 0);
			
			sprintf(convert, "convert -resize 1024x1024 frame0000.png squares%4.4d.jpg", t);
			printf("%s\n\n", convert);
			system(convert);
		}

	}
	
	
		
	return 0;
}

void box_class::start_job()
{
	// Here I'll control the whole job.
	

	int x,y,ci,i,j;
	int current;
	int c0x,c0y, Lcell;
	// Apply some contrast
	double xp,yp, op, a;
	int radius;
	int row, col;

	// A little check
	if ((startx > L) || (starty>L) || (Lwin>L))		
	{
		cout << "Something is bigger than L." << endl;
		exit(1);
	}
	
	op = resolution / (Lwin);
	
	if (op<1)
	{
		radius=0;
		resolution = (int) Lwin;
	}
	else
	{
		radius = lround(op);
		
		
		// I want to change the resolution so the radius is an exact number of squares
			resolution = lround(Lwin * radius);
		cout << resolution << " " << op << " " << radius << endl;
		// recalculate
		
		Lwin = ((double) resolution / ((double) radius));
		op = resolution / (Lwin);
		radius = lround(op);
	}
	
	// Work out the cell range
	c0x = (int) (Wcells*startx/L);
	c0y = (int) (Wcells*starty/L);
	
	Lcell = (int) (Wcells*(startx+Lwin+1)/L) - c0x + 1;
	if ((c0x + Lcell > Wcells) || (c0y + Lcell > Wcells))
		Lcell=Wcells;
	
	cout << "Resolution: " << resolution << ". Square Length = " << radius << " pixels. Unrounded=" << op << ", Lwins=" << Lwin << endl;
	
	pngwriter png(resolution,resolution,1.0,"squares.png");
	
//	png.filledcircle(resolution/2,resolution/2,resolution,0.8,0.9,1.0);
//	png.flood_fill(resolution/2,resolution/2,0.8,0.9,1.0);
	
	
	for (y=c0y;y<c0y+Lcell;y++)
	{
		for (x=c0x;x<c0x+Lcell;x++)
		{
			ci = (y%Wcells)*Wcells + (x%Wcells);
			
			current = cell[ci];
			while (current < MAXNP)
			{
			
				xp = (particle[current].r[0]-0.5 - startx);
				xp = xp/Lwin;
				col = lround(xp*resolution);
				
				yp = (particle[current].r[1]-0.5 - starty);
				yp = yp/Lwin;
				row = lround(yp*resolution);
				
				if (radius==0)
					png.plot(col,row,0.0,0.0,0.0);
				else
					png.filledsquare(col, row, col+radius-1, row+radius-1, 0.0,0.0,0.0);
			
				current = particle[current].next;
			}

		}
	}
	
	png.close();
	
	
}

box_class::box_class()
{
	int i;
	int Nper_slab;
	
	// We have no particles at the start
	Np=0;

	// Create cells
	Wcells = (int) (L/RCUT + EPS);
	if (Wcells < 3) Wcells = 3; // Not likely.
	
	Ncells = (int) (pow((double) Wcells,DIMENSION) + EPS);
	cell_length = L/Wcells;
		
	cell = new int [Ncells];
	cout << "Created " << Ncells << " cells. Initialising..."; cout.flush();
	for (i=0;i<Ncells;i++)
		cell[i] = MAXNP;	// Sets all linked lists to NULL
	cout << "done." << endl;
	
	cout << "Created " << Ncells << " cells" << endl;
	
	cout << "Done. Np=" << Np << endl;
	
	vol = pow(L,DIMENSION);
	
	
}


void box_class::add_particle(int in_particle, int ic)
{
	int current;

	// Put the particle in the cell

	if (cell[ic] > MAXNP-1)
	{
		// Top of the list
		cell[ic] = in_particle;
	}
	else
	{
		// Find the end and add it
		current = cell[ic];
		while (particle[current].next < MAXNP)
			current = particle[current].next;
		particle[current].next = in_particle;
	}

	// New particle is at the end of the list
	particle[in_particle].next = MAXNP;
					   
	// Increase the particle count
	Np++;
	//#pragma omp atomic
	//Nptot++;	// Careful accessing this
}


	
int box_class::assign_cell(int inpart)
{
	int x,y,z,ci;
	

	//which cell?
	x = (int) (particle[inpart].r[0] / cell_length);
	y = (int) (particle[inpart].r[1] / cell_length);
	
#if DIMENSION==2
	ci = x + y*Wcells;
#elif DIMENSION==3
	z = (int) (particle[inpart].r[2] / cell_length);
	ci = x + y*Wcells + z*Wcells*Wcells;
#endif

	return(ci);
}


void box_class::load_positions(ising_class *ising, double cx, double cy)
{
	uint64_t i;
	int inpart, ic, Npt;
	int x,y;
	
	
	
	// Create the particle array
	Npt = (int) (L*L*0.9);
	particle = new particle_class [Npt];
	
	cout << "Npt=" << Npt << flush;
	
	for (inpart=0, i=0;i<N;i++)
	{
		if (ising->getspin(i))
		{
			x = (i%WIDTH);
			y = i/WIDTH;
			
			// The ising model is so huge we can only keep particles near the centre
			if ((x > cx - L/2) && (x < cx + L/2) && (y > cy - L/2) && (y < cy + L/2)) 
			{
				// annoyingly the coordinates are reversed (up<->down etc)
				particle[inpart].r[0]= L/2 + x + 0.5 - cx;
				particle[inpart].r[1]= L/2 + y + 0.5 - cy;
				ic = assign_cell(inpart);
				add_particle(inpart, ic);
				inpart ++;
			}
		}
	}
	
	cout << ", Np=" << Np << endl;
}


///////////////// ISING CLASS METHODS

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
	if (block_length==0)
		block_length=1;
	// We might want to round Wrn up
	Wrn = lround(WIDTH * Lfac / block_length);
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
		pngwriter png(Wrn,Wrn,0.0,filename);
		
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
				av_spin = 0.5 - tanh(steep * (av_spin - 0.5)) / tanh(0.5 * steep) / 2.0;	// tanh
				
				png.plot(col, row, av_spin, av_spin, av_spin);
				
			}
		}
		
		png.close();
	}
	else
	{
		cout << "Haven't built PBC into this yet, choose different x,y" << endl;
		cout << "x=" << x << " y=" << y << " block_length=" << block_length << " Wrn=" << Wrn << endl;
	}
	
}

ising_class::ising_class()
{
	
	spin = new unsigned char [N/8];
	
}

inline bool ising_class::getspin(long index)
{
	int shift = (index & 7);
	return (!((spin[(index >> 3)] & (1 << shift)) >> shift));
}

void ising_class::load_config(const char * filename)
{
	ifstream lattice;
	
	lattice.open(filename,ifstream::binary);
	lattice.read((char *) spin,N/8);
	lattice.close();
}
