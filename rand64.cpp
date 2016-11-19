#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;
#include "/home/da246/MersenneTwister.h"

int main()
{
	long x,i,com;

	MTRand mt;
com=4294967297;
	for (i=0;i<1000000;i++)
	{
		x = mt.randInt(4294967295) ;
		x += mt.randInt(4096 - 1) * 4294967296;
	//	cout << x << endl;
		if (x>com)
			cout << x << endl;	
	}
	
	return 0;
}

