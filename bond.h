#ifndef BOND_H
#define BOND_H

#include <unistd.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include <time.h>
#include <string>
#include <vector>
#include <sstream>

#include <cmath>
#include <cassert>
#include <iomanip>

//int z0=10;

class bond {
		
	protected:
		int to, from;
		int * val;
		bool occ;
	
	public:
		// functions that return stuff:
		bond();
		bond(int);
		bond(int, int);
		bond(int, int, int*);
		bond(int, int, int);
		bond(int, int, int*, bool);
		bond(int, int, int, bool);
		//
		void init(int, int, int*, bool);
		//
		int getto();
		int getfrom();
		int getval();
		
		//int * getval();
		bool isocc();				// and this is, more or less, an optional pram. we can also use a range of values to determine this.
		std::pair<int, int> getfromto();
		//
		// functions that do stuff:
		void setto(int);
		void setfrom(int);
		void setval(int);
		void setval(int * x);
		void setocc();			// default to "True": setocc(1)
		void setocc(bool);
		void setocc(bool*);
		void setunocc();			// setocc(0)
		
	};


#endif
