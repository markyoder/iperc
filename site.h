#ifndef SITE_H
#define SITE_H

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

class site {
	// a lattice site.
	protected:
		//int index;	// effectively its position
		//std::vector<int> neighbors;
		//int* lattice;
		int* val;		// this will either be an independent int* pointer or will point to a lattice (array).
		std::vector<int> neighbors;
		//int nNeighbors;
	public:		
		// constructors:
		site();
		site(int*);
		site(int*, int);
		site(int*, int, int);
		~site();
		void init(int*, int, int); // z, index, num-neighbors.

		// functional bits:
		int  getindex();
		int* getz();	// return pointer val to z-value
		int getval();	// return z value (at val pointer)
		//
		//std::vector<int> getNeighbors();	// returns a vector of neighbor addresses.
		std::vector<int> getNeighbors();
		int getNeighbor(int);
		int getnNeighbors();
		//
		//void setNeighbors(std::vector<int>);
		void setNeighbor(int, int);
		void setz(int);
		void setNeighbors(std::vector<int>);
		void setNeighbors(int*, int);	//array*, length.
		//
	

	};


#endif
