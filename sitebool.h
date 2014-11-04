#ifndef SITEBOOL_H
#define SITEBOOL_H

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

// site object designed to be used in a site-bond hybrid percolation lattice.

class sitebool {
	// a lattice site.
	protected:
		// int * val // depricating val -> bool(isocc).
		// int index // (index will be implied by order in array to reduce memory footprint).
		//
		bool * occ;		// is occupied?
		std::vector<int> neighbors;	// each site has neighbors. the lattice will interpret bonds from <site, neighbor> pairs.
		//
	public:		
		// constructors:
		sitebool();
		sitebool(bool*);
		sitebool(bool*, int nneighbors);
		sitebool(bool*, std::vector<int> neighbs);
		//~site();
		void init(bool*, std::vector<int>); // occ, neighbors

		// functional bits:
		//int  getindex();
		bool * getocc();	// return pointer val to occ-state
		bool isocc();	// return z value (at val pointer)
		//
		//std::vector<int> getNeighbors();	// returns a vector of neighbor addresses.
		
		std::vector<int> getNeighbors();
		int getNeighbor(int);
		//std::pair<int, int> getNeighborBond(int);
      int getnNeighbors();
      int popNeighbor(int i);	// remove ith neighbor
		//
		void setNeighbors(std::vector<int>);
		void setNeighbor(int, int);
		void setocc();
		void setunocc();
		void setocc(bool);
		void setocc(bool*);	// set occ address
		
		void setNeighbors(int*, int);	//array*, length.
		//
	

	};


#endif
