#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <string>

#include "site.h"
#include "lattice.h"
// note: compile with site.o/cpp, lattice.o/cpp, etc.

int main() {
	int width=256;
	int randseed=time(NULL);
	int latgeom=0;
	int dlamode=1;
	int percmode=1;
   int drs=0;
   std::cout << "zero...\n";
   lattice thislat;
   std::cout << "first set..\n";
	//lattice::lattice lat = lattice::lattice(width, randseed, latgeom);
	//
	//
   //lattice thislat = lattice(width, randseed+drs, latgeom);
  	
  	thislat.init(width, randseed+drs, latgeom, 1);
	thislat.initializeCluster(width*width/2 + width/2);	// in the middle.
	thislat.makePercolatingCluster(dlamode, percmode);
	thislat.writegrid(std::string("gridout-square-dla.dat"));
	drs++;
   //
	std::cout << "second set..\n";
	thislat = lattice(width, randseed+drs, latgeom);
	thislat.initializeCluster(width*width/2 + width/2);	// in the middle.
	thislat.makePercolatingCluster(0, percmode);
	thislat.writegrid(std::string("gridout-square.dat"));
	drs++;
	//
	// one more...
	thislat = lattice(width, randseed+drs, 1);
	thislat.initializeCluster(width*width/2 + width/2);	// in the middle.
	thislat.makePercolatingCluster(dlamode, percmode);		// dla, perc-1 
	thislat.writegrid(std::string("gridout-diam-dla.dat"));
	drs++;
	
	for (int i=0; i<1; i++) {
	// one more...
	lattice newlat = lattice(width, randseed+drs, 2);
	newlat.initializeCluster(width*width/2 + width/2);	// in the middle.
	newlat.makePercolatingCluster(dlamode, percmode);		// dla, perc-1 
	std::cout << "writing file " << i << "\n";
	newlat.writegrid(std::string("data/gridout-tri-dla-%d.dat", i));
	drs++;
	//delete &newlat;
	}
	
	/*
	std::vector< std::vector<float> > L = thislat.getLatticeMap();
	for (int i=0; i<5; i++) {
		for (int j=0; j<4; j++) {
			std::cout << "i,j: " << i << ", " << j << ", " << L[i].size() << "\n";
			}
		}
	*/
	
	return 0;
	};
