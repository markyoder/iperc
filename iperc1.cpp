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
	int thislatlen=256;
	int randseed=time(NULL);
	int latgeom=0;
	//lattice::lattice lat = lattice::lattice(thislatlen, randseed, latgeom);
	//
	lattice thislat(thislatlen, randseed, latgeom);
	
	//printf(" some lattice bits: %d, %d, %d, %d", 5,6, 7, 8);
	printf(" some lattice bits : %d, %d, %d, %d\n", thislat.getlatsize(), thislat.getgeometry(), thislat.getsiteval(0), thislat.getsiteval(100));
	printf(" some lattice bits*: %d, %d, %d, %d\n", thislat.getlatsize(), thislat.getgeometry(), *thislat.getsite(0), *thislat.getsite(100));
	//, thislat.getsite(0), thislat.getsite(100));
	thislat.setval(72, 45);
	printf(" some lattice bits: %d, %d, %d, %d\n", thislat.getlatsize(), thislat.getgeometry(), *thislat.getsite(72), thislat.getwidth());
	//
        thislat.initializeCluster(std::vector<int>(1,2024));
        thislat.ipcluster.elements.push_back(2025);
        thislat.ipcluster.elements.push_back(2026);
        //
        thislat.initializeCluster(thislat.ipcluster.getelements());

        printf("initialized cluster: %d\n", thislat.ipcluster.getelement(0));
	thislat.writegrid(std::string("gridout0.dat"));
	//
	std::map<int, int> myset = thislat.getNeighbors(thislat.ipcluster.elements);
	//printf("Neighbors, as returned by lattice::getNeighbors().\n");
	/*
	for (std::map<int, int>::iterator it=myset.begin(); it!=myset.end(); ++it) {
		//printf("neighbors: %d, %d\n", x.first, x.second);
		std::cout << "neighbors: " << it->first << ": " << it->second << '\n';
		}
	printf("\n");
	*/
	/*
	std::map<int, int>* b1 = thislat.getboundaries();
	for (std::map<int, int>::iterator it=b1->begin(); it!=b1->end(); ++it) {
		//printf("neighbors: %d, %d\n", x.first, x.second);
		std::cout << "neighbors2: " << it->first << ": " << it->second << '\n';
		}
	printf("\n");
	*/
	//
	thislat.ipstep();
	//
	thislat.writegrid(std::string("gridout2.dat"));
	myset = thislat.getNeighbors(thislat.ipcluster.elements);
	//printf("Neighbors, as returned by lattice::getNeighbors().\n");
	for (std::map<int, int>::iterator it=myset.begin(); it!=myset.end(); ++it) {
		//printf("neighbors: %d, %d\n", x.first, x.second);
		std::cout << "neighbors: " << it->first << ": " << it->second << '\n';
		}
	printf("\n");
	
	thislat.ipstep();
	//
	thislat.writegrid(std::string("gridout2.dat"));
	myset = thislat.getNeighbors(thislat.ipcluster.elements);
	//printf("Neighbors, as returned by lattice::getNeighbors().\n");
	for (std::map<int, int>::iterator it=myset.begin(); it!=myset.end(); ++it) {
		//printf("neighbors: %d, %d\n", x.first, x.second);
		std::cout << "neighbors: " << it->first << ": " << it->second << '\n';
		}
	printf("\n");
	
	int maxI = thislatlen*25;
	std::cout<<"maxI: " << maxI << "\n";
	
	for (int i=0; i<(maxI); i++) {
		thislat.ipstep(1);
		}
	//
	myset = thislat.getNeighbors(thislat.ipcluster.elements);
	//printf("Neighbors, as returned by lattice::getNeighbors().\n");
	for (std::map<int, int>::iterator it=myset.begin(); it!=myset.end(); ++it) {
		//printf("neighbors: %d, %d\n", x.first, x.second);
		std::cout << "neighborses: " << it->first << ": " << it->second << " : " << thislat.nOcN(it->second) << "\n";
		}
	printf("\n");
	thislat.writegrid(std::string("gridout2.dat"));
	
	};
