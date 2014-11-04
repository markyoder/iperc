#ifndef CLUSTER_H
#define CLUSTER_H

#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

class cluster
{
protected:
	//std::vector<int> elements;
	//std::vector<int> perimeter;
	//int size;	// this may be used to track size, in the event that elements is pre-expanded and we can't get an accurate size from
					// the vector object. this will also be used as a checksum; size=nSteps=elements.sizeOf() (or something like that).

public:
    cluster();
	 cluster(int);
	 cluster(std::vector<int>);	//initialize with multiple initial sites
	// cluster(std::vector<int>&);
	 void init(std::vector<int>);
	 ~cluster();
	 //
	 std::vector<int> elements;
	 //std::vector<int> perimeter;	// probably best to make these public, rather than write aliases for all the vector member functs.
	 										// but in the end, we might not even use this perimeter, only the "boundaries" feature from lattice.
	 //std::vector<std::vector<int>> elements2;		// this will be like [ [sequence index, lattice index, zval], [], ...]	maybe.
	 																// nominally, we can get away with using the vector order, but it's not a totally
	 																// stable approach.
	 std::vector<int> displacedZs;			// displaced z-values (z-values of elements). (note: elements, displacedZs basically constitute
	 std::vector<std::vector<int> > edges;	// an <int, int> vector. edges = <val, from, to>. note edge index can be calculated
	 													//  from the "from" index.
	 													// for the time being, we'll keep these in chronological order.
	 //
	 // members:
	 //int size;
	 
	 // member functions:
	 void addelement(int);
	 void addperimeter(int);
	 void setelement(int, int);
	 void setelements(std::vector<int>&);		// set the whole enchilada
	 void setperimeters(std::vector<int>&);	// set the whole enchilada
	 void setperimeter(int, int);
	 void eraseelement(int, int);
	 void eraseperimeter(int, int);
	 void addDisplacedZ(int);
	 void addEdge(int, int);		// ifrom, ito
	 void addEdge(int, int, int);	// val, ifrom, ito
	 void addEdge(std::vector<int>&);
	 
	 //void addedge(int, int);
	 //void addedges(std::multimap<int, int>);


	 int getsize();
	 std::vector<int> getelements();
	 std::vector<int> getperimeters();	// return the whole perimeter/elements
	 std::vector<int> getDisplacedZs();
	 //std::multimap<int, int> getedges();
	 std::vector< std::vector<int> > getedges();	// use a plain vector structure, this way we can keep chronological order.
	 std::vector<int> getedge(int);
	 int getelement(int);
	 int getperimeter(int);					// get one element (but just a copy of it -- don't pop()
	 int popelement(int);
	 int popperimeter(int);
	 int getDisplacedZ(int);

};

#endif // CLUSTER_H


