#include "cluster.h"


cluster::cluster() {
	// initialize the 0 index.
        //init(std::vector<int>(1,0));
    init(std::vector<int>());
	}

cluster::cluster(int index) {
	init(std::vector<int>(1, index));
	}

cluster::cluster(std::vector<int> invec) {
	//
	init(std::vector<int>(invec.begin(), invec.end()));
	}

//cluster::cluster(std::vector<int>& invec) {
//	//
//	init(*invec);
//}

void cluster::init(std::vector<int> invec) {
	// note on vectors: to pass an existing vector, use the vector<>& input type (aka, pass the address of a vector):
	// int = bar(vector<>& invec) { blah, blah}; int foo = funct(existing_vector)
	// to pass a created vector, use vector<> input type: aka:
	// int bar(vector<> invec) { blah, blah}; foo = bar(vector<>()); (aka, send a copy).
	//
	//size = invec.size();	// be sure this does not include empty values, which should be accruate. note .capacity() indicates allocated memory.
	this->elements = std::vector<int>(invec.begin(), invec.end());
	this->displacedZs = std::vector<int>();
	//
	// note: cluster() is not a terribly smart class. it does not know its lattice geometry, so it will have to be told its perimeter.
	// it still makes a good place, however, to store the perimter (though it is equally a member of the lattice class as the cluster class).
	}

cluster::~cluster() {
/*
	delete &elements;
	delete &displacedZs;
	//delete &perimeter;
	delete &edges;
	*/
	}

void cluster::addelement(int val) {
	elements.push_back(val);
}

//void cluster::addperimeter(int val) {
//	perimeter.push_back(val);
//}

void cluster::setelement(int i, int val) {
	elements[i]=val;
}

//void cluster::setperimeter(int i, int val) {
//	perimeter[i]=val;
//}

void cluster::addDisplacedZ(int i) {
	this->displacedZs.push_back(i);
	}

void cluster::addEdge(int ifrom, int ito) {
	// note: we're using a default value=0
	this->addEdge(0, ifrom, ito);
	}
	
void cluster::addEdge(int val, int ifrom, int ito) {
	int tmpints[] = {val, ifrom, ito};
	std::vector<int> tmpvec(tmpints, tmpints+3);
	//this->addEdge(std::vector<int>(tmpints, tmpints+3));
	this->addEdge(tmpvec);
	}
	
void cluster::addEdge(std::vector<int>& invec) {
	while (invec.size()<3) {
		invec.insert(invec.begin(),0);
		}
	this->edges.push_back(invec);
	}

void cluster::eraseelement(int i, int j) {
	elements.erase(elements.begin()+i, elements.begin()+j);
}

//void cluster::eraseperimeter(int i, int j) {
//	perimeter.erase(elements.begin()+i, elements.begin()+j);
//}

int cluster::getsize() {
	return this->elements.size();
	}

std::vector<int> cluster::getelements() {
	return elements;
}

//std::vector<int> cluster::getperimeters() {
//	return perimeter;
//}

std::vector<int> cluster::getDisplacedZs() {
	return this->displacedZs;
}

int cluster::getelement(int i) {
	return elements[i];
}

std::vector<std::vector<int> > cluster::getedges() {
	return this->edges;
	};
std::vector<int> cluster::getedge(int i) {
	return this->edges[i];
	}

//int cluster::getperimeter(int i) {
//	return perimeter[i];
//}

int cluster::popelement(int i) {
	int relem = elements[i];
	elements.erase(elements.begin()+i);
	return relem;
}
/*
int cluster::popperimeter(int i) {
	int relem = perimeter[i];
	perimeter.erase(elements.begin()+i);
	return relem;
}
*/

int cluster::getDisplacedZ(int i) {
	return this->displacedZs[i];
	}


