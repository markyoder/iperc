#include "site.h"

//

site::site()
{
	 init(new int[0], 0, 0);
}

site::site(int* Z)
{
		  init(Z, 0, 0);
}

site::site(int* Z, int index)
{
		  init(Z, index,0);
}

site::site(int* Z, int index, int nNeighbors)
{
	 init(Z, index, nNeighbors);
}

//void init(int*, int, std::vector<int>);
void site::init(int* Z, int index0, int NN)
	{
		// note: index variable is dead. for now, just leave it in place in case we want to revive it.
		//
		this->val = Z;
		//this->index=index0;
		//neighbors.assign(neighbors0.begin(), neighbors0.end());
		//neighbors = neighbors0;   // or some syntax like this.
		//
		//this->nNeighbors = NN;
		//std::cout<<"initiate neighbors in site..." << index0 << " with " << nNeighbors << " neighbors \n";
		//neighbors = new int[Nneighbors];
		this->neighbors = std::vector<int>(NN);
		// and let's leave setting of neighbors to a separate call for now?
		//std::cout<<"initiated site..." << index0 << "\n";
	}

site::~site() {
	//delete &(this->neighbors);
	//delete &nNeighbors;
	//delete &index;
	//delete val;
	}

int site::getval() {
	return *val;	// return the value at the pointer val.
}

int* site::getz() {
	// return the pointer to the z-value.
	return val;
}

//int site::getindex() {
//    return index;
//}

std::vector<int> site::getNeighbors() {
	return neighbors;
}

int site::getNeighbor(int nindex) {
	 if (nindex>this->getnNeighbors()) {
		  return -1;
    }
    return neighbors[nindex];
}

int site::getnNeighbors() {
   // return this->nNeighbors;
   return int(this->neighbors.size());
}

void site::setz(int thisz){
	 //lattice[index] = thisz;
	 *val = thisz;
}

void site::setNeighbors(std::vector<int> invec) {
	neighbors = std::vector<int>(invec);
	};

void site::setNeighbors(int * inlist, int listlen) {
	neighbors = std::vector<int>(inlist, inlist+listlen);
	}
	

