#include "sitebool.h"
//
sitebool::sitebool() {
	 bool * z;
	 *z=false;
	 init(z, std::vector<int>());
	}

sitebool::sitebool(bool * z) {
		  init(z, std::vector<int>());
	}

sitebool::sitebool(bool * z, int n) {
		  init(z, std::vector<int>(n));
	}
	
sitebool::sitebool(bool * z, std::vector<int> neighbs) {
	 init(z, neighbs);
	}

//void init(int*, int, std::vector<int>);
void sitebool::init(bool* Z, std::vector<int> neighbs) {
	this->occ=Z;
	this->neighbors=neighbs;
	}


bool * sitebool::getocc() {
	return this->occ;
	}
	
bool sitebool::isocc() {
	// return the pointer to the z-value.
	return *(this->occ);
	}

std::vector<int> sitebool::getNeighbors() {
	// return indices of bonds
	return this->neighbors;
	}

int sitebool::getNeighbor(int i) {
	if (int(this->neighbors.size())>i) {
		  return -1;
    }
	return this->neighbors[i];
	}

int sitebool::getnNeighbors() {
    return this->neighbors.size();
	}

int sitebool::popNeighbor(int i) {
	int rval=this->neighbors[i];
	this->neighbors.erase(this->neighbors.begin() + i);
	return rval;
	}

void sitebool::setNeighbors(std::vector<int> neighbs) {
	this->neighbors=neighbs;
	}
	
void sitebool::setNeighbor(int i0, int val) {
	this->neighbors[i0]=val;
	}

void sitebool::setocc() {
	this->setocc(true);
	}
void sitebool::setunocc() {
	this->setocc(false);
	}
void sitebool::setocc(bool b) {
	*(this->occ)=b;
	}
void sitebool::setocc(bool *b) {
	this->occ=b;
	}

