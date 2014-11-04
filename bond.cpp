#include "bond.h"

//
bond::bond() {
	int * z = 0;
	init(0,0, z,0);
	}

bond::bond(int to) {
	int * z = 0;
	init(to, 0, z, 0);
	}

bond::bond(int to, int from) {
	int * z = 0;
	init(to, from, z, 0);
	}

bond::bond(int to, int from, int * z) {
	init(to, from, z, 0);
	}
bond::bond(int to, int from, int zval) {
	int * z;
	*z = zval;
	init(to, from, z, 0);
	}
bond::bond(int to, int from, int zval, bool isocc) {
	int * z;
	*z = zval;
	init(to, from, z, isocc);
	}
bond::bond(int to, int from, int * z, bool isocc) {
	init(to, from, z, isocc);
	}
	
//void init(int*, int, std::vector<int>);
//init(int, int, int*, bool);
void bond::init(int to, int from, int* val, bool isocc) 	{
	this->to=to;
	this->from=from;
	this->val = val;
	this->occ=isocc;
	//
	if (this->val==NULL) {
		*this->val=0;
		}
	if (this->occ==NULL) {
		this->occ=0;
		}
	}

int bond::getto() {
	return this->to;
	}

int bond::getfrom() {
	return this->from;
	}

int bond::getval() {
	return *(this->val);
	}

bool bond::isocc() {
	return this->occ;
	}


void bond::setto(int x) {
	this->to=x;
	}

void bond::setfrom(int x) {
	this->from=x;
	}

void bond::setval(int x) {
	*(this->val)=x;
	}

void bond::setval(int * x) {
	this->val = x;
	}

void bond::setocc() {
	this->setocc(true);
	}
void bond::setunocc() {
	this->setocc(false);
	}
void bond::setocc(bool x) {
	this->occ=x;
	}
void bond::setocc(bool * x) {
	&(this->occ=x);
	}


