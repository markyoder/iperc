//#include <unistd.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>

#include "site.h"
std::vector<int> myvec();
int main();
std::vector<int> nothervec(std::vector<int>&);


std::vector<int> myvec()
	{
	int myints[]={1,2,3,5};
	std::vector<int> rvec(myints, myints+4);
	//
	return rvec;
	}

std::vector<int> nothervec(std::vector<int>& invec) {
	
	//invec = std::vector<int>(10,4);
	printf("from nothervec: %d\n", invec[0]);
	return std::vector<int>(invec.begin(), invec.end());
	}

int main()
{
	printf("this is my C++ test program.\n");
	//
	int L, Ntotal;
	L = 128;
	Ntotal = L*L;
	//
	int thisindex=L+L/2;
	int * ary1 = new int[5];
	int dneighbors[4] = {-1, 1, -L, L};
	std::vector<int> thisneighbor (dneighbors, dneighbors+4);
	//
	printf("vector bit[0]: %d, %d\n", thisneighbor[0], thisneighbor[1]);
	
	site s1=site();
	//
	//printf("site1: %d, %d\n", s1.getindex(), s1.getz());
	std::vector<int> newvec = myvec();
	std::vector<int> newvec2(newvec);
	std::vector<int> newvec3;
	std::vector<int> newvec4(nothervec(newvec));
	//
	printf("newvec4: %d\n", newvec4[0]);
	
	std::vector<int> vec5(5);
	printf("vec5.size(): %d\n", vec5.size());
	
	//newvec3.assign(newvec.begin(), newvec.end());
	//newvec3.assign(newvec);
	
	//newvec=myvec();
	printf("new vector: %d, %d\n", newvec[0], newvec[1]);
	printf("new vector: %d, %d\n", newvec2[0], newvec2[1]);
	
	//printf("new vector3: %d, %d\n", newvec3[2], newvec3[3]);
	//
	// modulus tests:
	//printf("mods: 5%100= %d, 105%100= %d, -5%100= %d", (5%100), (105%100), (-5%100));
	unsigned int x0=100;
	std::cout << "mods: 5%100=" << 5%100 << ", 105%100= " << 105%100 << ", -5%100=" << -5%100 <<"\n";
	std::cout << "mods: 5%100=" << 5%x0 << ", 105%100= " << 105%x0 << ", -5%100=" << -5%x0 <<"\n";
	std::cout << "mods: squares=" << 144%12 << ", 105%100= " << 36%6 << ", -5%100=" << 37%6 <<"\n";
	//
	// casting:
	int a=144;
	float b;
	b = float (a);
	
	printf("a: %d, b: %f, %f, %f\n", a, b, pow(b,.5), pow(a, .5));
	
	
	
}


