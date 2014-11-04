#ifndef LATTICEBIP_H
#define LATTICEBIP_H

#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "site.h"
#include "lattice.h"
#include "cluster.h"

class lattice {
    // basically a collection of sites
    protected:
    // site *sites;
    std::vector<bool> siteVals;
    std::vector<int> bondvals;
    //
    std::vector<site> sites;
    std::vector<bond> bonds;
    //    
    std::map<int, int> boundaries;	// this will be <z*, index> (pointing to bonds).
    //int *clusters;
    int latsize;
    int width;
    int z0, randSeed;
    //
    // for now, put lattice type vars here. maybe we'll make a separate object later.
    //int nNeighbors;
    int geometry;
    //int* nDeltas;	// deltas for neighbors (define geometry)

    public:
    //
    lattice();
    lattice(int latlen0);
    lattice(int latlen, int randseed);
    // lattice(int latlen, int randseed, int* Zs);
    lattice(int latlen, int randseed, int latgeom);
    lattice(int latlen, int randseed, int latgeom, int dlamode);
    ~lattice();
    //
    void init(int latlen, int randseed, int latgeom);
    void init(int latlen, int randseed, int latgeom, int dlamode);
    //site* getsite(int index);
    bool getsiteval(int index);
    bool * getsite(int index);
    int getlatsize();
    int getgeometry();
    int getsite0();
    int getwidth();
    int getrandSeed();
    int randSeedPlus();
    int randSeedPlus(int);
    int getclustsize();
    int getdlamode();
    //
    //void setz(int, int);
    void setsite(int, bool);
    void setgeometry(int);
    void setz0(int);
    void initializeLattice();
    void initializeLattice(int rseed);
    void setrandSeed();
    void setrandSeed(int);
    void setdlamode(int);
    //
    void writegrid(std::string);
    //
    int dlamode;
    std::vector< std::vector<float> > getLatticeMap();		// <x,y,z, val>, noting that we need to change the current z-value convention.
    std::vector< std::vector<float> > getClusterMap();		// a more compact form of getLatticeMap()
    //
    // i-percolation specific (eventually, these should be ported to a subclass)
    cluster ipcluster;
    int zClust;		// default z value for cluster members.
    std::map<int, std::string> geometrylist;
    std::map<int, int> neighborcount;		// number. neighbors for geometries.
    //
    // get a full set of neighbors for one or many elements.
    std::map<int, int> getNeighbors(int k);
    std::map<int, int> getNeighbors(std::vector<int>&);
    std::map<int, int> getNeighbors(cluster);
    std::map<int, int> getNeighbors();	// assume ipcluster
    //
    std::map<int, int>* getboundaries();
	 //
    void initializeCluster(std::vector<int>);
    void initializeCluster(int csite);
    void initializeCluster();
    //
    int ipstep();
    int ipstep(int);	// int method: aggregation type (0=normal, 1=DLA, ...)
    int nOccNeighbors(int);	// number of occupied neighbors around index k
    int nOccNeighbors(int, int);
    int nOcN(int, int);
    int nOcN(int);				// shorthand for nOccNeighbors()
	 int makePercolatingCluster();
    int makePercolatingCluster(int);
    int makePercolatingCluster(int, int);	
    //int makePercolatingCluster(int, int, bool);	// aggregation type, percolation type (1,2,3,4 sides, etc.), refreshRand
	//
	};
#endif // LATTICE_H

