#ifndef LATTICE_H
#define LATTICE_H

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
#include "bond.h"

class lattice {
    // 
    protected:
    // site *sites;
    std::vector<site> sites;
    std::vector<bond> bonds;
    //
    std::vector<int> siteVals;
    std::vector<int> bondVals;
    //
    //int* siteVals;
    std::map<int, int> boundaries;	// this will be <z*, index>	// and now pointing to bonds, not sites.
    											// boundary sites can be pulled from bond.from values.
    //
    int latsize, width, z0, randSeed;
    //
    // for now, put lattice type vars here. maybe we'll make a separate object later.
    int nNeighbors;
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
    //~lattice();
    //
    void init(int latlen, int randseed, int latgeom);
    void init(int latlen, int randseed, int latgeom, int dlamode);
    //site* getsite(int index);
    int getsiteval(int index);
    int* getsite(int index);
    int getbondval(int index);
    int * getbond(int index);
    bond * getbondObj(int);
    int getlatsize();
    int getbondsize();
    int getnNN();	// num nearest neighbors (for lattice geometry)
    int getgeometry();
    int getsite0();
    int getwidth();
    int getrandSeed();
    int randSeedPlus();
    int randSeedPlus(int);
    int getclustsize();
    int getdlamode();
    int getBondIndex(int, int); 		// siteindex, neighborindex
    //
    //void setz(int, int);
    void setval(int, int);		// should eventually be chaged to setsiteval()
    void setbondval(int, int);
    void setgeometry(int);
    void setz0(int);
    void initializeLattice();
    void initializeLattice(int rseed);
    void initializeSites();
    void initializeSites(int);
    void initializeBonds();
    void initializeBonds(int);
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
    std::map<int, int> geomNN;	// geom-num neighbors
    //
    // get a full set of neighbors (sites) for one or many elements.
    std::map<int, int> getNeighbors(int k);
    std::map<int, int> getNeighbors(std::vector<int>&);
    std::map<int, int> getNeighbors(cluster);
    std::map<int, int> getNeighbors();	// assume ipcluster; neighbor sites?
    //
    // get a full set of neighbors (Edge) for one or many elements.
	 std::map<int, int> getNeighborBonds();	// assume ipcluster
    std::map<int, int> getNeighborBonds(int k);
    std::map<int, int> getNeighborBonds(std::vector<int>&);
    std::map<int, int> getNeighborBonds(cluster);
    std::vector< std::vector<int> > getTips();		// get ends of branches; sites with 1 neighbor.
    //
    std::map<int, int>* getboundaries();
	 //
    void initializeCluster(std::vector<int>);
    void initializeCluster(int csite);
    void initializeCluster();
    //
    int ipstep();
    int ipstep(int);	// int method: aggregation type (0=normal, 1=DLA, ...)
    int ipstepSite();
    int ipstepSite(int);	// "legacy" ipstep for pure site percolation.
	 int nOccBonds(int);
	 int nOccBonds(int, int);	// number of occupied edges...
    int nOccNeighbors(int);	// number of occupied neighbors around index k
    int nOccNeighbors(int, int);
    int nOcN(int, int);
    int nOcN(int);				// shorthand for nOccNeighbors()
	 int makePercolatingCluster();
    int makePercolatingCluster(int);
    int makePercolatingCluster(int, int);	
    //int makePercolatingCluster(int, int, bool);	// aggregation type, percolation type (1,2,3,4 sides, etc.), refreshRand

    //
    //std::map<int, int> boundaries;

};
#endif // LATTICE_H

