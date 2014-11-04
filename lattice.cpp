#include "lattice.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>


lattice::lattice() {
	init(1, 0, 0, 1);
}

lattice::lattice(int width0) {
	 init(width0, time(NULL), 0, 1);
}

lattice::lattice(int width0, int randseed ) {
	 init(width0, randseed,  0, 1);
}

lattice::lattice(int width0, int randseed, int latgeom) {
	 init(width0, randseed, latgeom, 1);
}

lattice::lattice(int width0, int randseed, int latgeom, int dlamode) {
	 init(width0, randseed, latgeom, dlamode);
}

void lattice::init(int width0, int randseed, int latgeom){
	init(width0, randseed, latgeom, 1);
	}

//void lattice::init(int latlen, int randseed, int* Zs, int latgeom){
void lattice::init(int width0, int randseed, int latgeom, int dlamode){
	//
	if (randseed==0) {
		randseed = time(NULL);
		}
	//
	this->dlamode=dlamode;
	this->latsize=width0*width0;
	//if (latgeom==3) this->latsize=width0*width0*width0; // cubic...
	if (latgeom==3) this->latsize*=width0;
	//	
	this->width=width0;
	this->randSeed=randseed;
	//srand(randseed);
	this->geometry = latgeom;		
	// std::map<int, std::string> geometrylist;
	this->geometrylist[0]="square";
	this->geometrylist[1]="diamond";
	this->geometrylist[2]="triangular";
	this->geometrylist[3]="cubic";
	//
	this->geomNN[0]=4;
	this->geomNN[1]=4;
	this->geomNN[2]=6;
	this->geomNN[3]=6;
	//
	this->dlamode=0;	// default... (for this project, we might make default=1, dla).
	this->z0=10;
	this->zClust = 1;	// value for cluster members.
	//
	//this->siteVals = new int[latsize];
	//this->sites = new site[latsize];
	this->siteVals = std::vector<int>(latsize);
	this->sites = std::vector<site>(latsize);
	this->boundaries = std::map<int, int>();	
	//
	int bsize=this->latsize*(this->geomNN[latgeom]);
	this->bondVals =  std::vector<int>(bsize);
	//
	//std::cout<<"bondvals initiated (" << bsize << ", " << (this->geomNN[latgeom]) << ")\n";
	this->bonds = std::vector<bond>();
	//
	this->initializeLattice(randseed);
	//
	// define geometry:
	this->setgeometry(latgeom);
	this->ipcluster = cluster();	// a cluster object.
	//this->initializeCluster();
	//
	}

/*
lattice::~lattice() {
	
	delete &sites;
	delete &siteVals;
   delete &boundaries;	// this will be <z*, index>
   delete &latsize;
   delete &width;
   delete &z0;
   delete &randSeed;
   delete &ipcluster;
    //
    // for now, put lattice type vars here. maybe we'll make a separate object later.
    delete &nNeighbors;
    delete &geometry;
	}
*/

void lattice::initializeLattice() {
	int newrseed=this->randSeedPlus();
	initializeLattice(newrseed);
	}

void lattice::initializeLattice(int rseed) {
	this->initializeSites(rseed);
	rseed=this->randSeedPlus();
	this->initializeBonds(rseed);
	this->boundaries = std::map<int, int>();
	}

void lattice::initializeSites() {
	//
	int newrseed=this->randSeedPlus();
	initializeSites(newrseed);
	}
	
void lattice::initializeSites(int rseed) {
	//int dummyvar=0;
	// it is quite important that both the indices and values for each lattice site are unique.
	// rather than be at the mercy of rand(), let's use a randomly ordered sequence of numbers.
	// use: z0 < z < latsize+z0 ordered randomly.
	srand(rseed);
	std::map<int, int> tmpZs;
	std::map<int, int>::iterator it;
	for (int i=0; i<this->latsize;i++) {
		// this map will be ordered by the first column (key), a random number.
		tmpZs.insert(std::pair<int, int>(rand(), this->z0+i));
		}
	//
	//for (int i=0; i<latsize;i++) {
	int i=0;
	for (it=tmpZs.begin(); it!=tmpZs.end(); ++it) {
		//siteVals[i]=z0+rand()%(RAND_MAX-z0);
		this->siteVals[i] = (*it).second;
		//std::cout << "this zval: " << siteVals[i] << "\n";
		i++;
		}
	//
	//this->boundaries = std::map<int, int>();
	}

void lattice::initializeBonds() {
	//
	int newrseed=this->randSeedPlus();
	initializeBonds(newrseed);
	}

void lattice::initializeBonds(int rseed) {
	//int dummyvar=0;
	// it is quite important that both the indices and values for each lattice site are unique.
	// rather than be at the mercy of rand(), let's use a randomly ordered sequence of numbers.
	// use: z0 < z < latsize+z0 ordered randomly.
	srand(rseed);
	std::map<int, int> tmpZs;
	std::map<int, int>::iterator it;
	//for (int i=0; i<this->getbondsize();i++) {
	for (int i=0; i<((this->getnNN()))*(this->latsize); i++) {
		// this map will be ordered by the first column (key), a random number.
		tmpZs.insert(std::pair<int, int>(rand(), this->z0+i));
		}
	//
	//for (int i=0; i<latsize;i++) {
	int i=0;
	for (it=tmpZs.begin(); it!=tmpZs.end(); ++it) {
		//siteVals[i]=z0+rand()%(RAND_MAX-z0);
		this->bondVals[i] = it->second;
		
		//this->bondVals.push_back(it->second);
		//std::cout << "this zval: " << it->second << "\n";
		//std::cout<< "bondVal:" << bondVals[i] << "\n";
		i++;
		}
	//
	}
	
void lattice::setrandSeed() {
	this->randSeed=time(NULL);
	}

void lattice::setrandSeed(int x) {
	this->randSeed=x;
	}

void lattice::setgeometry(int geom) {
	//
	this->geometry=geom;
	this->nNeighbors=geomNN[geom];
	//
	int w=width;
	int w2=width*width;
	//
	switch (geom) {
	case 0:{
      std:: cout << "setting square-geometry" << "\n";
      // square:
		//nNeighbors = 4;
		
		//
		// neighbors:
		int tmpNeighbors[4] = {-1, 1, -w, w};
		//
		//std::cout << "begin setting sites and neighborses\n";
		for (int i=0; i<latsize;i++) {
			// site::site(int* Z, int index, int nNeighbors)
			//std::cout << "create site #" << i << "//" << latsize << "\n";
			//this->sites[i] = site(&siteVals[i], i, nNeighbors);
			this->sites[i].init(&siteVals[i], i, nNeighbors);
			
			//std::cout<<"site #" << i << "/" << latsize << " set.\n";
			//std::cout << "created site #" << i << "\n";
			//
			int* thisneighbors = new int[4];
			for(int j=0; j<4; j++) {
				thisneighbors[j] = tmpNeighbors[j]+i;
				}
			if (i%width==0) {
				// left boundary
				thisneighbors[0]+=w;
				};
			if (i%width==(width-1)) {
				// right
				thisneighbors[1]-=w;
				};
			if (i/width==0) {
				//thisneighbors[2]+=latsize;
				thisneighbors[2]+=w2;
				};
			if (i/width == (width-1)) {
				//thisneighbors[3]-=latsize;
				thisneighbors[3]-=w2;
				};
			
			//std::cout << "itting through lattice. " << i << "\n";
			sites[i].setNeighbors(thisneighbors, nNeighbors);
			};
      //std::cout << "set sites and neighbors\n";
      break;
      }
		case 1: {
		// diamond (NNN):
      std::cout << "setting diamond-geometry" << "\n";
		//nNeighbors = 4;
		//
		// neighbors:
		int tmpNeighbors[4] = {-width-1, -width+1, width-1, width+1};
		//
		for (int i=0; i<latsize;i++) {
			// site::site(int* Z, int index, int nNeighbors)
			sites[i] = site(&siteVals[i], i, nNeighbors);
			//
			int thisneighbors[4];
			for(int j=0; j<nNeighbors; j++) {
				thisneighbors[j] = tmpNeighbors[j]+i;
				}
			
			// nominally, these can be interpreted a couple of different ways for diamond geometry, but this should be fine.
			if (i%width==0) {
				// left boundary
				thisneighbors[0]+=width;
				thisneighbors[2]+=width;
				};
			if (i%width==(width-1)) {
				// right
				thisneighbors[1]-=width;
				thisneighbors[3]-=width;
				};
			if (i/width==0) {
				thisneighbors[1]+=latsize;
				thisneighbors[0]+=latsize;
				};
			if (i/width == (width-1)) {
				thisneighbors[3]-=latsize;
				thisneighbors[2]-=latsize;
				};
			
			//cout << "itting through lattice. " << i << "\n";
			sites[i].setNeighbors(thisneighbors, nNeighbors);
			};
         break;
			}
		case 2: {
			// triangular...:
			//nNeighbors = 6;
		   std::cout << "setting tri-geometry" << "\n";
			//
			// neighbors:
			int tmpNeighbors[6] = {width, -width, -1, 1, width+1, -width+1};	// or width-1, -width-1 (aka, -(width+1), -(1-width) for odd/even).
			int parity;	// next-nearest neighbor terms.
			//
			for (int i=0; i<latsize;i++) {
				parity = (-1 + 2*((i/width)%2));				// note the fancy way of writing odd/even
				tmpNeighbors[4] = parity*(width+1);			// even rows "reach back"; odd rows should be shifted right (positive)
				tmpNeighbors[5] = parity*(-width+1);						
				//
				sites[i] = site(&siteVals[i], i, nNeighbors);
				//
				int thisneighbors[6];
				for(int j=0; j<nNeighbors; j++) {
					// a crude approach to fixing the boundary conditions.
					thisneighbors[j] = tmpNeighbors[j]+i;
					if (thisneighbors[j]<0) {
						thisneighbors[j]+=latsize;
						}
					thisneighbors[j]=thisneighbors[j]%latsize;
					}
				sites[i].setNeighbors(thisneighbors, nNeighbors);
				};
		      break;
		      }
		case 3: {
			// cubic:
			//this->nNeighbors = 6;
			std::cout<< "setting cubic geometry.\n";
			//
			int w3=width*width*width;
			int tmpNeighbors[6] = {-1, 1, -w, w, -w2, w2};
			//
			//std::cout << "begin setting sites and neighborses\n";
			for (int i=0; i<latsize;i++) {
				// site::site(int* Z, int index, int nNeighbors)
				//std::cout << "create site #" << i << "//" << latsize << "\n";
				//this->sites[i] = site(&siteVals[i], i, nNeighbors);
				this->sites[i].init(&siteVals[i], i, nNeighbors);
			
				//std::cout<<"site #" << i << "/" << latsize << " set.\n";
				//std::cout << "created site #" << i << "\n";
				//
				int* thisneighbors = new int[6];
				for(unsigned int j=0; j<6; j++) {
					thisneighbors[j] = tmpNeighbors[j]+i;
					}
				if (i%width==0) {
					// left boundary
					thisneighbors[0]+=w;
					};
				if (i%width==(width-1)) {
					// right
					thisneighbors[1]-=w;
					};
				if (i/width==0) {
					thisneighbors[2]+=w2;
					};
				if (i/width == (width-1)) {
					thisneighbors[3]-=w2;
					};
				if (i/(width*width)==0) {
					thisneighbors[4]+=w3;
					}
				if (i/(width*width)==(width*width-1)) {
					thisneighbors[5]-=w3;
					}
				//std::cout << "itting through lattice. " << i << "\n";
				sites[i].setNeighbors(thisneighbors, nNeighbors);
				};
      //std::cout << "set sites and neighbors\n";
      break;
			}
         default: {
				std::cout << "else-etting square geom." << "/n";
				this->setgeometry(0);
				break;
				}
		};
	// now, set edges (this second looping is probably not the most efficient, but as i recall there is
	// some trick to getting the CASE statement to fire correctly inside a loop)
	int nnn=this->getnNN();
	int k;
	//std::cout<<"sizes: " << this->bondVals.size() << ", " << this->latsize*nnn << "\n";
	// a better way to do this is to explicitly use the getNeighbors() function to return a list
	// of neighbors -- in the event a site is given a short list of neighbors. maybe...
	for (int i=0; i<latsize;i++) {
		// for each site:
		for (int j=0; j<nnn; j++) {
			// for each neighbor:
			k=nnn*i + j;
			//std::cout<<"ks: " << this->bondVals.size() << ", " << this->latsize*nnn << ", " << k << "\n";
			// int to, int from, int * z, bool isocc) {
			//this->bonds[k]=bond(this->sites[i].getNeighbor(j), i, &this->bondVals[k], false);
			//
			//std::cout<< "edges (" << latsize << "): " << this->sites[i].getNeighbor(j) << ", " << i << ", " << k << ", " << (this->bondVals[k]) << "\n";
			//
			this->bonds.push_back(bond(this->sites[i].getNeighbor(j), i, &(this->bondVals[k]), false));
			//std::cout << "edgobj[" << k << "]: " << this->bonds[k].getto() << ", " << this->bonds[k].getfrom() << ", " << this->bonds[k].getval() << "\n";
			}
		}
	};

void lattice::setz0(int thisz) {
		// z0 is the threshold value for random numbers.
		this->z0=thisz;
	}

void lattice::setval(int index, int val) {
	//setz(index, val);
	this->siteVals[index] = val;
	}

void lattice::setbondval(int index, int val) {
	this->bondVals[index]=val;
	}

//void lattice::setz(int index, int val) {
//	siteVals[index] = val;
//	}

void lattice::setdlamode(int md) {
    this->dlamode=md;
	}

int lattice::getdlamode() {
    return this->dlamode;
    }

void lattice::writegrid(std::string fname) {
	//
	std::ofstream myfile;
	myfile.open(fname.c_str());
	myfile << "#lattype\t" << this->geometry << "\t" << this->geometrylist[this->geometry] << "\t" << this->width << "\n";
	myfile << "#i\tx\ty\tz\n";
	//
	for (int i=0;i<latsize;i++) {
                //myfile << i << "\t" << i%width << "\t" << i/width << "\t" << siteVals[i] << "\t"<< sites[i].getNeighbor(0) << "\t"<< sites[i].getNeighbor(1) << "\t"<< sites[i].getNeighbor(2) << "\t"<< sites[i].getNeighbor(3) <<"\n";
                myfile << i << "\t" << i%width << "\t" << i/width << "\t" << siteVals[i];
                for (int j=0; j<sites[i].getnNeighbors(); j++) {
                    myfile  << "\t" << sites[i].getNeighbor(j);
                }
                myfile << "\n";
            }
	myfile.close();
	
	}

void lattice::initializeCluster(int csite) {
	if (csite==NULL) {
		this->initializeCluster();
		}
	else {
		initializeCluster(std::vector<int>(1, csite));
		}
	};

void lattice::initializeCluster() {
	int L = this->getwidth();
	int i0 = (L*L + L)/2;

   if (this->geometry==3) {
		//i0 += ((this->latsize)*(this->latsize)*(this->latsize))/2;
		i0 = (L*L*L + L*L + L)/2;
		}
      //
   initializeCluster(i0);
	}

void lattice::initializeCluster(std::vector<int> csites) {
    //this->ipcluster = cluster(csites);	// initialize cluster member with provided list
    //this->ipcluster.elements = csites;
    //
    this->ipcluster = cluster();
    // update the cluster???
    // now, update the lattice.
    std::cout << "doing initialize cluster: " << csites.size() << ", " << csites[0] << "\n";
    for (unsigned int i=0; i<csites.size(); i++) {
        this->ipcluster.addDisplacedZ(this->siteVals[csites[i]]);
        this->ipcluster.addelement(csites[i]);
        siteVals[csites[i]] = this->zClust;	// 1 (maybe) being default "cluster" value. note the random values are all z>z0=10.
        // note: we've initiated with one or more sites, but we've not defined the bonds. let's just leave it at that.
        // if we're initializing with multiple sites, there is acctually no reason to assume they're connected.
        }
	//
	// update perimeter:
    //this->boundaries = getNeighbors(csites);
    this->boundaries = this->getNeighborBonds(csites);
    }

int lattice::getsite0() {
	return this->z0;
}

int lattice::getsiteval(int i) {
	return siteVals[i];
}

int * lattice::getsite(int i) {
	// equivalently (i think): return (siteVals + i);
	return &siteVals[i];
	}

int lattice::getbondval(int i) {
	return bondVals[i];
	}

int * lattice::getbond(int i) {
	return &bondVals[i];
	}

bond * lattice::getbondObj(int i) {
	return &bonds[i];
	}

int lattice::getlatsize() {
	//return latsize;
	return this->siteVals.size();
}
int lattice::getbondsize() {
	return this->bondVals.size();
	}

int lattice::getnNN() {
	return this->geomNN[this->geometry];
	}

int lattice::getgeometry() {
	return geometry;
}

int lattice::getwidth() {
	return width;
	}

int lattice::getrandSeed() {
	return this->randSeed;
	}

int lattice::randSeedPlus() {
	return this->randSeedPlus(1);
	}

int lattice::randSeedPlus(int n) {
	this->randSeed+=n;
	return this->randSeed;
	}

int lattice::getclustsize() {
	return this->ipcluster.getsize();
	}

std::map<int, int> lattice::getNeighbors(int k) {
	std::vector<int> tmpVector(1, k);
	return getNeighbors(tmpVector);
	}

std::map<int, int> lattice::getNeighbors() {
	return getNeighbors(this->ipcluster.elements);
}

std::map<int, int> lattice::getNeighbors(cluster X) {
	return getNeighbors(X.elements);
}

std::map<int, int> lattice::getNeighbors(std::vector<int>& positions) {
	std::map<int, int> myneighbors;
	std::vector<int> thisneighbor(nNeighbors);
	//
   //std::vector<int>::iterator it;
	//
	for (unsigned int i=0; i<positions.size(); i++) {
		thisneighbor = sites[positions[i]].getNeighbors();		// from sites in vector<int> form
		//
		for (unsigned int j=0; j<thisneighbor.size(); j++) {
			myneighbors.insert (std::pair<int, int>(siteVals[thisneighbor[j]], thisneighbor[j]) );		// aka, [val, index]
			//myneighbors.insert(thisneighbor.begin(), thisneighbor.end());
			//myneighbors.insert(thisneighbor.begin(), thisneighbor.end());
			}
		}
   // now, clean up neighbors... remove all cluster-value keys (or less than z0?)
   myneighbors.erase(this->zClust); 
	//
	/*
	for (unsigned int i=0; i<positions.size(); i++) {
		// remove cluster site indeces from neighbors
		// note (we can use positions[] because we just initialized this->ipcluster with it).
		myneighbors.erase(positions[i]);
		}
	*/

	return myneighbors;
}

std::map<int, int> lattice::getNeighborBonds() {
	return getNeighborBonds(this->ipcluster.elements);	// assume ipcluster
	}

std::map<int, int> lattice::getNeighborBonds(int k) {
	std::vector<int> tmp(1,k);
	return getNeighborBonds(tmp);	
	}

std::map<int, int> lattice::getNeighborBonds(cluster X) {
	return getNeighborBonds(X.elements);	
	}

std::map<int, int> lattice::getNeighborBonds(std::vector<int>& positions) {
	// assume all edge values are unique (see initialization process).
	std::map<int, int> myneighbors;
	int k;
	//std::vector<int> thisneighbor(this->getnNN());
	//
	for (int i=0; i<positions.size(); i++) {
		//thisneighbor = sites[positions[i]].getNeighbors();		// from sites in vector<int> form
		for (int j=0; j<this->getnNN(); j++) {
			// for each site, get neighbor bonds.
			k=this->getBondIndex(positions[i], j);
			myneighbors.insert(std::pair<int, int>(this->getbondval(k), k) );
			}
		//
		/*
		for (unsigned int j=0; j<thisneighbor.size(); j++) {
			myneighbors.insert (std::pair<int, int>(siteVals[thisneighbor[j]], thisneighbor[j]) );		// aka, [val, index]
			//myneighbors.insert(thisneighbor.begin(), thisneighbor.end());
			//myneighbors.insert(thisneighbor.begin(), thisneighbor.end());
			}
		*/
		}
   // now, clean up neighbors... remove all cluster-value keys (or less than z0?)
   // (for optimization; note this step returns only unoccupied bonds).
   myneighbors.erase(this->zClust); 
	//
	return myneighbors;
	}
std::vector< std::vector<int> > lattice::getTips() {
	// return a list of branch terminal edges
	// note that both sites and edges are set to zClust when they are occupied.
	//
	int tosite;
	std::vector< std::vector<int> > tips;
	std::vector<int> newtip(3);
	std::vector<int> neighborSites;
	
	int ntolinks=0;
	int nfromlinks=0;		// links/bonds to/from another site.
	//bool hastolink=false;
	//bool hasfromlink=false;
	int nNN = this->getnNN();	
	//
	for (int i=0; i< (this->ipcluster.edges.size()); i++) {
		tosite=this->ipcluster.edges[i][2];	// edges are [val, from, to]
		//neighborSites = this->getNeighbors(tosite);
		neighborSites = this->sites[tosite].getNeighbors();
		//neighborBonds = this->getNeighborBonds(tosite);
		// now, we see 1)if any of our neighborBonds are occupied, 2) if any neighbor sites have bonds to us...
		// note that getNeighborBonds() excludes occupied bonds, so we can determine whether or not
		// we're connected from the list length or by looking FOR our bond (!=zClust)
		//
		//if neighborSites.size()<(this->getnNN()) {
		// this implies that one or more links were excluded because they're occupied. however,
		// if we implement some sort of general connectivity, in which we establish boundaries, etc. based
		// on neighbor lists, this will create problems, so let's do it another way. let's just get the bond sites manually.
		int bondindex = tosite*nNN;
		//
		//for (it=neighborBonds.begin(); it!=neighborBonds.end(); ++it) {
		//for (it=(this->bonds.begin() + bondindex); it!=neighborBonds.end(); ++it) {
		ntolinks=0;
		for (int j=bondindex; j<(bondindex+nNN); j++) {
			// this site/bond cannot link to any other sites/bonds.
			if (this->bonds[j].getval()==this->zClust) {
				 // this site links TO another site (should be 0 for an end).
				 ntolinks++;
				 break;		// even one is too many...
				 //std::cout << "bnds1 (" << tosite << "):: " << this->bonds[j].getto() << ", " << this->bonds[j].getval()  << ", " << ntolinks << ")\n";
				 
				 }
			}
		
		//
		// now, check neighbors to see if they're bonding to this site:
		//for (it=neighborSites.begin(); it!=neighborSites.end(); ++it) {
		for (int k=0; k<neighborSites.size(); k++) {
			//int siteindex=it->second;	// as yet, site occupation is not considered.
			int siteindex=neighborSites[k];
			//neighborBonds=this->getNeighborBonds(siteindex);
			//for(it2=neighborBonds.begin(); it2!=neighborBonds.end(); ++it2) {
			int edgeindex=siteindex*nNN;
			nfromlinks=0;
			for (int j=edgeindex; j<(edgeindex+nNN); j++) {
				// is tosite in the neighboring site's neighborbonds?
				if (this->bonds[j].getto()==tosite and this->bonds[j].getval()==zClust) {
					// this site is linked FROM another site. for tips, there should be one "from" link (equal to thisedge.getfrom())
					nfromlinks++;
					// break;
					}
				}
			if (nfromlinks>1) break;
			}
		//
		//std::cout << "links(" << tosite << "): " << ntolinks << ", " << nfromlinks << "\n";
		if (ntolinks==0 and nfromlinks<=1) {
			// note: this allows for an isolated site.
			tips.push_back(this->ipcluster.edges[i]);
			}
		
		}
	return tips;
	};

std::map<int, int>* lattice::getboundaries() {
	return &boundaries;
	}

int lattice::ipstep() {
	//return lattice::ipstep(0);
	return lattice::ipstep(this->dlamode);
	}

//
int lattice::ipstep(int method) {
	// 1) find new "invaded" bont/site k from perimeters[] (either lowest or highest value)
	// 2) set siteVals[k] -> clustVal;
	// 3) set bondVals[k'] -> clustVal;
	// 4) add k to this->ipcluster[]  (sites)
	// 5) add k' to this->ipcluster[] (edges)
	// 6) add (new) neighbors of k to perimeter
	//
	// largest z-value will be at end of perimeters. note: .end() gives pointer to after last element.
	int newbondindex, newbondval, newsitetoIndex, newsitetoVal, newsitefromIndex, newsitefromval;
	
	if (this->ipcluster.getsize()<1) {
		this->initializeCluster();
		}
	std::map<int, int>::iterator it=boundaries.end();
	it--;
	//
	newbondindex = it->second;
	newbondval = it->first;		// boundary sites are sorted by this value.
	//std::cout<<"new bond: " << newbondindex << ", " << newbondval << "\n";
	//
	//int newsiteindex = it->second;
	//int newsiteval = it->first;		// boundary sites are sorted by this value.
	//
	newsitefromIndex = this->bonds[newbondindex].getfrom();
	newsitetoIndex = this->bonds[newbondindex].getto();
	//
	// is the "to" site occupied?
	while ( (this->getsiteval(newsitetoIndex)<z0 or (method==1 and nOcN(newsitetoIndex, 1)>1) ) and it!=boundaries.begin() ) {
		// if the "to" site for this bond is already occupied...
		it--;
		newbondindex = it->second;
		//newbondval = it->first;		// boundary sites are sorted by this value.
		//newsitefromIndex = this->bonds[newbondindex].getfrom();
		newsitetoIndex   = this->bonds[newbondindex].getto();
		}
	newbondval = it->first;
	newsitefromIndex = this->bonds[newbondindex].getfrom();
	newsitetoVal = this->getsiteval(newsitetoIndex);
	//
	// default is method=0;, but we don't really do anything with that.
	//
	/*
	if (method==1){
		// DLA: each new member has exactly 1 occupied neighbor.
		// nOcN(int k, kmax)
		while (nOcN(newsiteindex, 1)>1 and it!=boundaries.begin()) {
			it--;
			newsiteindex = it->second;
			newsiteval = it->first;
			}
		}
	*/
	//	
	// do we have a legal element? have we saturated the lattice; are there un-invaded sites available?
	if (newsitetoVal<this->z0) {
		//return NULL;
		std::cout<<"bonded an occupied site.\n";
		return -1;  
		}
	// we've selected a leagal element	
	//
	//boundaries.erase(newsiteval);	// note: this requires that our random numbers are unique, since boundaries is a map().
	boundaries.erase(newbondval);		// this is a map() operation.
	this->setval(newsitetoIndex, this->zClust);		// set the lattice element value to indicate "occupied"
	this->setbondval(newbondindex, this->zClust);
	//
	// add site-data to ipcluster()
	this->ipcluster.addelement(newsitetoIndex);
	this->ipcluster.addDisplacedZ(newsitetoVal);
	//
	// add edge (bond) data to ipcluster()
	this->ipcluster.addEdge(newbondval, newsitefromIndex, newsitetoIndex);		// this is equvalent to to addelement()+addDispz() for sites.
	//
	//std::vector<int> siteneighbors = sites[newsiteindex].getNeighbors();		// vector<int> of neighbor indeces.
	std::map<int, int> bondneighbors = this->getNeighborBonds(newsitetoIndex);
	//for (unsigned int i=0; i<bondneighbors.size(); i++) {
	for (std::map<int, int>::iterator it = bondneighbors.begin(); it != bondneighbors.end(); ++it) {
		// skip bonds with occupied sites...
		//if (this->getsiteval(bonds[(bondneighbors.begin()+i).second()].getto())==this->zClust) {
		if (this->getsiteval(bonds[(*it).second].getto())==this->zClust) {
			// the target site of this bond is occupied, so we can skip adding it.
			continue;
			}
		//boundaries.insert(std::pair<int, int>(*getsite(siteneighbors[i]), siteneighbors[i]) );
		//boundaries.insert(bondneighbors.begin()+i);
		boundaries.insert(*it);
		}
	
	// note: this will build a time-series of clusters and bonds added to the cluster. a bond-site
	// plot can be easily constructed from the bonds[] array.
	
	return newsitetoIndex;	// or we might return the bond index, or whatever. this will do probably.
	}
	

int lattice::ipstepSite() {
	//return lattice::ipstep(0);
	return lattice::ipstepSite(this->dlamode);
	}
//
int lattice::ipstepSite(int method) {
	// 1) find new "invaded" site k from perimeters[] (either lowest or highest value)
	// 2) set siteVals[k] -> clustVal;
	// 3) add k to this->ipcluster[]
	// 4) add (new) neighbors of k to perimeter
	//
	// 5) eventually loop over n steps or do n simultaneous steps?
	//
	// largest z-value will be at end of perimeters. note: .end() gives pointer to after last element.
	std::map<int, int>::iterator it=boundaries.end();
	it--;
	// we're really talking about bonds now...
	int newsiteindex = it->second;
	int newsiteval = it->first;		// boundary sites are sorted by this value.
	//
	//
	// default is method=0;, but we don't really do anything with that.
	//
	if (method==1){
		// DLA: each new member has exactly 1 occupied neighbor.
		// nOcN(int k, kmax)
		while (nOcN(newsiteindex, 1)>1 and it!=boundaries.begin()) {
		//while (nOcN(newsiteindex, 1)>1) {
			// note: == 1 is also legitimate, and nominally a nice bit of redundancy.
			// note also: this condition should never happen.
			//std::cout << "soft-skip... " << nOcN(newsiteindex, 1) << " \n";
			it--;
			newsiteindex = it->second;
			newsiteval = it->first;
			}
		/*	
		if (nOcN(newsiteindex, 1)>1 and it==boundaries.begin()) {
			// we failed to find a legal site.
			std::cout<<"hard-skip.\n";
			return NULL;
			}
		*/
		}
	//	
	// do we have a legal element? have we saturated the lattice; are there un-invaded sites available?
	if (newsiteval<this->z0) {
		//return NULL;
		return -1;  
		}
	//
	// we've selected a legale elment:	
	//
	boundaries.erase(newsiteval);	// note: this requires that our random numbers are unique, since boundaries is a map().
	this->setval(newsiteindex, this->zClust);		// set the lattice element value to indicate "occupied"
	this->ipcluster.addelement(newsiteindex);
	this->ipcluster.addDisplacedZ(newsiteval);
	//this->ipcluster.addedge() // but now that i think of it, edges don't really make sense in iperc, do they?
	//
	std::vector<int> siteneighbors = sites[newsiteindex].getNeighbors();		// vector<int> of neighbor indeces.
	for (unsigned int i=0; i<siteneighbors.size(); i++) {
		boundaries.insert(std::pair<int, int>(*getsite(siteneighbors[i]), siteneighbors[i]) );
		}
	
	
	return newsiteindex;
	}

int lattice::nOccNeighbors(int k) {
	return this->nOccNeighbors(k, this->nNeighbors);
	}
	
int lattice::nOccNeighbors(int k, int nmax) {
	// how many neighbors are occupied?
	// kmax: we may only need to know if n>= nmax; this will be good for speed.
	int nocc = 0;
	std::vector<int> Ns = this->sites[k].getNeighbors();
	//
	for (unsigned int i=0; i<Ns.size(); i++) {
		if (this->getsiteval(Ns[i])==this->zClust) {
			// z-value of neighbor-site indicates occupied.
			nocc++;
			if (nocc>=nmax) { break;}	// make it go faster
			}
		}
	return nocc;
	}

int lattice::nOcN(int k) {
	//
	return this->nOccNeighbors(k);
	}

int lattice::nOcN(int k, int nmax) {
	//
	return this->nOccNeighbors(k, nmax);
	}

int lattice::getBondIndex(int siteindex, int neighborindex) {
	int eindex=siteindex*this->getnNN() + neighborindex;
	return eindex;
	}

int lattice::makePercolatingCluster() {
   int thisdlamode = this->dlamode;
	if (thisdlamode==NULL) thisdlamode=0;
	return this->makePercolatingCluster(thisdlamode,0);
	}
int lattice::makePercolatingCluster(int mode) {
	return this->makePercolatingCluster(mode, 0);
	}

int lattice::makePercolatingCluster(int mode, int ptype) {
	// mode: 0: default, 1: DLA, etc.
	if (mode<0 or mode>1) {mode=1;}
	if (ptype==0) {ptype=1;}
	//if (refreshRand == 1) this->randSeedPlus();	// but at this time, this is out of place since the lattice is already initted.
	//
	// we might just skip this percolation type stuff and stick to one specific side...
	// for now, just ignore perc-type parameter.
	/*
	short ptop=0;
	short pbottom=0;
	short pleft=0;
	short pright=0;
	short pup = 0;	//3d
	short pdown = 0;
	*/
	//int newindex = this->ipcluster.getelement(0);
	int newindex;
	bool doit=true;
	//
	int w=this->width;
	//while doit==1 {
	//while ((ptop + pbottom + pleft + pright + pup + pdown)<ptype) {
	//while (newindex%w!=0 and newindex%w!=(w-1) and newindex/w!=0 and newindex/w!=(w-1)){
	while (doit==true) {
		newindex=this->ipstep(mode);
		//
		if ( this->geometry!=3 and (newindex%w==0 or newindex%w==(w-1) or newindex/w==0 or newindex/w==(w-1)) or (this->geometry==3 and ( newindex/(w*w) == 0 or newindex/(w*w)==(w-1) or newindex%w==0 or newindex%w==(w-1) or (newindex/w)%w==0 or (newindex/w)%w==(w-1) ) ) ) {
			//std::cout<< "it's percolated, " << newindex%w << ", " << newindex/w << ", " << (this->geometry==3 and newindex/(w*w)) << " \n";
			doit=false;
			}

		if (newindex>(this->latsize) or newindex<0) {
			break;
			}
		//std::cout
		//
		//std::cout<< "newIndex: " << newindex << ", %width: " << newindex%(this->width) << ", /width: " << newindex/(this->width) << "\n";
		//
		/*
		if (newindex%this->width == 0) {pleft=1;};
		if (newindex%this->width == (width-1)) {pright=1;};
		if (((newindex/width)%width) == 0) {ptop=1;};
		if (((this->width)*(this->width)) == (this->latsize/this->width-1)) {pbottom=1;};	// note: this is valid for 2 or 3 D.
		//
		if (this->geometry==3) {
			if (newindex/((this->width)*(this->width)) == 0) {pup=1;};
			if (newindex/((this->width)*(this->width)) == ((this->width)*(this->width))-1) {pdown=1;};
			}
		if (newindex>(this->latsize) or newindex<0) {
			break;	// somehow, we've gone off the grid...
			}
		//
		*/
		
		}
	//std::cout << "percolated at " << newindex << "\n";
	return this->ipcluster.elements.size();
	}

std::vector< std::vector<float> > lattice::getLatticeMap() {
	//
	std::vector< std::vector<float> > outvec(this->latsize);
	float x=NULL;
	float y=NULL;
	float val=NULL;
	//int* val;
	float z=0;
			
	for (int i=0; i<this->latsize; i++) {
		//switch (this->geometry) {
		if (this->geometry==0 or this->geometry==1) {
			// square or diamond
			//
			x = float(i%(this->width));
			y = float(i/(this->width));
			}
		if (this->geometry==2) {
			// triangular
			y = float(i/(this->width));
			int parity=int(y)%2;
			// dx = .5*((datalen/width)%2) // assume this or -(this) parity... (or opposite, aka odd/even advance by .5).
			x = float(i%width) + .5*float(parity);
			}
		if (this->geometry==0 or this->geometry==1) {
			// square or diamond
			//
			x = float(i%(this->width));
			y = float((i/(this->width))%(this->width));
			z = float(i/((this->width)*(this->width)));
			}
		//		
		val = float(*(this->getsite(i)));
		float tmp[5] = {float(i),x,y,z,val};
		//
		//std::cout<< "rwves set: " << x << ", " << y << ", " << z << ", " << val  << "\n";
		//outvec.push_back(std::vector<float>(tmp, (tmp+4)));
		outvec[i] = std::vector<float>(tmp, (tmp+5));
		//outvec.push_back(rwvec);
		//std::cout << "rowsize: " << outvec[i].size() << "\n";
		};
	return outvec;
	
	}
	
std::vector< std::vector<float> > lattice::getClusterMap() {
	//
	std::vector< std::vector<float> > outvec(this->ipcluster.getsize());
	float x=NULL;
	float y=NULL;
	float val=NULL;
	//int* val;
	float z=0;
	int k;		
	for (int i=0; i<this->ipcluster.getsize(); i++) {
		k=this->ipcluster.getelement(i);		// k is the index of a lattice.zValues site
		//
		if (this->geometry==0 or this->geometry==1) {
			// square or diamond
			//
			x = float(k%(this->width));
			y = float(k/(this->width));
			}
		if (this->geometry==2) {
			// triangular
			y = float(k/(this->width));
			int parity=int(y)%2;
			// dx = .5*((datalen/width)%2) // assume this or -(this) parity... (or opposite, aka odd/even advance by .5).
			x = float(k%width) + .5*float(parity);
			}
		if (this->geometry==3) {
			// square or diamond
			//
			int w=this->width;
			x = float(k%w);
			y = float((k/w)%w);
			z = float(k/(w*w));
			}
		//
		val = float(*(this->getsite(k)));
		float tmp[5] = {float(k),x,y,z,val};
		//outvec.push_back(std::vector<float>(tmp, (tmp+4)));
		outvec[i]=std::vector<float>(tmp, (tmp+5));
		};
	return outvec;
	}
	
	
	
	
	
	
	
	
	
	
