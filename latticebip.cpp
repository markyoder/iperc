#include "latticebip.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>


latticebip::lattice() {
	init(0, 0, 0, 1);
}

latticebip::lattice(int width0) {
	 init(width0, time(NULL), 0, 1);
}

latticebip::lattice(int width0, int randseed ) {
	 init(width0, randseed,  0, 1);
}

latticebip::lattice(int latlen, int randseed, int latgeom) {
	 init(width0, randseed, latgeom, 1);
}

latticebip::lattice(int width0, int randseed, int latgeom, int dlamode) {
	 init(width0, randseed, latgeom, dlamode);
}

void latticebip::init(int width0, int randseed, int latgeom){
	init(width0, randseed, latgeom, 1);
	}

//void latticebip::init(int latlen, int randseed, int* Zs, int latgeom){
void latticebip::init(int width0, int randseed, int latgeom, int dlamode){
	//
	if (randseed==0) {
		randseed = time(NULL);
		}
	//
	this->dlamode=dlamode;
	this->latsize=width0*width0;
	if (latgeom==3) this->latsize=width0*width0*width0; // cubic...
	#
	this->width=width0;
	this->randSeed=randseed;
	//srand(randseed);
	this->geometry = latgeom;	
	//
	// std::map<int, std::string> geometrylist;
	this->geometrylist[0]="square";
	this->geometrylist[1]="diamond";
	this->geometrylist[2]="triangular";
	this->geometrylist[3]="cubic";
	// 
	this->neighborcount[0]=4; //"square";
	this->neighborcount[1]=4; //"diamond";
	this->neighborcount[2]=6; //"triangular";
	this->neighborcount[3]=6; //"cubic";
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
	this->initializeLattice(randseed);
	//
	// define geometry:
	this->setgeometry(latgeom);
	this->ipcluster = cluster();	// a cluster object.
	//this->initializeCluster();
	//
	}

/*
latticebip::~lattice() {
	
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

void latticebip::initializeLattice() {
	//initializeLattice(time(NULL));
	int newrseed=this->randSeedPlus();
	initializeLattice(this->randSeed);
	//delete &newrseed;
	}
	
void latticebip::initializeLattice(int rseed) {
	int dummyvar=0;
	// it is quite important that both the indices and values for each lattice site are unique.
	// rather than be at the mercy of rand(), let's use a randomly ordered sequence of numbers.
	// use: z0 < z < latsize+z0 ordered randomly.
	srand(rseed);
	std::map<int, int> tmpZs;
	std::map<int, int>::iterator it;
	for (int i=0; i<latsize;i++) {
		// this map will be ordered by the first column (key), a random number.
		tmpZs.insert(std::pair<int, int>(rand(), this->z0+i));
		}
	//
	//for (int i=0; i<latsize;i++) {
	int i=0;
	for (it=tmpZs.begin(); it!=tmpZs.end(); ++it) {
		//siteVals[i]=z0+rand()%(RAND_MAX-z0);
		this->siteVals[i] = it->second;
		//std::cout << "this zval: " << siteVals[i] << "\n";
		i++;
		}
	//
	this->boundaries = std::map<int, int>();
	}
	
void latticebip::setrandSeed() {
	this->randSeed=time(NULL);
	}

void latticebip::setrandSeed(int x) {
	this->randSeed=x;
	}

void latticebip::setgeometry(int geom) {
	//
	this->geometry=geom;
	int w=width;
	int w2=width*width;
	//
	switch (geom) {
	case 0:{
      std:: cout << "setting square-geometry" << "\n";
      // square:
		nNeighbors = 4;
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
		nNeighbors = 4;
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
			nNeighbors = 6;
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
			this->nNeighbors = 6;
			std::cout<< "setting cubic geometry.\n";
			//
			int w3=width*width*width;
			int tmpNeighbors[6] = {-1, 1, -w, w, -w2, w2};
			//
			//std::cout << "begin setting sites and neighborses\n";
			for (unsigned int i=0; i<latsize;i++) {
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
	};

void latticebip::setz0(int thisz) {
		// z0 is the threshold value for random numbers.
		this->z0=thisz;
	}

void latticebip::setsite(int index, bool val) {
	//setz(index, val);
	siteVals[index] = val;
	}

//void latticebip::setz(int index, int val) {
//	siteVals[index] = val;
//	}

void latticebip::setdlamode(int md) {
    this->dlamode=md;
}

int latticebip::getdlamode() {
    return this->dlamode;
    }

void latticebip::writegrid(std::string fname) {
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

void latticebip::initializeCluster(std::vector<int> csites) {
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
        siteVals[csites[i]] = 1;	// 1 (maybe) being default "cluster" value. note the random values are all z>z0=10.
        
        }
	//
	// update perimeter:
    this->boundaries = getNeighbors(csites);
    }

void latticebip::initializeCluster(int csite) {
	if (csite==NULL) {
		this->initializeCluster();
		}
	else {
		initializeCluster(std::vector<int>(1, csite));
		}
	};

void latticebip::initializeCluster() {
	int L = this->getwidth();
	int i0 = (L*L + L)/2;

   if (this->geometry==3) {
		//i0 += ((this->latsize)*(this->latsize)*(this->latsize))/2;
		i0 = (L*L*L + L*L + L)/2;
		}
      //
   initializeCluster(i0);
}

int latticebip::getsite0() {
	return this->z0;
}

int latticebip::getsiteval(int i) {
	return siteVals[i];
}

int* latticebip::getsite(int i) {
	// equivalently (i think): return (siteVals + i);
	return &siteVals[i];
	}

int latticebip::getlatsize() {
	return latsize;
}

int latticebip::getgeometry() {
	return geometry;
}

int latticebip::getwidth() {
	return width;
	}

int latticebip::getrandSeed() {
	return this->randSeed;
	}

int latticebip::randSeedPlus() {
	return this->randSeedPlus(1);
	}

int latticebip::randSeedPlus(int n) {
	this->randSeed+=n;
	return this->randSeed;
	}

int latticebip::getclustsize() {
	return this->ipcluster.getsize();
	}

std::map<int, int> latticebip::getNeighbors(int k) {
	std::vector<int> tmpVector(1, k);
	return getNeighbors(tmpVector);
	}

std::map<int, int> latticebip::getNeighbors() {
	return getNeighbors(this->ipcluster.elements);
}

std::map<int, int> latticebip::getNeighbors(cluster X) {
	return getNeighbors(X.elements);
}

std::map<int, int> latticebip::getNeighbors(std::vector<int>& positions) {
	std::map<int, int> myneighbors;
	std::vector<int> thisneighbor(nNeighbors);
	//
   //std::vector<int>::iterator it;
	//
	for (unsigned int i=0; i<positions.size(); i++) {
		thisneighbor = sites[positions[i]].getNeighbors();		// from sites in vector<int> form
		//
		for (unsigned int j=0; j<thisneighbor.size(); j++) {
			myneighbors.insert (std::pair<int, int>(siteVals[thisneighbor[j]], thisneighbor[j]) );
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

std::map<int, int>* latticebip::getboundaries() {
	return &boundaries;
	}

int latticebip::ipstep() {
	//return latticebip::ipstep(0);
	return latticebip::ipstep(this->dlamode);
	}
//
int latticebip::ipstep(int method) {
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
	int newsiteindex = it->second;
	int newsiteval = it->first;
	//printf("ipstep: %d, %d\n", it->first, it->second);
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
	if (newsiteval<this->z0) { return NULL; }
	//
	// we've selected a legale elment:	
	//
	boundaries.erase(newsiteval);	// note: this requires that our random numbers are unique...
	this->setsite(newsiteindex, this->zClust);
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

int latticebip::nOccNeighbors(int k) {
	return this->nOccNeighbors(k, this->nNeighbors);
	}
	
int latticebip::nOccNeighbors(int k, int nmax) {
	// how many neighbors are occupied?
	// kmax: we may only need to know if n>= nmax; this will be good for speed.
	int nocc = 0;
	std::vector<int> Ns = this->sites[k].getNeighbors();
	//
	for (unsigned int i=0; i<Ns.size(); i++) {
		if (this->getsiteval(Ns[i])==this->zClust) {
			// z-value of neighbor-site indicates occupied.
			nocc++;
			if (nocc>nmax) { break;}	// make it go faster
			}
		}
	return nocc;
	}

int latticebip::nOcN(int k) {
	//
	return this->nOccNeighbors(k);
	}

int latticebip::nOcN(int k, int nmax) {
	//
	return this->nOccNeighbors(k, nmax);
	}

int latticebip::makePercolatingCluster() {
   int thisdlamode = this->dlamode;
	if (thisdlamode==NULL) thisdlamode=0;
	return this->makePercolatingCluster(thisdlamode,0);
	}
int latticebip::makePercolatingCluster(int mode) {
	return this->makePercolatingCluster(mode, 0);
	}

int latticebip::makePercolatingCluster(int mode, int ptype) {
	// mode: 0: default, 1: DLA, etc.
	if (mode<0 or mode>1) {mode=1;}
	if (ptype==0) {ptype=1;}
	//if (refreshRand == 1) this->randSeedPlus();	// but at this time, this is out of place since the lattice is already initted.
	//
	// we might just skip this percolation type stuff and stick to one specific side...
	short ptop=0;
	short pbottom=0;
	short pleft=0;
	short pright=0;
	short pup = 0;	//3d
	short pdown = 0;
	int newindex;
	//
	//while doit==1 {
	while ((ptop + pbottom + pleft + pright + pup + pdown)<ptype) {
		newindex=this->ipstep(mode);
		//std::cout
		//
		//std::cout<< newindex << ", %width: " << newindex%(this->width) << ", /width: " << newindex/(this->width) << "\n";
		//
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
		//if ((short(ptop) + short(pbottom) + short(pleft) + short(pright))>=ptype) {
		//	// we've percolated
		//	doit=0;
		//	}
		}
	//std::cout << "percolated at " << newindex << "\n";
	return this->ipcluster.elements.size();
	}

std::vector< std::vector<float> > latticebip::getLatticeMap() {
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
	
std::vector< std::vector<float> > latticebip::getClusterMap() {
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
	
	
	
	
	
	
	
	
	
	
