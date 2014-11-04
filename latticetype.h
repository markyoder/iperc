#ifndef LATTICETYPE_H
#define LATTICETYPE_H

class latticetype
{
protected:
	int nNeighbors;
	int *nDeltas;		// array of delta-Ns to define neighboring sites.
	int geometry;		// geometry of lattice (triangular, square, etc.)

public:
	latticetype();
	latticetype(int);
	//
	void latticetype::init(int ltype)
	//
	int getnNeighbors();
	int getgeometry();
	//
	void setnNeighbors(int);		// but you'd want to use this carefully; nominally, this should be set with the geometry package.
	void setndeltas(int*);
	void setgeometry(int);
 }
#endif // LATTICETYPE_H
