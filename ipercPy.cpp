#include <Python.h>
#include <iostream>

#include "lattice.h"
#include "bond.h"


// declarations:
static PyObject * ipercPy_echo(PyObject *self, PyObject *args);
static PyObject * ipercPy_ipercObj(PyObject *self, PyObject *args);
static PyObject * ipercPy_makePercolatingCluster(PyObject *self, PyObject *args);
static PyObject * ipercPy_ipstep(PyObject *self, PyObject *args);
static PyObject * ipercPy_initializeLattice(PyObject *self, PyObject *args);
static PyObject * ipercPy_initializeCluster(PyObject *self, PyObject *args);
static PyObject * ipercPy_randSeedPlus(PyObject *self, PyObject *args);
static PyObject * ipercPy_init(PyObject *self, PyObject *args);

static PyObject * ipercPy_setrandSeed(PyObject *self, PyObject *args);
static PyObject * ipercPy_setdlamode(PyObject *self, PyObject *args);

static PyObject * ipercPy_getrandSeed(PyObject *self, PyObject *args);
static PyObject * ipercPy_getwidth(PyObject *self, PyObject *args);
static PyObject * ipercPy_getgeom(PyObject *self, PyObject *args);
static PyObject * ipercPy_getLatSize(PyObject *self, PyObject *args);
static PyObject * ipercPy_getClusterSites(PyObject *self, PyObject *args);
static PyObject * ipercPy_getDisplacedZs(PyObject *self, PyObject *args);
static PyObject * ipercPy_getClusterBonds(PyObject *self, PyObject *args);
static PyObject * ipercPy_getClustSize(PyObject *self, PyObject *args);
static PyObject * ipercPy_getLatticeSites(PyObject *self, PyObject *args);
static PyObject * ipercPy_getLatticeMap(PyObject *self, PyObject *args);
static PyObject * ipercPy_getClusterMap(PyObject *self, PyObject *args);
static PyObject * ipercPy_getdlamode(PyObject *self, PyObject *args);
static PyObject * ipercPy_getLatticeBonds(PyObject *self, PyObject *args);
static PyObject * ipercPy_getboundaries(PyObject *self, PyObject *args);
static PyObject * ipercPy_getnNN(PyObject *self, PyObject *args);
static PyObject * ipercPy_getTips(PyObject *self, PyObject *args);

//
lattice globLat;

// ipercPy_initializeLattice   setrandSeed

// code:
// ipercPy_getClusterBonds

static PyMethodDef ipercMethods[] = {
    {"ipecho",  ipercPy_echo, METH_VARARGS, "Echo the hard way."},
    {"ipercObj", ipercPy_ipercObj, METH_VARARGS, "Initialize ipercolation object"},
    {"makePercolatingCluster", ipercPy_makePercolatingCluster, METH_VARARGS, "Make a percolating cluster"},    
    {"ipstep", ipercPy_ipstep, METH_VARARGS, "Do an invasion percolation step; return site or edge or something."},
    {"initializeLattice", ipercPy_initializeLattice, METH_VARARGS, "(Re)initialize lattice (set new random numbers)."},
    {"initializeCluster", ipercPy_initializeCluster, METH_VARARGS, "(Re)initialize cluster object."},
    {"setrandSeed", ipercPy_setrandSeed, METH_VARARGS, "Set new random number seed."},
    {"randSeedPlus", ipercPy_randSeedPlus, METH_VARARGS, "Increment random number seed by +1"},
    {"getrandSeed", ipercPy_getrandSeed, METH_VARARGS, "Return random number seed"},
    {"getClusterSites", ipercPy_getClusterSites, METH_VARARGS, "Return indeces of cluster sites."},
    {"getLatticeSites", ipercPy_getLatticeSites, METH_VARARGS, "Return indeces of cluster sites."},
	 {"getLatticeBonds", ipercPy_getLatticeBonds, METH_VARARGS, "Return lattice bonds/edges."},
    {"getDisplacedZs", ipercPy_getDisplacedZs, METH_VARARGS, "Return displaced Zs."},
    {"getClusterBonds", ipercPy_getClusterBonds, METH_VARARGS, "Return displaced Zs."},
    {"getLatSize", ipercPy_getLatSize, METH_VARARGS, "Return lattice size (int)."},
    {"getwidth", ipercPy_getwidth, METH_VARARGS, "Return lattice width (int)."},
    {"getgeom", ipercPy_getgeom, METH_VARARGS, "Return lattice geometry (int)."},
    {"getClustsize", ipercPy_getClustSize, METH_VARARGS, "Return lattice width (int)."},
    {"getLatticemap", ipercPy_getLatticeMap, METH_VARARGS, "Return [x,y,z,val] from lattice."},
    {"getClustermap", ipercPy_getClusterMap, METH_VARARGS, "Return [x,y,z,val] from lattice."},
    {"init", ipercPy_init, METH_VARARGS, "Initialize lattice (execute init())."},
    {"getdlamode", ipercPy_getdlamode, METH_VARARGS, "Initialize lattice (execute init())."},
    {"setdlamode", ipercPy_setdlamode, METH_VARARGS, "Initialize lattice (execute init())."},
    {"getboundaries", ipercPy_getboundaries, METH_VARARGS, "return boundary values."},
    {"getnNN", ipercPy_getnNN, METH_VARARGS, "return number of neighbors (for bond-site conversions)."},
    {"getTips", ipercPy_getTips, METH_VARARGS, "return a list of terminal edges (branch tips)."},
    {NULL, NULL, 0, NULL} 
    };

//     
//    

static PyObject * ipercPy_initializeLattice(PyObject *self, PyObject *args) {
	globLat.initializeLattice();
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject * ipercPy_initializeCluster(PyObject *self, PyObject *args) {
	//const int *x;
	/*
	if (!PyArg_ParseTuple(args, "s", &x)) {
        return NULL;
        }
   */
   int i0 = NULL;
   int ok;
   ok = PyArg_ParseTuple(args, "|l", &i0);
   if (i0==NULL) {
		globLat.initializeCluster();
		}
	else {
		globLat.initializeCluster(i0);
		};
	//
	return Py_None;
	}
	
static PyObject * ipercPy_setrandSeed(PyObject *self, PyObject *args) {
	const int *x;
	if (!PyArg_ParseTuple(args, "s", &x)) {
        return NULL;
        }
	globLat.setrandSeed(*x);
	//
	return Py_None;
	}

static PyObject * ipercPy_setdlamode(PyObject *self, PyObject *args) {
        //const int *x;
        int ok;
        int dlamode;
        ok = PyArg_ParseTuple(args, "|l", &dlamode);
        globLat.setdlamode(dlamode);
        //
        return Py_None;
        }
static PyObject * ipercPy_getdlamode(PyObject *self, PyObject *args) {
    int x = globLat.getdlamode();
    return Py_BuildValue("i", x);
}
// ipercPy_getnNN
static PyObject * ipercPy_getnNN(PyObject *self, PyObject *args) {
    int x = globLat.getnNN();
    return Py_BuildValue("i", x);
}

static PyObject * ipercPy_randSeedPlus(PyObject *self, PyObject *args) {
	int x = globLat.randSeedPlus();
	return Py_BuildValue("i", x);
	}

static PyObject * ipercPy_getrandSeed(PyObject *self, PyObject *args) {
	int x = globLat.getrandSeed();
	return Py_BuildValue("i", x);
	}

static PyObject * ipercPy_echo(PyObject *self, PyObject *args)
{
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        // parse args, expect "s"tring format, to address of command ???
        return NULL;
    //sts = system(command);
    sts=42;
    std::cout<<command << ", " << sts << "\n";
    //
    /*
    // more sophisticated error handling:
    if (sts < 0) {
        PyErr_SetString(SpamError, "System command failed");
        return NULL;
    }
    return PyLong_FromLong(sts);
    */
    
    return Py_BuildValue("i", sts);
}

static PyObject * ipercPy_init(PyObject *self, PyObject *args) {
	// for now, just return a fixed prammed class object?
	//
	int ok;
	int z=0;
	int thiswidth=256;
	int randseed=time(NULL);
	int latgeom=0;
	int dlamode=1;
   int percmode=0;
   int drs=0;
   //
   // width0, int randseed, int latgeom
   //
   ok = PyArg_ParseTuple(args, "ll|l", &thiswidth, &latgeom, &randseed);
   if (randseed==0) randseed=time(NULL);
   //
   std::cout << "width: " << thiswidth << ", geom: " << latgeom << ", rseed: " << randseed << "\n";
    globLat.init(thiswidth, randseed, latgeom);
    globLat.setdlamode(dlamode);
	
    return Py_BuildValue("i", z);

	}

static PyObject * ipercPy_ipercObj(PyObject *self, PyObject *args) {
	// for now, just return a fixed prammed class object?
	//
	int ok;
	int z=0;
	int thislatlen=256;
	int randseed=time(NULL);
	int latgeom=0;
	int dlamode=1;
   int drs=0;
   //
   ok = PyArg_ParseTuple(args, "lll|l", &thislatlen, &latgeom, &dlamode, &randseed);
   if (randseed==0) randseed=time(NULL);
   //
   std::cout << "width: " << thislatlen << ", geom: " << latgeom << ", dla: " << dlamode << ", rseed: " << randseed << "\n";
   //
   //lattice* thislat;
   //thislat = &globLat;
	//lattice::lattice lat = lattice::lattice(thislatlen, randseed, latgeom);
	//
   globLat = lattice(thislatlen, randseed+drs, latgeom);
  // thislat->init(thislatlen, randseed+drs, latgeom);
  int i0 = (thislatlen*thislatlen + thislatlen)/2;
  if (latgeom==3) {i0+=thislatlen*thislatlen*thislatlen/2;}
   
   globLat.initializeCluster(i0);
   globLat.setdlamode(dlamode);
	//globLat = thislat;
	//z = globLat.makePercolatingCluster(dlamode, percmode);
	//globLat.writegrid(std::string("gridout-square-dla2.dat"));
	
	//z = globLat.makePercolatingCluster(dlamode, percmode);
	//globLat.writegrid(std::string("gridout-square-dla2.dat"));
	
   return Py_BuildValue("i", z);
}

static PyObject * ipercPy_makePercolatingCluster(PyObject *self, PyObject *args) {
    // presently, this function (correctly) always returns the same cluster. in order to fix this, we need
    // to reinitialize (or at least reset the values of the lattice).
    //
    int ok;
    int dlamode=1;
    //
    ok = PyArg_ParseTuple(args, "|l", &dlamode);
    //if (randseed==0) randseed=time(NULL);
    std::cout<<"making percolationg cluster, dlamode= " << dlamode << "\n";
    //
    int z=globLat.makePercolatingCluster(dlamode);
    //globLat.writegrid(std::string("gridout-square-dla2.dat"));
    return Py_BuildValue("i", z);
    }
static PyObject * ipercPy_ipstep(PyObject *self, PyObject *args) {
    // presently, this function (correctly) always returns the same cluster. in order to fix this, we need
    // to reinitialize (or at least reset the values of the lattice).
    //
    int ok;
    int dlamode=1;
    //
    ok = PyArg_ParseTuple(args, "|l", &dlamode);
    //if (randseed==0) randseed=time(NULL);
    std::cout<<"doing ipstep, dlamode= " << dlamode << "\n";
    //
    int z=globLat.ipstep(dlamode);
    //globLat.writegrid(std::string("gridout-square-dla2.dat"));
    return Py_BuildValue("i", z);
    }
    
static PyObject * ipercPy_getClusterSites(PyObject *self, PyObject *args) {
	
	PyObject *pylist;
	PyObject *item;
	std::vector<int> C = globLat.ipcluster.getelements();
	pylist = PyTuple_New(C.size());
	for (unsigned int i=0; i < C.size(); i++) {
		item=PyInt_FromLong(C.at(i));
		PyTuple_SetItem(pylist, i, item);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;
	
	//return Py_BuildValue("i", 0);
	};

// ipercPy_getboundaries
static PyObject * ipercPy_getboundaries(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *pyrow;
	PyObject *item;
	std::map<int, int> * M = globLat.getboundaries();
	
	pylist = PyTuple_New(M->size());
	int i=0;
	for (std::map<int, int>::iterator it=M->begin(); it!= M->end(); ++it) {
		pyrow = PyTuple_New(2);
		PyTuple_SetItem(pyrow, 0, PyInt_FromLong((*it).second));	// index
		PyTuple_SetItem(pyrow, 1, PyInt_FromLong((*it).first));	// val
		PyTuple_SetItem(pylist, i, pyrow);
		i++;
		//
		}
	return pylist;
	
	//return Py_BuildValue("i", 0);
	};
	
static PyObject * ipercPy_getClusterBonds(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *pyrow;
	PyObject *item;
	std::vector<std::vector<int> > C = globLat.ipcluster.getedges();
	std::vector<int> row();
	int rowsize=3;	// assume this is right...
	//
	pylist = PyTuple_New(C.size());
	for (unsigned int i=0; i < C.size(); i++) {
		//int rowsize=C[i].size();
		//if (rowsize==0) continue;
		//std::cout << "rowsize: " << rowsize << "\n";
		pyrow  = PyTuple_New(rowsize);
		//
		// from, to, val.
		PyTuple_SetItem(pyrow, 0, PyInt_FromLong(C[i][1]));
		PyTuple_SetItem(pyrow, 1, PyInt_FromLong(C[i][2]));
		PyTuple_SetItem(pyrow, 2, PyInt_FromLong(C[i][0]));
		//
		PyTuple_SetItem(pylist, i, pyrow);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;
	
	//return Py_BuildValue("i", 0);
	};

static PyObject * ipercPy_getClusterBond(PyObject *self, PyObject *args) {
	PyObject *pyrow;
	PyObject *item;
   int ok;
   int i=0;
   //
   ok = PyArg_ParseTuple(args, "|l", &i);
   std::vector<int> R = globLat.ipcluster.getedge(i);
   pyrow = PyTuple_New(R.size());
   for (unsigned int j=0; j<R.size(); j++) {
   	item = PyInt_FromLong(R[j]);
   	PyTuple_SetItem(pyrow, j, item);
   	}
   return pyrow;
   }


static PyObject * ipercPy_getLatticeBonds(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject * pyrow;
	PyObject *item;
	bond * B;
	int bsize = globLat.getbondsize();
	//
	pylist = PyTuple_New(bsize);
	for (int i=0; i < bsize; i++) {
		pyrow = PyTuple_New(3);
		B=globLat.getbondObj(i);
		//
		PyTuple_SetItem(pyrow, 0, PyInt_FromLong(B->getfrom()));
		PyTuple_SetItem(pyrow, 1, PyInt_FromLong(B->getto()));
		PyTuple_SetItem(pyrow, 2, PyInt_FromLong(B->getval()));
		//
		
		PyTuple_SetItem(pylist, i, pyrow);
		}
	return pylist;
	
	//return Py_BuildValue("i", 0);
	};

static PyObject * ipercPy_getLatticeSites(PyObject *self, PyObject *args) {
	
	PyObject *pylist;
	PyObject *item;
	int * start = globLat.getsite(0);
	pylist = PyTuple_New(globLat.getlatsize());
	for (int i=0; i < globLat.getlatsize(); i++) {
		item=PyInt_FromLong(*(start + i));
		PyTuple_SetItem(pylist, i, item);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;
	
	//return Py_BuildValue("i", 0);
	};
	
static PyObject * ipercPy_getDisplacedZs(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *item;
	std::vector<int> C = globLat.ipcluster.getDisplacedZs();
	pylist = PyTuple_New(C.size());
	int csize=C.size();
	for (unsigned int i=0; i < csize; i++) {
		item=PyInt_FromLong(C.at(i));
		PyTuple_SetItem(pylist, i, item);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;
	};
	
static PyObject * ipercPy_getwidth(PyObject *self, PyObject *args) {
	int w = globLat.getwidth();
	return Py_BuildValue("i", w);
	};
	
static PyObject * ipercPy_getgeom(PyObject *self, PyObject *args) {
	int g = globLat.getgeometry();
	return Py_BuildValue("i", g);
	};
	
static PyObject * ipercPy_getLatSize(PyObject *self, PyObject *args) {
	int w = globLat.getlatsize();
	return Py_BuildValue("i", w);
	};


static PyObject * ipercPy_getClustSize(PyObject *self, PyObject *args) {
	int w = globLat.ipcluster.getsize();
	return Py_BuildValue("i", w);
	};

// percPy_getLatticeMap
static PyObject * ipercPy_getLatticeMap(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *pyrow;
	PyObject *item;
	std::vector< std::vector<float> > C = globLat.getLatticeMap();
	int csize=C.size();
	pylist = PyTuple_New(csize);
	for (unsigned int i=0; i < csize; i++) {
		int rowsize=C[i].size();
		if (rowsize==0) continue;
		//std::cout << "rowsize: " << rowsize << "\n";
		pyrow  = PyTuple_New(rowsize);	// should be 4...
		for (unsigned int j=0; j < rowsize; j++) {
			item=PyFloat_FromDouble(double(C[i][j]));
			PyTuple_SetItem(pyrow, j, item);
			}
		//
		//PyTuple_SetItem(pylist, i, item);
		PyTuple_SetItem(pylist, i, pyrow);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;
	//return pyrow;	
	}

static PyObject * ipercPy_getClusterMap(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *pyrow;
	PyObject *item;
	std::vector< std::vector<float> > C = globLat.getClusterMap();
	int csize=C.size();
	pylist = PyTuple_New(csize);
	for (unsigned int i=0; i < csize; i++) {
		int rowsize=C[i].size();
		if (rowsize==0) continue;
		pyrow  = PyTuple_New(rowsize);	// should be 4...
		for (unsigned int j=0; j < rowsize; j++) {
			item=PyFloat_FromDouble(double(C[i][j]));
			PyTuple_SetItem(pyrow, j, item);
			}
		//
		//PyTuple_SetItem(pylist, i, item);
		PyTuple_SetItem(pylist, i, pyrow);
		//Py_DECREF(C.at(i));	// i have no idea what this does.
		}
	return pylist;	
	//return pyrow;
	}

static PyObject * ipercPy_getTips(PyObject *self, PyObject *args) {
	PyObject *pylist;
	PyObject *pyrow;
	PyObject *item;
	std::vector<std::vector<int> > M = globLat.getTips();
	//
	pylist = PyTuple_New(M.size());
	//int i=0;
	//for (std::map<int, int>::iterator it=M->begin(); it!= M->end(); ++it) {
	for (int i=0; i<M.size(); i++) {
		pyrow = PyTuple_New(3);
		PyTuple_SetItem(pyrow, 0, PyInt_FromLong(M[i][0]));	// from
		PyTuple_SetItem(pyrow, 1, PyInt_FromLong(M[i][1]));	// to
		PyTuple_SetItem(pyrow, 2, PyInt_FromLong(M[i][2]));	// val
		PyTuple_SetItem(pylist, i, pyrow);
		//
		}
	return pylist;
	}

/*
PyMODINIT_FUNC
initspam(void)
{
    PyObject *m;

    m = Py_InitModule("spam", SpamMethods);
    if (m == NULL)
        return;

    SpamError = PyErr_NewException("spam.error", NULL, NULL);
    Py_INCREF(SpamError);
    PyModule_AddObject(m, "error", SpamError);
}
*/


PyMODINIT_FUNC
initipercPy(void)
{
    (void) Py_InitModule("ipercPy", ipercMethods);
}

int main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initipercPy();
    }

