# since don't use a default MAKE filename, execute via:
# make -f makeperc
###
# set variables:
#
CXX=g++ -Wall
#INCLUDES=-I/home/myoder/Documents/Source -I/home/myoder/Documents/gnuplot_i++
INCLUDES=-I/home/myoder/Documents/Source
#LDFLAGS=-L/usr/libs
LIBS=-L/usr/libs 
CXXFLAGS=$(INCLUDES) 

all: iperc1 iperc2

iperc1.o: iperc1.cpp site.o lattice.o
	$(CXX) iperc1.cpp site.o lattice.o -c

iperc2.o: iperc2.cpp site.o lattice.o
	$(CXX) iperc2.cpp site.o lattice.o -c


site.o: site.cpp site.h cluster.o
	$(CXX) $(FLAGS) site.cpp -c

lattice.o: lattice.cpp lattice.h site.o cluster.o
	$(CXX) $(FLAGS) lattice.cpp site.o cluster.o -c
	
cluster.o: cluster.cpp cluster.h
	$(CXX) $(FLAGS) cluster.cpp -c
	
iperc1: iperc1.o site.o lattice.o cluster.o
	$(CXX) iperc1.o site.o lattice.o cluster.o -o iperc1 $(LIBS) $(INCLUDES)
	
iperc2: iperc2.o site.o lattice.o cluster.o
	$(CXX) iperc2.o site.o lattice.o cluster.o -o iperc2 $(LIBS) $(INCLUDES)	

clean:
	rm *.o
	rm iperc1
	rm iperc2
