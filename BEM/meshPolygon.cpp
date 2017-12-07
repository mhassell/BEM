// generate a mesh for observing potentials

#include "meshPolygon.hpp"
#include <Eigen/Dense>
#include <cassert>

/*

input: 	geometry g
		double* box: [xmin xmax ymin ymax]
		double h: spacing from the boundary
		int nx, ny : number of points in each direction

output: double** mesh of points (x,y)

*/

// constructor: just copy the essential fields from geometry & expand
mesh::mesh(const geometry& g){

	

}

mesh::meshPolygon(const geometry& g, double* box, double h, int nx, int ny){

	// check that mesh params and the box are reasonable
	assert(h!=0);
	assert(xmax>xmin && ymax>ymin);

	double xmin = box[0];
	double xmax = box[1];
	double ymin = box[2];
	double ymax = box[3];

	double* xpts = new double[nx+1];
	double* ypts = new double[ny+1];

	double xstep = xmax - xmin;
	double ystep = ymax - ymin;

	for(int i = 0; i < nx + 1; i++){
		xpts[i] = xmin + i*xstep;
	}

	for(int i = 0; i < ny + 1; i++){
		ypts[i] = ymin + i*ystep;
	}

	for(int i = 0; i < nElts; i++){
		
	}

}

void mesh::rayCasting(){

}
