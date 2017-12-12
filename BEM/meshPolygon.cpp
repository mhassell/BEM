// generate a mesh for observing potentials

#include "meshPolygon.hpp"
#include "inPolygon.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include "omp.h"

// constructor: just copy the essential fields from geometry & expand
mesh::mesh(const geometry& g){

	elements = g.elements;
	coordinates = g.coordinates;
	normals = g.normals;

}

/*
input: 	geometry g
		double* box: [xmin xmax ymin ymax]
		double h: spacing from the boundary
		int nx, ny : number of points in each direction

output: double** mesh of points (x,y) in the mesh class
*/

void mesh::meshPolygon(double* box, double h, int nx, int ny){

	double xmin = box[0];
	double xmax = box[1];
	double ymin = box[2];
	double ymax = box[3];

	// check that mesh params and the box are reasonable
	assert(h!=0);
	assert(xmax>xmin && ymax>ymin);

	double* xpts = new double[nx+1];
	double* ypts = new double[ny+1];

	double xstep = xmax - xmin;
	double ystep = ymax - ymin;

	// linspace the box
	for(int i = 0; i < nx + 1; i++){
		xpts[i] = xmin + i*xstep;
	}

	for(int i = 0; i < ny + 1; i++){
		ypts[i] = ymin + i*ystep;
	}

	// now tensorize the box
	double **Xgrid = new double*[nx+1];
	for(int i = 0; i < nx + 1; i++){
		Xgrid[i] = new double[ny+1];
	}

	double **Ygrid = new double*[ny+1];
	for(int i = 0; i < ny + 1; i++){
		Ygrid[i] = new double[nx+1];
	}

	for(int i = 0; i < nx+1; i++){
		for(int j = 0; j < ny+1; j++){
			Xgrid[i][j] = xpts[i];
		}
	}

	for(int i = 0; i < ny+1; i++){
		for(int j = 0; j < nx+1; j++){
			Ygrid[i][j] = ypts[i];
		}
	}

	// move the coordinates in the normal direction by h
	for(int i = 0; i < elements.rows(); i++){
		coordinates(elements(i,0),0) *= (1+h)*normals(elements(i,0),0);
		coordinates(elements(i,1),0) *= (1+h)*normals(elements(i,1),0);
		coordinates(elements(i,0),1) *= (1+h)*normals(elements(i,0),1);
		coordinates(elements(i,1),1) *= (1+h)*normals(elements(i,1),1);
	}

	// make a polygon as an array of points
	Point* polygon = new Point[elements.rows()];
	int nElts = elements.rows();	
	int count = 0;
	int next = 0;

	while(count != nElts){
		polygon[next].x = coordinates(elements(next,0),0);
		polygon[next].y = coordinates(elements(next,0),1);
		next = elements(next,1);
		std::cout << next << std::endl;		
		count++;
	}

}




