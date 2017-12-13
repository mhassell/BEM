// generate a mesh for observing potentials

#include "meshPolygon.hpp"
#include "inPolygon.hpp"
#include "matrixRoutines.hpp"

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

	double xstep = (xmax - xmin)/(double)nx;
	double ystep = (ymax - ymin)/(double)ny;

	std::cout << xstep << std::endl;

	// linspace the box
	for(int i = 0; i < nx + 1; i++){
		xpts[i] = xmin + i*xstep;
	}

	for(int i = 0; i < ny + 1; i++){
		ypts[i] = ymin + i*ystep;
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
	int nElts = elements.rows();	// nElts = nNode, but be careful!
	int count = 0;
	int next = 0;

	while(count != nElts){
		polygon[next].x = coordinates(elements(next,0),0);
		polygon[next].y = coordinates(elements(next,0),1);
		
		// find where my next node is
		for(int i = 0; i < nElts; i++){
			if(elements(next,1)==elements(i,0)){
				next = i;
				break;
			}
		}
		count++;
	}

	// now check for points in the polygon
	Point p;
	 
	for(int i = 0; i < nx+1; i++){
		p.x = xpts[i];
		for(int j = 0; j < ny+1; j++){
			p.y = ypts[j];
			// std::cout << p.x << " " << p.y << std::endl;
			if(isInside(polygon, nElts, p)){
				// std::cout << "inside" << std::endl;	
			}
		}
	}

	delete[] xpts;
	delete[] ypts;
	delete[] polygon;

}




