// generate a mesh for observing potentials

#include "meshPolygon.hpp"
#include "inPolygon.hpp"
#include "matrixRoutines.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include <fstream>

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

	double* xpts = new double[nx];
	double* ypts = new double[ny];

	double xstep = (xmax - xmin)/(double)nx;
	double ystep = (ymax - ymin)/(double)ny;

	// linspace the box
	for(int i = 0; i < nx; i++){
		xpts[i] = xmin + i*xstep;
	}

	for(int i = 0; i < ny; i++){
		ypts[i] = ymin + i*ystep;
	}

	// move the coordinates in the normal direction by h
	/*
	for(int i = 0; i < elements.rows(); i++){
		coordinates(elements(i,0),0) += h*normals(elements(i,0),0);
		coordinates(elements(i,1),0) += h*normals(elements(i,1),0);
		coordinates(elements(i,0),1) += h*normals(elements(i,0),1);
		coordinates(elements(i,1),1) += h*normals(elements(i,1),1);	
	}
	*/

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
	for(int i = 0; i < nx; i++){
		p.x = xpts[i];
		for(int j = 0; j < ny; j++){
			p.y = ypts[j];
			Xpts.push_back(p.x);
			Ypts.push_back(p.y);
			if(isInside(polygon, nElts, p)){
				Inside.push_back(true);
			}
			else{
				Inside.push_back(false);
			}
		}
	}

	delete[] xpts;
	delete[] ypts;
	delete[] polygon;

}

// eventually put a string argument for the filename
void mesh::saveMesh(){

	// make sure we don't error and fail to close the file
	assert(Xpts.size()==Ypts.size());

	std::ofstream file("mesh.csv");
	for(int i = 0; i < Xpts.size(); i++){
		file << Xpts[i] << "," << Ypts[i] << "," << Inside[i] << std::endl;
	}

	file.close();

}

