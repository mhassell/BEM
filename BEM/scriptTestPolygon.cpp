// script to test the mesh class and polygon functionality

#include <iostream>
#include <Eigen/Dense>

#include "meshPolygon.hpp"
#include "geometry.hpp"

int main(){

	// make the geometry
	Eigen::MatrixXd coords(4,2);
    Eigen::MatrixXi elts(4,2);

	coords << 	0, 	0,
				1, 	0,
			  	0.8, 0.8, 
			  	0.2, 0.6;

	elts << 	0, 	1,
				2,	3,
				1, 	2,
				3, 	0;
	
	geometry g(coords, elts);

	mesh myMesh(g);

	double box[4] = {-2, 2, -2, 2};
	double h = 0.01;
	int nx = 100;
	int ny = 100;

	myMesh.meshPolygon(box,h,nx,ny);

}
