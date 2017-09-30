#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#include "geometry.hpp"
#include "quadTables.hpp"
#include "testPotentials.hpp"

double kernel(double);

int main(){

	// polynomial degree
	int k = 0;

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
	//g.refine();

	// quadrature
	Eigen::MatrixXd q1d = tableGauss(5);	  

	// kernel
	double (*ker)(double) = &kernel;

	// obs pts
	Eigen::MatrixXd obs(5,2);

	obs << 1, 2,
		   2, 3,
		  -1.5, 0.8,
		  -1.1, -0.5,
		   1.5, -0.3;	

	// Eigen::MatrixXd SL = testPotentialXh(g, ker, obs, k, q1d);
	Eigen::MatrixXd DL = testPotentialYh(g, ker, obs, k, q1d);

}

double kernel(double x){

	return log(x);


}
