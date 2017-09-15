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
	// g.refine();

	// quadrature
	Eigen::MatrixXd q1d = tableGauss(63);	  

	double (*ker)(double) = &kernel;

}

double kernel(double x){

	return log(x);

}
