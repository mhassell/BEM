#include <iostream>
#include <Eigen/Dense>
#include <math.h>


#include "geometry.hpp"
#include "testsAndProjections.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"
#include "quadTables.hpp"


double f(double, double);

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

	// quadrature
	Eigen::MatrixXd q1d = tableGauss(63);	  

	// test against Xh
	double (*fp)(double,double) = &f;
	Eigen::MatrixXd fh;
	
	std::cout << "Here" << std::endl;

	fh = testXh(g,fp,k,q1d);

	std::cout << "Result for testXh: " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}		

}

double f(double x1, double x2){

	return pow(x1,2) + 3*x2;

}
