#include <iostream>
#include <Eigen/Dense>

#include "OperatorsAndPotentials.hpp"
#include "geometry.hpp"
#include "quadTables.hpp"

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
	Eigen::MatrixXd q1d = tableGauss(63);	  

	Eigen::MatrixXd MYhYh = massMatrixYhYh(g, k, q1d);

}
