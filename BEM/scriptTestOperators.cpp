#include<iostream>

#include "quadrature.hpp"
#include "quadTables.hpp"
#include "matrixRoutines.hpp"
#include "Operators.hpp"
#include "geometry.hpp"
#include <math.h>

double ker(double);

int main(){

	// polynomial degree
	int k = 2;

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

	allQuads qds = allQuadrature(k,0);

	Eigen::MatrixXd regular = qds.regular;
	Eigen::MatrixXd point = qds.point;
	Eigen::MatrixXd diagonal = qds.diagonal; 		 
	
	double (*kernel)(double) = &ker;
	
	Eigen::MatrixXd K = WeaklySingularXh(g, kernel, k, regular, point, diagonal);	
	
	printMatrix(K);
	
}

double ker(double x){
	return -log(x)/(2*M_PI);
}
