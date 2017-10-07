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

	Eigen::MatrixXd q1d = tableGauss(63);
	Eigen::MatrixXd qlog = tableLogGauss(39);

	Eigen::MatrixXd quadf = tensorize(q1d,q1d);
	preparedQuads qds = prepareQuad(q1d, qlog);
	
	Eigen::MatrixXd F2dsing = qds.F2dsing;
	Eigen::MatrixXd F2dssing = qds.F2dssing;

	double (*kernel)(double) = &ker;
	
	Eigen::MatrixXd K = WeaklySingularXh(g, kernel, k, quadf, F2dsing, F2dssing);	
	
	printMatrix(K);
}

double ker(double x){
	return -log(x)/(2*M_PI);
}
