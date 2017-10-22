#include<iostream>

#include "quadrature.hpp"
#include "quadTables.hpp"
#include "matrixRoutines.hpp"
#include "Operators.hpp"
#include "geometry.hpp"
#include <math.h>

double ker(double);
double kerDL(double);

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

	allQuads qds = allQuadrature(k,1);

	Eigen::MatrixXd regular = qds.regular;
	Eigen::MatrixXd point = qds.point;
	Eigen::MatrixXd diagonal = qds.diagonal; 
	Eigen::MatrixXd pole = qds.pole;		 
	
	double (*kernel)(double) = &ker;
	double (*kernelDL)(double) = &kerDL;	
	
	// Eigen::MatrixXd V = WeaklySingularXh(g, kernel, k, regular, point, diagonal);	
	// printMatrix(V);

	Eigen::MatrixXd K = DipoleXhYh(g, kernelDL, k, regular, pole);
	printMatrix(K);


}

double ker(double x){
	return -log(x)/(2*M_PI);
}

double kerDL(double x){
	return 1.0/(2*M_PI*pow(x,2));
}
