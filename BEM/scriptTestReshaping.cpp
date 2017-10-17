// a script to test solving the Laplace equation

#include <Eigen/Dense>
#include <math.h>

#include "geometry.hpp"
#include "quadrature.hpp"
#include "testsAndProjections.hpp"

double r1(double, double);
double r2(double, double);
double u(double, double);
double v1(double, double);
double v2(double, double);
double kerSL(double);
double kerDL(double);
double fone(double);

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

	// number of refinements
	int Nlev = 7;
	
	// observation points 
	Eigen::MatrixXd z(4,2);

	z << -1, -1, 
		 -1,  0,
		 -1,  1,
		 -1,  2;
	
	// quadrature things 
	allQuads quads = allQuadrature(k,1);
	Eigen::MatrixXd q1d = quads.q1d;
	Eigen::MatrixXd regular = quads.regular;
	Eigen::MatrixXd point = quads.point;
	Eigen::MatrixXd diagonal = quads.diagonal;
	Eigen::MatrixXd pole = quads.pole;
		
	double (*kerSLref)(double) = &kerSL;
	double (*kerDLref)(double) = &kerDL;

	double (*ur)(double, double) = &u;
	
	double (*v1r)(double, double) = &v1;
	double (*v2r)(double, double) = &v2;

	Eigen::MatrixXd beta0tmp = testYh(g, v1r, v2r, k, q1d);
	
	Eigen::MatrixXd beta0 = Eigen::MatrixXd::Zero((k+1)*g.nElts,1);

	for(size_t i = 0; i < g.nElts; i++){
		beta0(i) = beta0tmp(0,i);
	}

	beta0tmp.block(1,0,k,g.nElts).resize(k*g.nElts,1); // no good

	printMatrix(beta0tmp);

}

double r1(double x1,double x2){

	double tmp = pow(x1-0.3,2) + pow(x2-0.2,2);
	return pow(tmp,0.5);

}

double r2(double x1,double x2){

	double tmp = pow(x1-0.4,2) + pow(x2-0.3,2);
	return pow(tmp,0.5);

}

double u(double x1, double x2){

	return log(r1(x1,x2)/r2(x1,x2));

}

double v1(double x1, double x2){

	return (x1-0.3)/pow(r1(x1,x2),2) - (x2-0.4)/pow(r2(x1,x2),2);	

}

double v2(double x1, double x2){

	return (x1-0.2)/pow(r1(x1,x2),2) - (x2-0.3)/pow(r2(x1,x2),2);	

}

double kerSL(double x){

	return -log(x)/(2*M_PI);

}

double kerDL(double x){

	return 1.0/(2*M_PI*pow(x,2));

}

double fone(double x){
	
	return 1.0;	

}
