// a script to time some of the individual functions I've written to figure out where the bottleneck is

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <bench/BenchTimer.h>

#include "geometry.hpp"
#include "quadrature.hpp"
#include "Operators.hpp"
#include "OperatorsAndPotentials.hpp"
#include "testsAndProjections.hpp"
#include "testPotentials.hpp"

double r1(double,double);
double r2(double,double);
double u(double,double);
double v1(double,double);
double v2(double,double);
double kerSL(double);
double kerDL(double);
double fone(double,double);

int main(){


	Eigen::BenchTimer t;
	t.start();	

	// number of refinements
	int Nlev = 7;

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

	for(size_t i = 0; i < Nlev; i++){
		g.refine();
	}

	Eigen::MatrixXd V = WeaklySingularXh(g, kerSLref, k, regular, point, diagonal);
	
	t.stop();
	std::cout << t.value(Eigen::REAL_TIMER) << std::endl;

	return 0;

	Eigen::MatrixXd K = DipoleXhYh(g, kerDLref, k, regular, pole);
	Eigen::MatrixXd D = differentiationMatrix(g,k);
	Eigen::MatrixXd W = D.transpose()*V*D;
	Eigen::MatrixXd M = massMatrixXhYh(g,k,q1d);
	M.transposeInPlace();

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

	return std::log(r1(x1,x2)/r2(x1,x2));

}

double v1(double x1, double x2){

	double r1 = pow(x1-0.3,2) + pow(x2-0.2,2);
	double r2 = pow(x1-0.4,2) + pow(x2-0.3,2);
	double tmp = (x1-0.3)/r1 - (x1-0.4)/r2;
	return tmp;	

}

double v2(double x1, double x2){
	double r1 = pow(x1-0.3,2) + pow(x2-0.2,2);
	double r2 = pow(x1-0.4,2) + pow(x2-0.3,2);
	double tmp = (x2-0.2)/r1 - (x2-0.3)/r2;
	return tmp;	

}

double kerSL(double x){

	return -std::log(x)/(2*M_PI);

}

double kerDL(double x){

	return 1.0/(2*M_PI*pow(x,2));

}

double fone(double x1, double x2){
	
	return 1.0;	

}