#include <iostream>
#include <Eigen/Dense>
#include <math.h>


#include "geometry.hpp"
#include "testsAndProjections.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"
#include "quadTables.hpp"
#include "OperatorsAndPotentials.hpp"

double f(double, double);
double f1(double, double);
double f2(double, double);
double v1(double, double);
double v2(double, double);

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
	g.refine();
	
	// quadrature
	Eigen::MatrixXd q1d = tableGauss(63);	  

	// test against Xh
	double (*fp)(double,double) = &f;
	double (*fp1)(double,double) = &f1;
	double (*fp2)(double,double) = &f2;
	double (*v1r)(double,double) = &v1;
	double (*v2r)(double,double) = &v2;

	Eigen::MatrixXd fh;
	fh = testXh(g,fp,k,q1d);

	std::cout << "Result for testXh: " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}

	// testing against Yh
	fh.setZero();
	fh = testYh(g,fp,k,q1d);

	std::cout << "Result for testYh: " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}	 	

	// projecting into Xh
	fh.setZero();
	fh = projectXh(g,fp,k,q1d);

	std::cout << "Result for projectXh (scalar case): " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}

	// test vector functions against Yh	 	
	fh.setZero();
	fh = testYh(g,v1r,v2r,k,q1d);	

	std::cout << "Result for testYh (vector case): " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}

	// project vector functions into Xh	 	
	fh.setZero();
	fh = projectXh(g,fp1,fp2,k,q1d);	

	std::cout << "Result for projectXh (vector case): " << std::endl;
	for(size_t i = 0; i < k+1; i++){
		for(size_t j = 0; j < g.nElts; j++){
			std::cout << fh(i,j) << "    ";
		}
		std::cout << std::endl;
	}

	// project scalar functions into Yh
	fh.setZero();
	fh = projectYh(g,f,k,q1d);
	
	std::cout << "Result for projectYh (scalar case): " << std::endl;
	for(size_t i = 0; i < fh.size(); i++){
			std::cout << fh(i) << std::endl;
	}
	

}

double f(double x1, double x2){

	return pow(x1,2) + 3*x2;

}

double f1(double x1, double x2){

	return sin(x1) + cos(x2);

}

double f2(double x1, double x2){

	return x1*x2;

}

double r1(double x1,double x2){

	double tmp = pow(x1-0.3,2) + pow(x2-0.2,2);
	return pow(tmp,0.5);

}

double r2(double x1,double x2){

	double tmp = pow(x1-0.4,2) + pow(x2-0.3,2);
	return pow(tmp,0.5);

}

double v1(double x1, double x2){

	double tmp = (x1-0.3)/pow(r1(x1,x2),2);
	tmp -= (x2-0.4)/pow(r2(x1,x2),2);
	return tmp;	

}

double v2(double x1, double x2){

	double tmp = (x1-0.2)/pow(r1(x1,x2),2);
	tmp -= (x2-0.3)/pow(r2(x1,x2),2);
	return tmp;	

}
