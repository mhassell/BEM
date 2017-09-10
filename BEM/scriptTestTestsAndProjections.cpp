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
	fh = testYh(g,fp1,fp2,k,q1d);	

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
