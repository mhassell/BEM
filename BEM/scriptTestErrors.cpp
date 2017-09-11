#include <iostream>
#include <Eigen/Dense>
#include <math.h>


#include "geometry.hpp"
#include "testsAndProjections.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"
#include "quadTables.hpp"
#include "OperatorsAndPotentials.hpp"
#include "errors.hpp"

double f(double, double);
double f1(double, double);
double f2(double, double);

int main(){
	
	// levels of refinement
	int nLevels = 5;

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

	// function pointers
	double (*fp)(double,double) = &f;
	double (*fp1)(double,double) = &f1;
	double (*fp2)(double,double) = &f2;

	Eigen::MatrixXd fh;
/*
	// projecting into Xh
	std::cout << "Result for projectXh (scalar case): " << std::endl;
	for(int i = 0; i < nLevels; i++){

		g.refine();		
	
		fh.setZero();
		fh = projectXh(g,fp,k,q1d);

		double errorXhScalar = errorXh(g,fp,fh,k,q1d);
		std::cout << errorXhScalar << std::endl;

	}
	
	// project vector functions into Xh	 	
	fh.setZero();

	std::cout << "Result for projectXh (vector case): " << std::endl;
	for(int i = 0; i < nLevels; i++){

		g.refine();		
	
		fh.setZero();
		fh = projectXh(g,fp1,fp2,k,q1d);

		double errorXhVector = errorXh(g,fp1,fp2,fh,k,q1d);
		std::cout << errorXhVector << std::endl;

	}

*/

	// project scalar functions into Yh
	std::cout << "Result for projectYh (scalar case): " << std::endl;
	for(int i = 0; i < nLevels; i++){

		g.refine();		
	
		fh.setZero();
		fh = projectYh(g,f,k,q1d);

		double errorYhScalar = errorYh(g,f,fh,k,q1d);
		std::cout << errorYhScalar << std::endl;

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
