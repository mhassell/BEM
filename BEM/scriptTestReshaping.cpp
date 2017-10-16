#include <iostream>
#include <Eigen/Dense>

#include "matrixRoutines.hpp"
#include "testsAndProjections.hpp"

double v1(double x1, double x2){

	return (x1-0.3)/pow(r1(x1,x2),2) - (x2-0.4)/pow(r2(x1,x2),2);	

}

double v2(double x1, double x2){

	return (x1-0.2)/pow(r1(x1,x2),2) - (x2-0.3)/pow(r2(x1,x2),2);	

}

int main(){

	double (*v1r) = &v1;
	double (*v2r) = &v2;

	Eigen::MatrixXd beta0tmp = testYh(g, v1r, v2r, k, q1d);

	

}
