// a script to test solving the Laplace equation

#include <math.h>
#include "geometry.hpp"

double r1(double, double);
double r2(double, double);
double u(double, double);
double v1(double, double);
double v2(double, double);

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
	
	// exact solutions	
	Eigen::VectorXd c0;
	Eigen::VectorXd d0;

	c0 << 0.3, 0.2;
	d0 << 0.4, 0.3;

	
		

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





