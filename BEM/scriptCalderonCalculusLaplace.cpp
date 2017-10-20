// a script to test solving the Laplace equation

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sys/time.h>

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
double get_wall_time();

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
	
	// number of refinements
	int Nlev = 2;

	for(size_t i = 0; i < Nlev; i++){
		g.refine();
	}

		double start = get_wall_time();

		Eigen::MatrixXd V = WeaklySingularXh(g, kerSLref, k, regular, point, diagonal);
		Eigen::MatrixXd K = DipoleXhYh(g, kerDLref, k, regular, pole);
		Eigen::MatrixXd D = differentiationMatrix(g,k);
		Eigen::MatrixXd W = D.transpose()*V*D;
		Eigen::MatrixXd M = massMatrixXhYh(g,k,q1d);
		M.transposeInPlace();

		double (*ur)(double,double) = &u;

		Eigen::MatrixXd beta0 = testXh(g, ur, k, q1d);
		beta0.resize((k+1)*g.nElts,1);
	
		double (*v1r)(double,double) = &v1;
		double (*v2r)(double,double) = &v2;

		Eigen::MatrixXd beta1tmp = testYh(g, v1r, v2r, k, q1d);
	
		Eigen::MatrixXd beta1 = Eigen::MatrixXd::Zero((k+1)*g.nElts,1);

		for(size_t i = 0; i < g.nElts; i++){
			beta1(i) = beta1tmp(0,i);
		}

		for(size_t i = 1; i < k; i++){
			for(size_t j = 0; j < g.nElts; j++){
				beta1(i*g.nElts + j) = beta1tmp(i,j);
			}
		}

		Eigen::MatrixXd SL = testPotentialXh(g, kerSLref, z, k, q1d);
		Eigen::MatrixXd DL = testPotentialYh(g, kerDLref, z, k, q1d);

		double (*fonep)(double,double) = &fone;
	
		Eigen::MatrixXd ints = testXh(g, fonep, k, q1d);
		ints.resize(ints.rows()*ints.cols(),1);
		
		Eigen::MatrixXd C = testYh(g, fonep, fonep, k, q1d);

		Eigen::MatrixXd Ct = C.transpose();

		C = kron(C,Ct);

		Eigen::MatrixXd projU = projectYh(g, ur, k, q1d);
	
		Eigen::MatrixXd projV1 = projectXh(g, v1r, k, q1d);
		Eigen::MatrixXd projV2 = projectXh(g, v2r, k, q1d);

		Eigen::MatrixXd projUn = Eigen::MatrixXd::Zero(k+1,g.nElts);

		for(size_t i = 0; i < k+1; i++){
			for(size_t j = 0; j < g.nElts; j++){
				projUn(i,j) = projV1(i,j)*g.normals(j,0) + projV2(i,j)*g.normals(j,1);
			}
		} 

		// First kind indirect Dirichlet

		Eigen::MatrixXd lambda = solve(V,beta0);
		Eigen::MatrixXd uh = SL*lambda;

		double error = 0.0;		
		double diff = 0.0;

		for(size_t i = 0; i < z.rows(); i++){

			diff = std::abs(u(z(i,0),z(i,1))-uh(i));
			error = std::max(error, diff);

		}
		
		std::cout << "First kind indirect Dirichlet error: " << std::endl;
		std::cout << error << std::endl;

		/*
		double end = get_wall_time();
		double time_spent = (end-start);
		std::cout << "Total time: " << time_spent << std::endl;
		*/
		
		// First kind indirect with decaying SL		
		Eigen::MatrixXd Vlam(V.rows()+1, V.cols()+1);
		Vlam.block(0, 0, V.rows(), V.cols()) = V;

		for(size_t i = 0; i < ints.rows(); i++){
			Vlam(i,Vlam.cols()-1) = ints(i,0);
			Vlam(Vlam.rows()-1,i) = ints(i,0);
		}

		Vlam(Vlam.rows()-1, Vlam.cols()-1) = 0;

		// copy beta0 into beta02
		Eigen::MatrixXd beta02(beta0.rows()+1,1);
		for(size_t i = 0; i < beta0.rows(); i++){
			beta02(i) = beta0(i);
		}
		beta02(beta02.rows()-1,0) = 0;
		
		Eigen::MatrixXd lambda2 = solve(Vlam,beta02);
		
		uh = SL*lambda2.block(0,0,beta02.rows()-1,1);

		double cInf = lambda2(lambda2.rows()-1,0);

		diff = 0;
		error = 0;

		for(size_t i = 0; i < uh.rows(); i++){
			uh(i) += cInf;
			diff = std::abs(u(z(i,0),z(i,1))-uh(i));
			error = std::max(error, diff);

		}

		std::cout << "First kind indirect Dirichlet error (with decaying potential): " << std::endl;
		std::cout << error << std::endl;


		// indirect DL Neumann

		std::cout << W.rows() << " " << W.cols() << " " << C.rows() << " " << C.cols() <<  std::endl;		
	
		Eigen::MatrixXd Wprime = W+C;
		std::cout << Wprime.rows() << Wprime.cols() << beta1.rows() << std::endl;
		Eigen::MatrixXd phi = solve(-Wprime,beta1);

		std::cout << "here" << std::endl;
				
		uh = DL*phi;

		diff = 0;
		error = 0;

		for(size_t i = 0; i < uh.rows(); i++){
			diff = std::abs(u(z(i,0),z(i,1))-uh(i));
			error = std::max(error, diff);
		}

		std::cout << "Indirect DL Neumann error: " << std::endl;
		std::cout << error << std::endl;

		// Second kind direct Neumann
		phi.setZero();
		Eigen::MatrixXd RHS = 0.5*M*projUn + K.transpose()*projUn;
		phi = -solve(W+C,RHS);

		uh = DL*phi - SL*projUn;

		diff = 0;
		error = 0;

		for(size_t i = 0; i < uh.rows(); i++){
			diff = std::abs(u(z(i,0),z(i,1))-uh(i));
			error = std::max(error, diff);
		}

		std::cout << "Direct Neumann error: " << std::endl;
		std::cout << error << std::endl;
		

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

	return (x1-0.3)/pow(r1(x1,x2),2) - (x2-0.4)/pow(r2(x1,x2),2);	

}

double v2(double x1, double x2){

	return (x1-0.2)/pow(r1(x1,x2),2) - (x2-0.3)/pow(r2(x1,x2),2);	

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

double get_wall_time(){

    struct timeval time;
	gettimeofday(&time, NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
