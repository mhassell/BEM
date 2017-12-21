// a script to plot solutions to the Laplace equation

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <fstream>

#include "geometry.hpp"
#include "quadrature.hpp"
#include "Operators.hpp"
#include "OperatorsAndPotentials.hpp"
#include "testsAndProjections.hpp"
#include "testPotentials.hpp"
#include "meshPolygon.hpp"

double r1(double,double);
double r2(double,double);
double u(double,double);
double v1(double,double);
double v2(double,double);
double kerSL(double);
double kerDL(double);
double fone(double,double);

// source for the exact solution
double a = 0.3;
double b = 0.4;
double c = 0.2;

int main(){

	// number of refinements
	int Nlev = 5;

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
	mesh myMesh(g);
	double box[4] = {-1, 1, -1, 1};
	double h = 0.1;
	int nx = 500;
	int ny = 500;	

	myMesh.meshPolygon(box, h, nx, ny);
	int meshSize = myMesh.Xpts.size();

	// this should be pulled into the classs
	assert(myMesh.Xpts.size()==myMesh.Ypts.size());
	
	Eigen::MatrixXd z(meshSize,2);

	for(int i = 0; i < meshSize; i++){
		z(i,0) = myMesh.Xpts[i];
		z(i,1) = myMesh.Ypts[i];
	}

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
		
	Eigen::MatrixXd Cprov = testYh(g, fonep, k, q1d);
	Eigen::MatrixXd C((k+1)*g.nElts,1);

	// put the first row of Cprov in the first chunk of C
	for(size_t i = 0; i < g.nElts; i++){
		C(i,0) = Cprov(0,i);
	}

	// now resize the rest of Cprov into C
	for(size_t i = 1; i < k+1; i++){
		C.block(i*g.nElts,0,g.nElts,1) = Cprov.block(i,0,1,g.nElts).transpose();
	}

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

	// flip the sign on beta0 for an exterior problem
	/*	
	for(int i = 0; i < beta0.size(); i++){
		beta0(i) *= -1;
	}
	*/

	projUn.resize(projUn.rows()*projUn.cols(),1);
	Eigen::MatrixXd lambda = solve(V,beta0);
	Eigen::MatrixXd uh = SL*lambda;

	// write the solution to a CSV
	std::ofstream file("solution.csv");
	for(int i = 0; i < meshSize; i++){
		file << myMesh.Xpts[i] << "," << myMesh.Ypts[i] << "," << myMesh.Inside[i]  
				<< "," << uh(i) << std::endl;
	}	
	file.close();	
		
}

double r1(double x1,double x2){

	double tmp = pow(x1-a,2) + pow(x2-c,2);
	return pow(tmp,0.5);

}

double r2(double x1,double x2){

	double tmp = pow(x1-b,2) + pow(x2-c,2);
	return pow(tmp,0.5);

}

double u(double x1, double x2){

	return std::log(r1(x1,x2)/r2(x1,x2));

}

double v1(double x1, double x2){

	double r1 = pow(x1-a,2) + pow(x2-c,2);
	double r2 = pow(x1-b,2) + pow(x2-a,2);
	double tmp = (x1-a)/r1 - (x1-b)/r2;
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

