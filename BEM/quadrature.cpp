#include <Eigen/Dense>
#include <iostream>
#include "matrixRoutines.hpp"
#include "quadrature.hpp"

Eigen::MatrixXd tensorize(const Eigen::MatrixXd& f, const Eigen::MatrixXd& g){
	
	size_t M = f.rows();
	size_t N = g.rows();

	Eigen::VectorXd xpts = f.block(0,0,f.rows(),1);
	Eigen::VectorXd ypts = g.block(0,0,g.rows(),1);
	grid gg = meshgrid(xpts,ypts);	
	
	Eigen::MatrixXd formula(N*M,3);

	gg.X.resize(N*M,1);
	gg.Y.resize(N*M,1);

	formula.block(0,0,N*M,1) = gg.X;
	formula.block(0,1,N*M,1) = gg.Y;

	Eigen::MatrixXd W(N,M);
	
	for(size_t i = 0; i < N; i++){
		for(size_t j = 0; j < M; j++){
			W(i,j) = f(j,1)*g(i,1);
		}
	}

	W.resize(N*M,1);

	formula.block(0,2,N*M,1) = W;

	return formula;

}

preparedQuads prepareQuad(const Eigen::MatrixXd &smoothf, const Eigen::MatrixXd &singf){

	Eigen::MatrixXd GxS = tensorize(smoothf,singf);

	Eigen::VectorXd X = GxS.block(0,0,GxS.rows(),1);
	Eigen::VectorXd Y = GxS.block(0,1,GxS.rows(),1);
	Eigen::VectorXd W = GxS.block(0,2,GxS.rows(),1);

	// point singularity
	Eigen::MatrixXd F2dsing = Eigen::MatrixXd::Zero(2*X.rows(),3);
	
	size_t nPts = X.rows();

	for(size_t i = 0; i < X.rows(); i++){
		F2dsing(i,0) = 2*X(i)*Y(i)-1;
		F2dsing(i,1) = 2*Y(i)-1;
		F2dsing(i,2) = 4*W(i)*Y(i);
	}

	for(size_t i = 0; i < X.rows(); i++){
		F2dsing(i+nPts,0) = 2*Y(i)-1;
		F2dsing(i+nPts,1) = 2*X(i)*Y(i)-1;
		F2dsing(i+nPts,2) = 4*W(i)*Y(i);
	}

	// diagonal singularity
	Eigen::MatrixXd F2dssing = Eigen::MatrixXd::Zero(2*X.rows(),3);

	for(size_t i = 0; i < X.rows(); i++){
		F2dssing(i,0) = 2*X(i)*(1-Y(i))-1;
		F2dssing(i,1) = 2*(Y(i) + X(i)*(1-Y(i)))-1;
		F2dssing(i,2) = 4*(1-Y(i))*W(i);
	}

	for(size_t i = 0; i < X.rows(); i++){
		F2dssing(i+nPts,0) = 2*(Y(i) + X(i)*(1-Y(i)))-1;
		F2dssing(i+nPts,1) = 2*X(i)*(1-Y(i))-1;
		F2dssing(i+nPts,2) = 4*(1-Y(i))*W(i);
	}

	preparedQuads quads = {F2dsing,F2dssing}; //should this be new'd to put it on the heap?

	return quads;

}

allquads allQuadrature(int k, bool overkill){

	if(overkill){
		Eigen::MatrixXd q1d = tableGauss(63);
		Eigen::MatrixXd qsing = tableLogGuass(35);
	}
	else{
		Eigen::MatrixXd q1d = tableGauss(k+2);
		Eigen::MatrixXd qsing = tableLogGuass(k+2);
	}

	Eigen::MatrixXd qreg = q1d;
	
	for(size_t i = 0; i < qreg.rows(); i++){
		
	}
	
	

}

















