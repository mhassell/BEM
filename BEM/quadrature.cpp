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

	Eigen::MatrixXd W(M,N);
	
	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < N; j++){
			W(i,j) = f(i,1)*g(j,1);
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

	

}


















