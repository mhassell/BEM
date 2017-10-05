#include <Eigen/Dense>
#include <iostream>

#include "matrixRoutines.hpp"
#include "geometry.hpp"
#include "legendrebasis.hpp"

Eigen::MatrixXd WeaklySingularXh(const geometry& g, double (*kernel)(double), int k, const Eigen::MatrixXd& quadf, const Eigen::MatrixXd& quads, const Eigen::MatrixXd& quadss){

	std::cout << "here" << std::endl;


	size_t Nelt = g.nElts;
	size_t Nqd = quads.rows();

    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);

	Eigen::MatrixXd P1tau = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2tau = Eigen::MatrixXd(Nqd,Nelt);

	for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - quadf(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadf(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - quadf(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadf(i,0))*g.coordinates(g.elements(j,1),1);
            P1tau(i,j) = 0.5*(1 - quadf(i,1))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadf(i,1))*g.coordinates(g.elements(j,1),0);
            P2tau(i,j) = 0.5*(1 - quadf(i,1))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadf(i,1))*g.coordinates(g.elements(j,1),1);
        }
    }

	Eigen::MatrixXd Polt = Eigen::MatrixXd::Zero(Nqd,Nelt);
	Eigen::MatrixXd Poltau = Eigen::MatrixXd::Zero(Nqd,Nelt);

	Eigen::VectorXd x1(Nqd);
	Eigen::VectorXd x2(Nqd);	
	
	for(size_t i = 0; i < Nqd; i++){
		x1(i) = quadf(i,0);
		x2(i) = quadf(i,1);
	}

	legendrebasis(k,x1,0,Polt);
	legendrebasis(k,x2,0,Poltau);

	Eigen::MatrixXd lengthlength = Eigen::MatrixXd::Zero(Nelt,Nelt);
	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < Nelt; j++){
			lengthlength(i,j) = 0.25*g.lengths(i)*g.lengths(j);
		}
	}

	// emulate the neighbor sparse matrix here by hand (csc format)
	for(size_t i = 0; i < Nelt; i++){
		lengthlength(i,i) = 0;
		lengthlength(i,g.next(i)) = 0;
		lengthlength(i,g.prev(i)) = 0;
	}

	printMatrix(lengthlength);
	
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero((k+1)*Nelt,(k+1)*Nelt);
	
	
	
}
