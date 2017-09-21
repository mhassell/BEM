#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#include "geometry.hpp"
#include "legendrebasis.hpp"

Eigen::MatrixXd testPotentialXh(const geometry& g, double (*kernel)(double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d){

	size_t Nelt = g.nElts;
	size_t Nqd = q1d.rows();
	size_t Nobs = obs.rows();

	// these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelt);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// reshape P1t and P2t
	Eigen::MatrixXd tmp(1,Nqd*Nelt);
	
	P1t.transposeInPlace();
	P1t.resize(1,Nqd*Nelt);
	
	P2t.transposeInPlace();
	P2t.resize(1,Nqd*Nelt);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd*Nelt);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd*Nelt);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){		
			Z1minusY1(i,j) = obs(i,0) - P1t(0,j);
			Z2minusY2(i,j) = obs(i,1) - P2t(0,j);
		}
	}

	Z1minusY1.resize(Nobs*Nelt,Nqd);
	Z2minusY2.resize(Nobs*Nelt,Nqd);
	
	Eigen::MatrixXd r(Nobs*Nelt, Nqd);

	for(size_t i = 0; i < Nobs*Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			r(i,j) = pow(Z1minusY1(i,j),2) + pow(Z2minusY2(i,j),2);
			r(i,j) = pow(r(i,j), 0.5);
		}
	}

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	std::cout << "here" << std::endl;
	
	Eigen::MatrixXd lb;
	legendrebasis(k,x,1,lb);

	std::cout << "here" << std::endl;

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < lb.cols(); j++){
			lb(i,j) *= q1d(i,1);
		}
	}

	Eigen::MatrixXd SL;

	for(size_t i = 0; i < r.rows(); i++){
		for(size_t j = 0; j < r.cols(); j++){
			SL(i,j) = kernel(r(i,j));		
		}	
	}

	SL = SL*lb;

	for(size_t i = 0; i < SL.rows(); i++){
		for(size_t j = 0; j < SL.cols(); j++){
			std::cout << SL(i,j) << " ";
		}	
		std::cout << std::endl;
	}	

	return SL;

}

