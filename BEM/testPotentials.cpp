#include <Eigen/Dense>
#include <iostream>

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
	/*
	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			tmp(0,i*Nqd + j) = P1t(i,j);
		}
	}
	*/
	//P1t = tmp;

	for(size_t i = 0; i < P1t.rows(); i++){
		for(size_t j = 0 ; j < P1t.cols(); j++){
			std::cout << P1t(i,j) << " ";		
		}
		std::cout << std::endl;
	}
	tmp.setZero();

	P2t.transposeInPlace();
	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			tmp(0,i*Nqd + j) = P2t(i,j);
		}
	}
	P2t = tmp;
	tmp.setZero();

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd*Nelt);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd*Nelt);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){
			
			Z1minusY1(i,j) = obs(i,0) - P1t(0,j);
			Z2minusY2(i,j) = obs(i,1) - P2t(0,j);
		}
	}

	/*
	for(size_t i = 0; i < Z2minusY2.rows(); i++){
		for(size_t j = 0; j < Z2minusY2.cols(); j++){
			std::cout << Z2minusY2(i,j) << " ";
		}
		std::cout << std::endl;
	}
	*/	

	// Nelt x 1 blocks of (Nobs x Nqd)
	tmp.resize(Nobs*Nelt,Nqd);
	tmp.setZero();
	
	Z2minusY2.transposeInPlace();
	Z2minusY2.resize(Nobs*Nelt,Nqd);

	/*
	for(size_t i = 0; i < Z2minusY2.rows(); i++){
		for(size_t j = 0; j < Z2minusY2.cols(); j++){
			std::cout << Z2minusY2(i,j) << " ";
		}
		std::cout << std::endl;
	}
	*/	
	
	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nelt; j++){
			for(size_t k = 0; k < Nqd; k++){
				// std::cout << i << " " << j << " " << k << std::endl;
				tmp(i*Nelt+j,k) = Z1minusY1(i,k*Nelt+j);		
			}
		}	
	}	

	/*
	for(size_t i = 0; i < tmp.rows(); i++){
		for(size_t j = 0; j < tmp.cols(); j++){
			std::cout << tmp(i,j) << " ";
		}
		std::cout << std::endl;
	}	
	*/
	Eigen::MatrixXd SL;

	return SL;

}


