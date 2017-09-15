#include <Eigen/Dense>

#include "geometry.hpp"
#include "legendrebasis.hpp"

Eigen::MatrixXd testPotentialXh(const geometry& g, double (*kernel)(double,double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d){

	size_t Nelt = g.nElt;
	size_t Nqd = q1d.rows();
	size_t Nobs = obs.rows();

	// these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelts);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelts);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelts);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	P1t.transposeInPlace();
	P2t.transposeInPlace();

	P1t.reshape(1,Nqd*Nelt);
	P2t.reshape(1,Nqd*Nelt);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){
			Z1minuxY1(i,j) = z(i,0) - P1t(0,j);
			Z2minuxY2(i,j) = z(i,1) - P2t(0,j);
		}
	}

	Z1minusY1.reshape(Nobs*Nelt, Nqd);
	Z2minusY2.reshape(Nobs*Nelt, Nqd);	

	Eigen::MatrixXcd SL;

	return SL;

}


