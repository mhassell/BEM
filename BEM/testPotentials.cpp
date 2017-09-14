#include <Eigen/Dense>
#include <complex>

#include "geometry.hpp"
#include "legendrebasis.hpp"

Eigen::MatrixCf testPotentialXh(const geometry& g, std::complex (*kernel)(double,double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d, std::complex s){

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

	P1t.reshape(1,Nqd);
	P2t.reshape(1,Nqd);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd);

	

}
