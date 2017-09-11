#include <Eigen/Dense>
#include <tr1/cmath>
#include <iostream>

#include "geometry.hpp"
#include "quadTables.hpp"
#include "legendrebasis.hpp"
#include "errors.hpp"

// compute the error between a scalar function and a given array in Xh
double errorXh(const geometry& g, double (*f)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d){

	int Nelt = g.nElts;
	int Nqd = q1d.rows();

	Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// exact solution
	Eigen::MatrixXd F(Nqd,Nelt);
	F.setZero();
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			F(i,j) = f(P1t(i,j),P2t(i,j));
		}
	}

	Eigen::MatrixXd Polt(Nqd, k+1);
	legendrebasis(k, x, 1, Polt);

	if(fh.cols() == 1){
		fh.resize(k+1,Nelt);
	}

	Eigen::MatrixXd fhh = Polt*fh;

	Eigen::MatrixXd tmp(Nqd, Nelt);
	tmp.setZero();
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			tmp(i,j) = pow(F(i,j)-fhh(i,j),2)*q1d(i,1)*g.lengths(j);		
		}
	}

	Eigen::VectorXd tmp2(Nelt);
	tmp2.setZero();

	for(size_t j = 0; j < Nelt; j++){
		for(size_t i = 0; i < Nqd; i++){
			tmp2(j) += tmp(i,j);		
		}
		tmp2(j) *= 0.5;	
	}

	double error = 0;
	
	for(size_t i = 0; i < Nelt; i++){
		error += tmp2(i);
	}

	error = pow(error,0.5);
	
	return error;

}

// compute the error between a vector function (dotted with the normal) and a given array in Xh
double errorXh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d){

	int Nelt = g.nElts;
	int Nqd = q1d.rows();

    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelt);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
            F(i,j) = f1(P1t(i,j),P2t(i,j))*g.normals(j,0) + f2(P1t(i,j),P2t(i,j))*g.normals(j,1);
        }
    }

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	

	Eigen::MatrixXd Polt(Nqd, k+1);
	legendrebasis(k, x, 1, Polt);

	if(fh.cols() == 1){
		fh.resize(k+1,Nelt);
	}

	Eigen::MatrixXd fhh = Polt*fh;

	Eigen::MatrixXd tmp(Nqd, Nelt);
	tmp.setZero();
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			tmp(i,j) = pow(F(i,j)-fhh(i,j),2)*q1d(i,1)*g.lengths(j);		
		}
	}

	Eigen::VectorXd tmp2(Nelt);
	tmp2.setZero();

	for(size_t j = 0; j < Nelt; j++){
		for(size_t i = 0; i < Nqd; i++){
			tmp2(j) += tmp(i,j);		
		}
		tmp2(j) *= 0.5;	
	}

	double error = 0;
	
	for(size_t i = 0; i < Nelt; i++){
		error += tmp2(i);
	}

	error = pow(error,0.5);
	
	return error;

}

// compute the error between a scalar function and a given array in Yh
double errorYh(const geometry& g, double (*f)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d){

	int Nelt = g.nElts;
	int Nqd = q1d.rows();

    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
	
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// exact solution
	Eigen::MatrixXd F(Nqd, Nelt);
	F.setZero();
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			F(i,j) = f(P1t(i,j),P2t(i,j));
		}
	}

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	Eigen::MatrixXd Psi(Nqd, k+2);
	legendrebasis(k, x, 2, Psi);

	if(fh.cols() == 1){
		fh.resize(k+1,Nelt);
	}

	Eigen::MatrixXd Fh(k+2,Nelt);

	
	// disassemble nodal DOFs
	for(size_t i = 0; i < Nelt; i++){
		Fh(0,i) = fh(0,g.elements(i,0));
		Fh(1,i) = fh(0,g.elements(i,1));
	}

	for(size_t i = 2; i < k+2; i++){
		for(size_t j = 0; j < Nelt; j++){
			Fh(i,j) = fh(i,j);		
		}
	}

	Eigen::MatrixXd fhh = Psi*Fh;

	Eigen::MatrixXd tmp(Nqd, Nelt);
	tmp.setZero();

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			tmp(i,j) = pow(F(i,j)-fhh(i,j),2)*q1d(i,1)*g.lengths(j);		
		}
	}

	Eigen::VectorXd tmp2(Nelt);
	tmp2.setZero();

	for(size_t j = 0; j < Nelt; j++){
		for(size_t i = 0; i < Nqd; i++){
			tmp2(j) += tmp(i,j);		
		}
		tmp2(j) *= 0.5;	
	}

	double error = 0;
	
	for(size_t i = 0; i < Nelt; i++){
		error += tmp2(i);
	}

	error = pow(error,0.5);
	
	return error;

}
