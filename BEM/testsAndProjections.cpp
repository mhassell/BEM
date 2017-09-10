//
//  testsAndProjections.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/22/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <Eigen/Dense>

#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "geometry.hpp"
#include "matrixRoutines.hpp"
#include "OperatorsAndPotentials.hpp"

// test a scalar function against Xh
Eigen::MatrixXd testXh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d)
{
    size_t Nqd = q1d.rows();
    size_t Nelts = g.nElts;
	
	Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(k+1,Nelts);

    // these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelts);
	Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelts);
	
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelts);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
            F(i,j) = f(P1t(i,j), P2t(i,j));
        }
    }

    // polynomial bits
    
    // basis on physical element
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);    
	for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }


    Eigen::MatrixXd P = Eigen::MatrixXd(Nqd, k+1);
    legendrebasis(k, x, 1, P);

    for(size_t i = 0; i < P.rows(); i++){
        for(size_t j = 0; j < P.cols(); j++){
            P(i,j) = P(i,j)*q1d(i,1);
        }
    }
    
    P.transposeInPlace();   // don't do P = P.transpose()
    fh = P*F;				// nice, overloaded operators!


    for(size_t i = 0; i < fh.rows(); i++){
        for(size_t j = 0; j < fh.cols(); j++){
			fh(i,j) = fh(i,j)*0.5*g.lengths(j);
        }
    }

	return fh;

}

// test a scalar valued function against Yh
Eigen::MatrixXd testYh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d)
{
    
    size_t Nqd = q1d.rows();
    size_t Nelts = g.nElts;

	Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(k+1,Nelts);
    
    // these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelts);
	Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelts);
	
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelts);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
            F(i,j) = f(P1t(i,j), P2t(i,j));
        }
    }
    
    // polynomial bits
    
    // basis on physical element
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }

    Eigen::MatrixXd Psi = Eigen::MatrixXd(Nqd, k+2);
    legendrebasis(k+1, x, 2, Psi);
    
    // attach to quad weights
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Psi.cols(); j++){
            Psi(i,j) = Psi(i,j)*q1d(i,1);
        }
    }
    
    // transpose and multiply
    Psi.transposeInPlace();
    Eigen::MatrixXd V = Psi*F;
    
    // scale by element lengths
    for(size_t i = 0; i < V.rows(); i++){
        for(size_t j = 0; j < V.cols(); j++){
            V(i,j) *= 0.5*g.lengths(j);
        }
    }
    
    // assemble the nodal DoFs
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements(j,0)) = V(0,j);
    }
    
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements(j,1)) += V(1,j);
    }
    
    // dump what's left in V into fh
    for(size_t i = 1; i < V.rows()-1; i++){
        for(size_t j = 0; j < V.cols(); j++){
            fh(i,j) = V(i+1,j);
        }
    }

	return fh;
    
}

// test vector valued functions (dotted with n) against Yh
Eigen::MatrixXd testYh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const Eigen::MatrixXd& q1d)
{
    
    size_t Nqd = q1d.rows();
    size_t Nelts = g.nElts;

	Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(k+1,Nelts);
    
    // these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelts);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelts);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelts);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
            F(i,j) = f1(P1t(i,j),P2t(i,j))*g.normals(j,0) + f2(P1t(i,j),P2t(i,j))*g.normals(j,1);
        }
    }


    // polynomial bits
    
    // basis on physical element
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }    

    Eigen::MatrixXd Psi = Eigen::MatrixXd(Nqd, k+2);
    legendrebasis(k+1, x, 2, Psi);
    
    // attach to quad weights
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Psi.cols(); j++){
            Psi(i,j) = Psi(i,j)*q1d(i,1);
        }
    }
    
    // transpose and multiply
    Psi.transposeInPlace();
    Eigen::MatrixXd V = Psi*F;
    
    // scale by element lengths
    for(size_t i = 0; i < V.rows(); i++){
        for(size_t j = 0; j < V.cols(); j++){
            V(i,j) *= 0.5*g.lengths(j);
        }
    }
    
    // assemble the nodal DoFs
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements(j,0)) = V(0,j);
    }
    
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements(j,1)) += V(1,j);
    }
    
    // dump what's left in V into fh
    for(size_t i = 1; i < V.rows()-1; i++){
        for(size_t j = 0; j < V.cols(); j++){
            fh(i,j) = V(i+1,j);
        }
    }

	return fh;
    
}

// project a scalar onto Xh
Eigen::MatrixXd projectXh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d)
{
 
    size_t Nqd = q1d.rows();
	size_t Nelts = g.nElts;    

	Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(k+1,Nelts);

    // basis on physical element
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }

	Eigen::MatrixXd P = Eigen::MatrixXd(Nqd, k+1);
    
    legendrebasis(k, x, 1, P);
    
    Eigen::MatrixXd Pt = P.transpose();
    
    // attach quad weights
    for(size_t i = 0; i < P.rows(); i++){
        for(size_t j = 0; j < P.cols(); j++){
            P(i,j)*=q1d(i,1)/2;
        }
    }

	// got up to here
    
    Eigen::MatrixXd PP = Pt*P;
    Eigen::MatrixXd ffh = Eigen::MatrixXd(k+1, g.nElts);
    
    ffh = testXh(g, f, k, q1d);
    
    fh = solve(PP, ffh);
   
    for(size_t i = 0; i < fh.rows(); i++){
        for(size_t j = 0; j < fh.cols(); j++){
            fh(i,j) /= g.lengths(j);
        }
    }

	return fh;
    
}


// project a vector field along the normal vector into Xh
Eigen::MatrixXd projectXh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const Eigen::MatrixXd& q1d)
{

    size_t Nqd = q1d.rows();
	size_t Nelts = g.nElts;

	Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(k+1, Nelts);

    Eigen::MatrixXd P = Eigen::MatrixXd(Nqd, k+1);
    
    // basis on physical element
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }
    
    legendrebasis(k, x, 1, P);
    
    Eigen::MatrixXd Pt = P.transpose();
    
    // attach quad weights
    for(size_t i = 0; i < P.rows(); i++){
        for(size_t j = 0; j < P.cols(); j++){
            P(i,j)*=q1d(i,1)/2;
        }
    }
    
    Eigen::MatrixXd PP = Pt*P;
    Eigen::MatrixXd fhx = Eigen::MatrixXd(k+1, g.nElts);
    Eigen::MatrixXd fhy = Eigen::MatrixXd(k+1, g.nElts);
    
    fhx = testXh(g, f1, k, q1d);
    fhy = testXh(g, f2, k, q1d);
    
    fhx = solve(PP, fhx);
    fhy = solve(PP, fhy);
    
    for(size_t i = 0; i < fh.rows(); i++){
        for(size_t j = 0; j < fh.cols(); j++){
            fh(i,j) = (g.normals(j,0)/g.lengths(j))*fhx(i,j)
                + (g.normals(j,1)/g.lengths(j))*fhy(i,j);
        }
    }

	return fh;
    
}

// vanilla Yh projection for a scalar function
Eigen::MatrixXd projectYh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d)
{
    
	size_t Nelt = g.nElts;
	size_t Nnd = g.nElts;

    Eigen::MatrixXd M = massMatrixYhYh(g,k,q1d);
	Eigen::MatrixXd b = testYh(g,f,k,q1d);

	Eigen::VectorXd rhs((k+1)*Nelt);
	rhs.setZero();
	
	// nodal DOFs first
	for(size_t i = 0; i < Nnd; i++){
		rhs(i) = b(0,i);	
	}

	// now internal DOFs, by element
    for(size_t j = 1; j < k+1; j++){
		for(size_t i = 0; i < Nelt; i++){
			rhs(Nelt+k*i+j-1) = b(j,i);		// careful with the indexing!
		}
	}

	Eigen::MatrixXd fh = solve(M,rhs);

	return fh;
	
}
