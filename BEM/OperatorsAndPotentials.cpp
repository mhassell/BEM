//
//  OperatorsAndPotentials.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

// Includes methods for computing: Mass matrices, BEM matrices, and discretized layer potentials

#include <Eigen/Dense>
#include <iostream>

#include "OperatorsAndPotentials.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"

/*

Eigen::MatrixXd massMatrixXhXh(const geometry& g, int k, const Eigen::MatrixXd& q1d)
{
    
}

void massMatrixXhYh(const geometry& g, int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
}

*/

Eigen::MatrixXd massMatrixYhYh(const geometry& g, int k, const Eigen::MatrixXd& q1d)
{
    
    size_t Nelt = g.nElts;
    size_t Nqd = q1d.rows();
    Eigen::MatrixXd Psi = Eigen::MatrixXd(Nqd, k+2);

    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);   // dump quad nodes into x
    }
    
    // generate the legendre basis and compute the inner products on the reference element
    legendrebasis(k+1, x, 2, Psi);
	
    Eigen::MatrixXd PsiPsi = Eigen::MatrixXd(k+2,k+2);
    Eigen::MatrixXd PsiQd = Eigen::MatrixXd(Nqd, k+2);
    
    // attach quad weights
    for(size_t i = 0; i < PsiQd.rows(); i++){
        for(size_t j = 0; j < PsiQd.cols(); j++){
            PsiQd(i,j) = Psi(i,j)*q1d(i,1);
        }
    }

    // and inner products are done
    PsiQd.transposeInPlace();
    PsiPsi = PsiQd*Psi;
    
    // which is more appropriate? mapped_matrix or compressed_matrix?
    Eigen::MatrixXd lengths = Eigen::MatrixXd(g.nElts,g.nElts);
	lengths.setZero();
    
    // map to physical elements with a sparse matrix of weights and a kronecker product
    for(size_t i = 0; i < g.nElts; i++){
        lengths(i,i) = 0.5*g.lengths(i);
    }
    
    Eigen::MatrixXd M = kron(lengths, PsiPsi);

	// good up to here - 9/6/17
    
    // nodalDOF array
    Eigen::MatrixXd nodalDOF = Eigen::MatrixXd(2,Nelt);
    
    for(size_t i = 0; i < Nelt; i++){
        nodalDOF(0,i) =     (k+2)*i;
        nodalDOF(1,i) = 1 + (k+2)*i;
    }

    // need to reshape this for easier access down below
    Eigen::VectorXd nodalDOFVector = Eigen::VectorXd(2*nodalDOF.cols());

    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < nodalDOF.size(); j++){
            nodalDOFVector(i+j*nodalDOF.size()) = nodalDOF(i,j);
        }
    }

    // internalDOF array
    Eigen::VectorXd internalDOF = Eigen::VectorXd((k+2)*Nelt - 2*nodalDOF.size());
    
	for(size_t i = 0; i < nodalDOF.size(); i++){
		// if the index is not in nodalDOF, then put it into internalDOF
	}
   
    // assemble M in to MM by rows
    Eigen::MatrixXd MM = Eigen::MatrixXd((k+1)*Nelt, (k+2)*Nelt);
    
	for(size_t i = 0; i < Nelt; i++){
        for(size_t j = 0; j < Nelt; j++){
            
        }
    }
    
    //assemble MM into fh by columns

	return MM;   
    
}
