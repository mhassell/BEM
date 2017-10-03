//
//  matrixRoutines.cpp
//  BEM
//
//  Created by Matthew Hassell on 7/2/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

// a place for routines not in boost/STL or to simplify other built-ins.

// functions:
//     kronecker product of two boost arrays
//     LU solvers (just a wrapper around the boost methods)
//     meshgrid for two Eigen::VectorXd vectors
//     printMatrix (print a matrix nicely)

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "matrixRoutines.hpp"

// kronecker product of two (double precision) Eigen dense matrices
Eigen::MatrixXd kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
{
    //kron of an MxN with a PxQ is an MN x PQ matrix
    
    size_t M,N,P,Q, i,j,k,l;
    M = A.rows();
    N = A.cols();
    P = B.rows();
    Q = B.cols();
    
    Eigen::MatrixXd C = Eigen::MatrixXd(M*P,N*Q);
    
    for(i=0; i<M; i++){
        for(j=0; j<N; j++){
            for(k=0; k<P; k++){
                for(l=0; l<Q; l++){
                    C(P*i+k, Q*j+l)=A(i,j)*B(k,l);
                }
            }
        }
    }
    
    return C;
    
}


// solve AX = B for a matrix X
Eigen::MatrixXd solve(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B){
    
    // size_t M = A.rows();
    size_t N = B.rows();
	size_t P = B.cols();
    
    Eigen::MatrixXd C = Eigen::MatrixXd(N,P);
    
    C = A.partialPivLu().solve(B);
    
    return C;
    
}


// solve Ax = b for a vector x 
Eigen::VectorXd solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b){
    
    size_t M = A.rows();
    
    Eigen::VectorXd C = Eigen::VectorXd(M);
    
    C = A.partialPivLu().solve(b);
    
    return C;
    
}

grid meshgrid(const Eigen::VectorXd &x, const Eigen::VectorXd &y){

	// return copies of x as rows, copies of y as columns

	size_t N = x.rows();
	size_t M = y.rows();

	Eigen::MatrixXd X(M,N);
	Eigen::MatrixXd Y(M,N);

	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < N; j++){
			X(i,j) = x(j);
			Y(i,j) = y(i);
		}
	}	

	grid g = {X,Y};  // should this be new'd to put it on the heap?

	return g;

}


void printMatrix(const Eigen::MatrixXd &A){

	for(size_t i = 0; i < A.rows(); i++){
		for(size_t j = 0; j < A.cols(); j++){
			printf("%5f   ", A(i,j));
		}
		std::cout << std::endl;
	}

}




