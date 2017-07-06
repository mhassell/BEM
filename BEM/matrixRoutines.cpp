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

#include "matrixRoutines.hpp"

// kronecker product of two (double precision) boost matrices
boost::numeric::ublas::matrix<double> kron(const boost::numeric::ublas::matrix<double>& A, const boost::numeric::ublas::matrix<double>& B)
{
    //kron of an MxN with a PxQ is an MN x PQ matrix
    
    size_t M,N,P,Q, i,j,k,l;
    M = A.size1();
    N = A.size2();
    P = B.size1();
    Q = B.size2();
    
    boost::numeric::ublas::matrix<double> C(M*P,N*Q);
    
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

// kronecker product of a (double precision) sparse with a vanilla boost matrix
boost::numeric::ublas::matrix<double> kron(boost::numeric::ublas::mapped_matrix<double>& A, const boost::numeric::ublas::matrix<double>& B)
{
    //kron of an MxN with a PxQ is an MN x PQ matrix
    
    size_t M,N,P,Q, i,j,k,l;
    M = A.size1();
    N = A.size2();
    P = B.size1();
    Q = B.size2();
    
    boost::numeric::ublas::matrix<double> C(M*P,N*Q);
    
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

// solve Ax = b or AX = B for a vector x or a matrix X
boost::numeric::ublas::matrix<double> solve(const boost::numeric::ublas::matrix<double>& A, const boost::numeric::ublas::matrix<double>& B){
    
    size_t M = B.size1();
    size_t N = B.size2();
    
    boost::numeric::ublas::matrix<double> C(M,N);
    
    
    
    return C;
    
}
