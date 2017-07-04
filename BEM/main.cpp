//
//  main.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include "geometry.hpp"
#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "quadTables.hpp"
#include "matrixRoutines.hpp"

int main(){
    
    size_t M, N, P, Q;
    M = 2;
    N = 2;
    P = 3;
    Q = 3;
    boost::numeric::ublas::matrix<double> A(M,N);
    boost::numeric::ublas::matrix<double> B(P,Q);
    for(size_t i = 0; i < M; i ++){
        for(size_t j = 0; j < N; j++){
            A(i,j) = 2*j + 3*i;
            std::cout << A(i,j) << "   ";
        }
        std::cout << '\n';
    }
    
    for(size_t i = 0; i < P; i++){
        for(size_t j = 0; j < Q; j++){
            B(i,j) = i*i + j+1;
            std::cout << B(i,j) << "   ";
        }
        std::cout << '\n';
    }
    
    boost::numeric::ublas::matrix<double> C = kron(A,B);
    
    for(size_t i = 0; i < C.size1(); i++){
        for(size_t j = 0; j < C.size2(); j++){
            std::cout << C(i,j) << "        ";
        }
        std::cout << '\n';
    }
    
}
