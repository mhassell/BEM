//
//  matrixRoutines.hpp
//  BEM
//
//  Created by Matthew Hassell on 7/2/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef matrixRoutines_hpp
#define matrixRoutines_hpp

#include <boost/numeric/ublas/matrix.hpp>

// kronecker product of two (double) matrices
boost::numeric::ublas::matrix<double> kron(const boost::numeric::ublas::matrix<double> A, const boost::numeric::ublas::matrix<double> B);

// solve Ax = b or AX = B for a vector x or a matrix X
boost::numeric::ublas::matrix<double> solve(const boost::numeric::ublas::matrix<double> A, const boost::numeric::ublas::matrix<double> B);

#endif /* matrixRoutines_hpp */
