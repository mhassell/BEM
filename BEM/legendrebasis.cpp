//
//  legendrebasis.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/18/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include "legendrebasis.hpp"

/*
 
 void legendreBasis(int n, std::vector<double> &x, boost::numeric::ublas::vector<double> &y, int type)
 
 Input: 
 
 n  : the degree of Legendre polynomial desired
 &x : reference to a std::vector<double> containing points in (-1,1) to evaluate the polynomials at
 &y : reference to a boost::numeric::ublas::vector<double> containing the values of the polynomials at the points in x
case : integer for the type of polynomial we compute

 */

void legendreBasis(int n, std::vector<double> &x, boost::numeric::ublas::matrix<double> &y, int type){
    
    // number of evaluation points
    size_t nPoints = y.size1();
    size_t nDegs = y.size2();
    
    // non-normalized legendre functions
    if(type == 0){
        for(size_t i = 0; i < nPoints; i++){
            for(size_t j = 0; j < nDegs; j++){
                y(i,j) = boost::math::legendre_p((int) j, x[i]);
            }
        }
    }
    
    
    
}
