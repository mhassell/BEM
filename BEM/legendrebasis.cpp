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
#include <math.h>

#include "legendrebasis.hpp"

/*
 
 void legendreBasis(int n, std::vector<double> &x, boost::numeric::ublas::vector<double> &y, int type)
 
 Input: 
 
 n  : the degree of Legendre polynomial desired
 &x : reference to a std::vector<double> containing points in (-1,1) to evaluate the polynomials at
 &y : reference to a boost::numeric::ublas::vector<double> containing the values of the polynomials at the points in x
type : integer for the type of polynomial we compute

 */

void legendreBasis(int n, std::vector<double> &x, boost::numeric::ublas::matrix<double> &y, int type){
    
    // number of evaluation points
    size_t nPoints = y.size1();
    size_t nDegs = y.size2();
    
    // type = 0 for orthonormal, type = 1 for orthogonal
    if(type == 0 || type == 1){
        for(size_t i = 0; i < nPoints; i++){
            for(size_t j = 0; j < nDegs; j++){
                y(i,j) = boost::math::legendre_p((int) j, x[i]);
                if(type == 0){
                    double val = 2/(2*(double)j+1);
                    y(i,j) = y(i,j)/sqrt(val);
                }
            }
        }
    }
    // lobotto functions
    if(type == 2){
        // make the first two lobotto functions
        for(int i = 0; i < nPoints; i++){
            y(i,0) = (1 - x[i])/2;
            y(i,1) = (1 + x[i])/2;
        }
        if(nDegs == 1){
            return;
        }
        else{
            for(int j = 2; j < nDegs; j++){
                double val = 2*(2* (double) j - 1);
                for(int i = 0; i < nPoints; i++){
                    y(i,j) = (boost::math::legendre_p(j+1,x[i]) - boost::math::legendre_p(j-1,x[i]))/sqrt(val);
                }
            }
        }
    }
}
