//
//  legendrebasis.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/18/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "legendrebasis.hpp"

/*
 
 void legendreBasis(int n, std::vector<double> &x, std::vector<double> &y, ...)
 
 Input: 
 
 n  : the degree of Legendre polynomial desired
 &x : reference to a std::vector<double> containing points in (-1,1) to evaluate the polynomials at
 &y : reference to a boost::numeric::ublas::vector<double> containing the values of the polynomials at the points in x
 variadics:
 
 
 
 
 */

void legendreBasis(int n, std::vector<double> &x, boost::numeric::ublas::vector<double> &y, ...){
    
    // number of evaluation points
    size_t nPoints = x.size();
    
    
    
    
    
}
