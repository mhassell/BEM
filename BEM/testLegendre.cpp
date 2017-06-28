//
//  testLegendre.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/27/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include "testLegendre.hpp"
#include "legendrebasis.hpp"

// some functions to test the generation of the Legendre basis

void testLegendre(int k, std::vector<double> pts, int type){
    
    // generate legendre basis of the prescribed type and print out the points
    size_t Nqd = pts.size();
    boost::numeric::ublas::matrix<double> y(Nqd, k+1);
    
    legendreBasis(k, pts, type, y);
    
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < k+1; j++){
            std::cout << y(i,j) << "         ";
        }
        std::cout << '\n';
    }
    
}
