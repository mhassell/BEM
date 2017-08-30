//
//  legendrebasis.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/18/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include "legendrebasis.hpp"
#include <tr1/cmath>
#include <Eigen/Dense>

/*
 
 void legendreBasis(int n, Eigen::VectorXd &x, int type, Eigen::MatrixXd &y)
 
 Input: 
 
 n  : the degree of Legendre polynomial desired
 &x : reference to an eigen::MatrixXd containing points in (-1,1) to evaluate the polynomials at
 &y : reference to an eigen::MatrixXd containing the values of the polynomials at the points in x (must be initialized)
type : integer for the type of polynomial we compute

 */

void legendrebasis(int n, Eigen::VectorXd &x, int type, Eigen::MatrixXd &y){
    
    // number of evaluation points
    size_t nPoints = x.size();
    size_t nDegs = n+1;
    
    // type = 0 for orthonormal, type = 1 for orthogonal
    if(type == 0 || type == 1){
        for(size_t i = 0; i < nPoints; i++){
            for(size_t j = 0; j < nDegs; j++){
                y(i,j) = std::tr1::legendre((int) j, x(i));
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
            y(i,0) = (1 - x(i))/2;
            y(i,1) = (1 + x(i))/2;
        }
        if(nDegs == 1){
            return;
        }
        else{
            for(int j = 2; j < nDegs; j++){
                double val = 2*(2* (double) j - 1);
                for(int i = 0; i < nPoints; i++){
                    y(i,j) = (std::tr1::legendre(j,x(i)) - std::tr1::legendre(j-2,x(i)))/sqrt(val);
                }
            }
        }
    }
}
