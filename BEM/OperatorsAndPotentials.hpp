//
//  OperatorsAndPotentials.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef OperatorsAndPotentials_hpp
#define OperatorsAndPotentials_hpp

#include <Eigen/Dense>

#include "legendrebasis.hpp"
#include "quadTables.hpp"
#include "geometry.hpp"
#include "matrixRoutines.hpp"

// mass matrices
/*
Eigen::MatrixXd massMatrixXhXh(const geometry& g, int k, const Eigen::MatrixXd& q1d);

Eigen::MatrixXd massMatrixXhYh(const geometry& g, int k, const Eigen::MatrixXd& q1d);
*/

Eigen::MatrixXd massMatrixYhYh(const geometry& g, int k, const Eigen::MatrixXd& q1d);

// differentiation matrix

Eigen::MatrixXd differentiationMatrix(const geometry &g, int k);

// operators

// potentials

#endif /* OperatorsAndPotentials_hpp */
