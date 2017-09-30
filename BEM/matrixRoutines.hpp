//
//  matrixRoutines.hpp
//  BEM
//
//  Created by Matthew Hassell on 7/2/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef matrixRoutines_hpp
#define matrixRoutines_hpp

#include<Eigen/Dense>

struct grid{
	Eigen::MatrixXd X,Y;
};

// kronecker product of two (double precision) matrices
Eigen::MatrixXd kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

// solve Ax = b or AX = B for a vector x or a matrix X
Eigen::MatrixXd solve(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

Eigen::VectorXd solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

grid meshgrid(const Eigen::VectorXd &x, const Eigen::VectorXd &y);

#endif /* matrixRoutines_hpp */
