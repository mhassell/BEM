//
//  testsAndProjections.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/22/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef testsAndProjections_hpp
#define testsAndProjections_hpp

#include "geometry.hpp"
#include <Eigen/Dense>

// test against the BEM spaces

Eigen::MatrixXd testXh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d);

Eigen::MatrixXd testYh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d);

Eigen::MatrixXd testYh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const Eigen::MatrixXd& q1d);

// project into the BEM spaces

Eigen::MatrixXd projectXh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d);

Eigen::MatrixXd projectXh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const Eigen::MatrixXd& q1d);

Eigen::MatrixXd projectYh(const geometry& g, double (*f)(double,double), int k, const Eigen::MatrixXd& q1d);

#endif /* testsAndProjections_hpp */
