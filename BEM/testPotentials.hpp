#ifndef testPotentials_hpp
#define testPotentials_hpp

#include <Eigen/Dense>

#include "geometry.hpp"
#include "quadTables.hpp"

Eigen::MatrixXd testPotentialXh(const geometry& g, double (*kernel)(double,double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d);

#endif
