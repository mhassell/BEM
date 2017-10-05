#ifndef operators_hpp
#define operators_hpp

#include <Eigen/Dense>

#include "geometry.hpp"

Eigen::MatrixXd WeaklySingularXh(const geometry&, double (*kernel)(double), int k, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&);

#endif
