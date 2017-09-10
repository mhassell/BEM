#ifndef errors_hpp
#define errors_hpp

#include "geometry.hpp"
#include "quadTables.hpp"
#include "legendrebasis.hpp"

double errorXh(const geometry& g, double (*f)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d);

#endif
