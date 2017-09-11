#ifndef errors_hpp
#define errors_hpp

#include "geometry.hpp"
#include "quadTables.hpp"
#include "legendrebasis.hpp"

double errorXh(const geometry& g, double (*f)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d);

double errorXh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d);

double errorYh(const geometry& g, double (*f)(double,double), Eigen::MatrixXd fh, int k, const Eigen::MatrixXd& q1d);

#endif
