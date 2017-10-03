#ifndef quadrature_hpp
#define quadrature_hpp

#include <Eigen/Dense>

struct preparedQuads{

	Eigen::MatrixXd F2dsing;
	Eigen::MatrixXd F2dssing;

};

Eigen::MatrixXd tensorize(const Eigen::MatrixXd& f, const Eigen::MatrixXd& g);

preparedQuads prepareQuad(const Eigen::MatrixXd &smoothf, const Eigen::MatrixXd &singf);

#endif
