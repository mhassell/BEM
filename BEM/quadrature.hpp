#ifndef quadrature_hpp
#define quadrature_hpp

#include <Eigen/Dense>

struct preparedQuads{

	Eigen::MatrixXd point;
	Eigen::MatrixXd diagonal;

};

Eigen::MatrixXd tensorize(const Eigen::MatrixXd& f, const Eigen::MatrixXd& g);

preparedQuads prepareQuad(const Eigen::MatrixXd &smoothf, const Eigen::MatrixXd &singf);

#endif
