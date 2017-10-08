#ifndef quadrature_hpp
#define quadrature_hpp

#include <Eigen/Dense>

struct preparedQuads{

	Eigen::MatrixXd F2dsing;
	Eigen::MatrixXd F2dssing;

};

struct allQuads{

	Eigen::MatrixXd q1d;
	Eigen::MatrixXd regular;
	Eigen::MatrixXd point;
	Eigen::MatrixXd diagonal;
	Eigen::MatrixXd pole;

};

Eigen::MatrixXd tensorize(const Eigen::MatrixXd& f, const Eigen::MatrixXd& g);

preparedQuads prepareQuad(const Eigen::MatrixXd &smoothf, const Eigen::MatrixXd &singf);

allQuads allQuadrature(int k, bool overkill);

#endif
