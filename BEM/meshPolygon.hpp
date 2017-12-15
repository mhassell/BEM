#ifndef MESH_POLYGON_HPP
#define MESH_POLYGON_HPP

#include <Eigen/Dense>
#include <vector>

#include "geometry.hpp"

class mesh{

public:

	mesh(const geometry& g);
	void meshPolygon(double* box, double h, int nx, int ny);
	void saveMesh();
	std::vector<double> Xpts;
	std::vector<double> Ypts;

private:

	Eigen::MatrixXi elements;
    Eigen::MatrixXd coordinates;
    Eigen::MatrixXd normals;

};

#endif
