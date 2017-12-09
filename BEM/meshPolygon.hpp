#ifndef MESH_POLYGON_HPP
#define MESH_POLYGON_HPP

class mesh{

public:

	mesh(const geometry& g);
	void meshPolygon(double* box, double h, int nx, int ny);

private:

	void rayCasting();
	Eigen::MatrixXi elements;
    Eigen::MatrixXd coordinates;
    Eigen::MatrixXd normals;

};

#endif
