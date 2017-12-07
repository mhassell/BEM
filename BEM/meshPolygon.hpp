#ifndef MESH_POLYGON_HPP
#define MESH_POLYGON_HPP

class mesh{

public:

	mesh(const geometry& g);
	void meshPolygon(double* box, double h);

private:

	void rayCasting();

};

#endif
