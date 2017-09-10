// need to test the geometry functions, seem to have a problem

#include <iostream>
#include "geometry.hpp"

int main(){

	// make the geometry
	Eigen::MatrixXd coords(4,2);
    Eigen::MatrixXi elts(4,2);

	coords << 	0, 	0,
				1, 	0,
			  	0.8, 0.8, 
			  	0.2, 0.6;

	elts << 	0, 	1,
				2,	3,
				1, 	2,
				3, 	0;
	
	geometry g(coords, elts);
	g.refine();

	std::cout << "Coordinates" << std::endl;
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.coordinates(i,0) << "  " << g.coordinates(i,1) << std::endl;
	}
	
	std::cout << "Elements" << std::endl;
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.elements(i,0) << "  " << g.elements(i,1) << std::endl;
	}

	std::cout << "Normal vectors:" << std::endl;  
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.normals(i,0) << "  " << g.normals(i,1) << std::endl;
	}

	std::cout << "Next element: " << std::endl;
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.next(i) << std::endl;
	}

	std::cout << "Prev element: " << std::endl;
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.prev(i) << std::endl;
	}

	std::cout << "Lengths" << std::endl;
	for(size_t i = 0; i < g.nElts; i++){
		std::cout << g.lengths(i) << std::endl;
	}
	
}
