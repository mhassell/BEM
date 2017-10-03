#include <iostream>

#include "quadTables.hpp"
#include "quadrature.hpp"
#include "matrixRoutines.hpp"

int main(){

	Eigen::MatrixXd qdX = tableGauss(5);
	Eigen::MatrixXd qdY = tableGauss(7);
	
	Eigen::MatrixXd qdd = tensorize(qdX, qdY); 

	printMatrix(qdd);

	/*
	Eigen::MatrixXd smoothf = tableGauss(5);
	Eigen::MatrixXd singf = tableLogGauss(7);

	preparedQuads quads = prepareQuad(smoothf, singf);

	std::cout << "Point singularity" << std::endl;
	printMatrix(quads.F2dsing);
	
	std::cout << std::endl;

	std::cout << "Diagonal singularity" << std::endl;
	printMatrix(quads.F2dssing);
	*/

}
