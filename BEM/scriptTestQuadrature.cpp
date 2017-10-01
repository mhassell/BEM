#include <iostream>

#include "quadTables.hpp"
#include "quadrature.hpp"
#include "matrixRoutines.hpp"

int main(){

	Eigen::MatrixXd qdX = tableGauss(5);
	Eigen::MatrixXd qdY = tableGauss(7);
	
	Eigen::MatrixXd qdd = tensorize(qdX, qdY); 

	printMatrix(qdd);

}
