#include <iostream>
#include <Eigen/Dense>
#include "quadTables.hpp"

int main(){

	int k = 2;

	std::cout << "Checking Gaussian quad for k = " << k << std::endl;

	Eigen::MatrixXd q1d = tableGauss(k);

	std::cout << "Number of rows and columns: " << std::endl;

	std::cout << q1d.rows() << "   " << q1d.cols() << std::endl;

	std::cout << "Nodes and weights: " << std::endl;

	for(size_t i = 0; i < q1d.rows(); i++){
		std::cout << q1d(i,0) << "   "  << q1d(i,1) << std::endl;
	}

	std::cout << "\n \n \n" << std::endl;

	std::cout << "Now checking log quadrature" << std::endl;	

	q1d = tableLogGauss(k);

	std::cout << "Number of rows and columns: " << std::endl;

	std::cout << q1d.rows() << "   " << q1d.cols() << std::endl;

	std::cout << "Nodes and weights: " << std::endl;

	for(size_t i = 0; i < q1d.rows(); i++){
		std::cout << q1d(i,0) << "   "  << q1d(i,1) << std::endl;
	}

}
