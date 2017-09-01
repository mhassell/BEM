#include <Eigen/Dense>
#include "legendrebasis.hpp"
#include <iostream>

// testing the legendre basis functions

int main(){
	
	int nPoints = 10;
	int degree = 2;

	Eigen::VectorXd points = Eigen::VectorXd(nPoints);
	
	// make evenly spaced points in [-1,1]
	double h = (double) 2/(nPoints-1);	
	
	for(size_t i = 0; i < nPoints; i++){		
		points(i) = -1 + i*h;
	}  

	// make the array to hold the output
	Eigen::MatrixXd y = Eigen::MatrixXd(nPoints, degree+1);	
	
	// first, test the orthonormal polynomials (opt = 0)
	legendrebasis(degree, points, 0, y);
	
	std::cout << "Orthonormal Legendre polynomials" << std::endl;
	for(size_t i = 0; i < nPoints; i++){
		for(size_t j = 0; j < degree+1; j++){
			std::cout << y(i,j) << "          ";		
		}
		std::cout << std::endl;
	}

	// Now check the orthogonal ones
	for(size_t i = 0; i < y.rows(); i++){
		for(size_t j = 0; j < y.cols(); j++){
			y(i,j) = 0;
		}
	}

	legendrebasis(degree, points, 1, y);

	std::cout << "Orthogonal Legendre polynomials" << std::endl;
	for(size_t i = 0; i < nPoints; i++){
		for(size_t j = 0; j < degree+1; j++){
			std::cout << y(i,j) << "          ";		
		}
		std::cout << std::endl;
	}

	// and finally the lobatto functions
	for(size_t i = 0; i < y.rows(); i++){
		for(size_t j = 0; j < y.cols(); j++){
			y(i,j) = 0;
		}
	}
	
	legendrebasis(degree, points, 2, y);

	std::cout << "Lobatto polynomials" << std::endl;
	for(size_t i = 0; i < nPoints; i++){
		for(size_t j = 0; j < degree+1; j++){
			std::cout << y(i,j) << "          ";		
		}
		std::cout << std::endl;
	}

	
}
