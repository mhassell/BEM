#include <iostream>
#include <Eigen/Dense>
#include "matrixRoutines.hpp"

int main(){
	
	// first, test the kron
	size_t rowsA = 5;
	size_t colsA = 3;

	size_t rowsB = 2;
	size_t colsB = 2;

	Eigen::MatrixXd A = Eigen::MatrixXd(rowsA,colsA);
	Eigen::MatrixXd B = Eigen::MatrixXd(rowsB,colsB);
	Eigen::MatrixXd C = Eigen::MatrixXd(rowsA*rowsB, colsA*colsB);

	A << 1,2,3,
		 6,7,8,
		 7,6,5,
		 2,1,0,
		 3,4,5;

	B << 2, 1,
		 1, 2;

	C = kron(A,B);
	
	std::cout << "Result of kronecker product" << std::endl;
	for(size_t i = 0; i < C.rows(); i++){
		for(size_t j = 0; j < C.cols(); j++){
			std::cout << C(i,j) << "    ";
		}
		std::cout << std::endl;
	}

	// now, test the solver for matrix-matrix equations
	
	rowsA = 3;
	colsA = 3;

	rowsB = 3;
	colsB = 2;

	A.resize(rowsA,colsA);
	B.resize(rowsB,colsB);

	C.resize(rowsB,colsB);

	// zero them out
	for(size_t i = 0; i < rowsA; i++){
		for(size_t j = 0; j < colsA; j++){
			A(i,j) = 0;
		}		
	}	


	for(size_t i = 0; i < rowsB; i++){
		for(size_t j = 0; j < colsB; j++){
			B(i,j) = 0;
			C(i,j) = 0;
		}		
	}

	A << 1,2,3,
	     5,7,9,
		 1,2,1;

	B << 1,2,
		 3,4,
         5,6;

	C = solve(A,B);

	std::cout << "Solution to the matrix equation:" << std::endl;		
	for(size_t i = 0; i < rowsB; i++){
		for(size_t j = 0; j < colsB; j++){
			std::cout << C(i,j) << "   ";
		}		
		std::cout << std::endl;
	}

	// and finally test the regular solver
	std::cout << "Solution to the linear equation:" << std::endl;
	
	Eigen::VectorXd b = Eigen::VectorXd(3);
	
	b << 1,
		 3,
		 5;		

	Eigen::VectorXd x = solve(A,b);

	for(size_t i = 0; i < b.rows(); i++){
		std::cout << x(i) << std::endl;
	}	
	
}
