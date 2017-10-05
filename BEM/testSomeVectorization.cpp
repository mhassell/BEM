#include <iostream>
#include <Eigen/Dense>
#include "matrixRoutines.hpp"

int main(){

	Eigen::MatrixXd column(5,2);
	column << 1, 2,
			  3, 4,
			  5, 6,
			  7, 8,
			  9, 0;

	Eigen::MatrixXd row(2,5);
	row << 1, 2, 3, 4, 5,
		   4, 3, 2, 1, 0;

	Eigen::MatrixXd ans = column.block(0,0,5,1)*row.block(0,0,1,5);

	printMatrix(ans);

}
