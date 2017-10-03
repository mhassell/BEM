#include <iostream>
#include <Eigen/Dense>

int main(){

	struct grid{
		Eigen::MatrixXd X,Y;	
	};

	Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,3);
	Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(5,7);
	

	grid g = {X, Y};

	std::cout << g.X.rows() << std::endl;
	std::cout << g.Y.rows() << std::endl;

	for(size_t i = 0; i < X.rows(); i++){
		for(size_t j = 0; j < X.cols(); j++){
			g.X(i,j) = i+j;  // assigns to g by <it> value </it>
		}	
	}

	for(size_t i = 0; i < X.rows(); i++){
		for(size_t j = 0; j < X.cols(); j++){
			std::cout << g.X(i,j) << "  ";
		}	
		std::cout << std::endl;
	}

}
