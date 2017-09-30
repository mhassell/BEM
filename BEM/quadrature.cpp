#include <Eigen/Dense>
#include "matrixRoutines.hpp"

Eigen::MatrixXd tensorize(const Eigen::MatrixXd& f, const Eigen::MatrixXd& g){
	
		Eigen::VectorXd xpts = f.slice(0,0,f.rows(),1);
		Eigen::VectorXd ypts = g.slice(0,0,g.rows(),1);
		grid g = meshgrid(xpts,ypts);

		for(size_t i = 0; i < g.xpts.rows(); i++){
			for(size_t j = 0; j < g.xpts.cols(); j++){
				std::cout << g.xpts(i,j) << " ";
			}
			std::cout << std::endl;
		}

}

