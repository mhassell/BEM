#include <Eigen/Dense>
#include <iostream>

// a script to figure out 3d reshaping/permuting

int main(){

	Eigen::MatrixXd test(1,27);
	
	for(size_t i = 0; i < test.cols(); i++){
		test(i) = i+1;
	}

	test.resize(3,9);

	for(size_t i = 0 ; i < test.rows(); i++){
		for(size_t j = 0; j < test.cols(); j++){
			std::cout << test(i,j) << " ";
		}	
		std::cout << std::endl;
	}

	double a[3][3][3];

	for(size_t k = 0; k < 3; k++){ // the slowest index will be the index cards
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				a[k][i][j] = test(i,j+3*k);
			}
		} 
	}

	std::cout << std::endl;
	
	// this looks okay
	for(size_t k = 0; k < 3; k++){ 
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				std::cout << a[k][i][j] << " ";
			}
			std::cout << std::endl;
		} 
	std::cout << std::endl;
	std::cout << std::endl;

	}
	
	/*	
	// now to permute the last two indices
	
	double ap[3][3][3];
	for(size_t k = 0; k < 3; k++){
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				ap[j][i][k] = a[k][i][j];
			}
		}	
	}
	// this works
	for(size_t k = 0; k < 3; k++){ 
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				std::cout << ap[k][i][j] << " ";
			}
			std::cout << std::endl;
		} 
	std::cout << std::endl;
	std::cout << std::endl;
	}
	*/

	// now with new and delete!
	const int dim1 = 3;
	const int dim2 = 3;
	const int dim3 = 3;

	/*
	int ***array = new int**[dim1];
	for(size_t i = 0; i < dim1; i++){
		array[i] = new int*[dim2];
		for(size_t j = 0; j < dim2; j++){	
			array[i][j] = new int[dim3];
		}
	}	
	for(size_t k = 0; k < dim1; k++){
		for(size_t i = 0; i < dim2; i++){
			for(size_t j = 0; j < dim3; j++){
				array[k][i][j] = ap[k][i][j];
			}
		}
	}
	for(size_t i = 0; i < dim1; i++){
		for(size_t j = 0; j < dim2; j++){
			for(size_t k = 0; k < dim3; k++){
				std::cout << array[i][j][k] << " ";
			}
		std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	*/

	// and this is what we need

	int ***arrayp = new int**[dim1];

	for(size_t i = 0; i < dim1; i++){
		arrayp[i] = new int*[dim3];
		for(size_t j = 0; j < dim3; j++){	
			arrayp[i][j] = new int[dim2];
		}
	}

	for(size_t i = 0; i < dim1; i++){
		for(size_t j = 0; j < dim2; j++){
			for(size_t k = 0; k < dim3; k++){
				arrayp[k][i][j] = a[k][i][j];
			}
		}
	}	
		
	delete [] arrayp;
	// delete []array;

}
