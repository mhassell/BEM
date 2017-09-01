#include <tr1/cmath>
#include <iostream>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>

// checking out the cmath legendre functions and comparing with boost libraries

int main(){

	size_t length = 11;

	int degree = 7;

	std::vector<double> boost_vals(length);
	std::vector<double> cmath_vals(length);

	for(size_t i = 0; i < length; i++){

		double point = (double) i/length;
		boost_vals[i] = boost::math::legendre_p(degree, point);
		cmath_vals[i] = std::tr1::legendre(degree, point);

	}	

	std::cout << "Boost value || cmath value" << std::endl;

	for(size_t i = 0; i < length; i++){

		std:: cout << boost_vals[i] << "    " << cmath_vals[i] << std::endl;
		 
	}	

	// testing the sqrt thing

	std::cout << sqrt(8) << std::endl;

}
