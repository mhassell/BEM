#include <tr1/cmath>
#include <iostream>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>

// checking out the cmath legendre functions and comparing with boost libraries

int main(){

	size_t length = 10;

	int degree = 3;

	std::vector<double> boost_vals(length);
	std::vector<double> cmath_vals(length);

	for(size_t i = 0; i < length; i++){
		boost_vals[i] = boost::math::legendre_p(degree, (double) i/length);
		cmath_vals[i] = std::tr1::assoc_legendre(degree, 0, (double) i/length);
	}	

	std::cout << "Boost value || cmath value" << std::endl;

	for(size_t i = 0; i < length; i++){
		std:: cout << boost_vals[i] << "    " << cmath_vals[i] << std::endl;
		 
	}	


}
