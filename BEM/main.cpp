//
//  main.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "geometry.hpp"
#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "quadTables.hpp"
#include "matrixRoutines.hpp"
#include "OperatorsAndPotentials.hpp"

int main(){

	int k = 2;
	
	Eigen::MatrixXd q1d = tableGauss(k);

	std::cout << q1d.rows() << "   " << q1d.cols() << std::endl;

	for(size_t i = 0; i < q1d.rows(); i++){
		std::cout << q1d(i,0) << "   "  << q1d(i,1) << std::endl;
	}

    
}

/*
    
    int k = 0;
	
    Eigen::MatrixXd coords(4,2);
    Eigen::MatrixXi elts(4,2);

	coords << 	0, 	0,
				1, 	0,
			  	0.8, 0.8, 
			  	0.2, 0.6;

	elts << 	0, 	1,
				2,	3,
				1, 	2,
				3, 	0;
	
	geometry g(coords, elts);	  
	
	g.refine();

	for(size_t i = 0; i < 8; i++){
		std::cout << g.coordinates(i,0) << ' ' << g.coordinates(i,1) << std::endl;
	}  

	for(size_t i = 0; i < 8; i++){
		std::cout << g.elements(i,0) << ' ' << g.elements(i,1) << std::endl;
	}

*/
