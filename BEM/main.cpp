//
//  main.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <vector>
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/matrix_sparse.hpp>
#include <Eigen/Dense>
#include <math.h>
#include "geometry.hpp"
#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "quadTables.hpp"
#include "matrixRoutines.hpp"
#include "OperatorsAndPotentials.hpp"

int main(){
    
    int k = 0;
	
    Eigen::MatrixXd coords(4,2);
    Eigen::MatrixXd elts(4,2);

	coords << 	0, 	0,
				1, 	0,
			  	0.8, 0.8, 
			  	0.2, 1;

	elts << 	0, 	1,
				2,	3,
				1, 	2,
				3, 	0;

	Eigen::VectorXd test(5);

	test << 1, 2, 3, 4, 5;

	for(size_t i = 0; i < 5; i++){
		std::cout << test(i) << std::endl;
	}
	//geometry g(coords, elts);	    
    
}

 /* 
    std::vector<std::vector<double> > coords;
    std::vector<std::vector<int> > elts;
    coords = {{0,0},{1,0},{0.8,0.8},{0.2,1}};
    elts = {{0,1},{2,3},{1,2},{3,0}};
    geometry g(coords,elts);
    
    std::vector<std::vector<double> > q1d = tableGauss(3);
    
    boost::numeric::ublas::matrix<double> fh((k+1)*g.nElts, (k+1)*g.nElts);
    massMatrixYhYh(g, k, q1d, fh);
	*/
	
    /*
    for(size_t i = 0; i < fh.size1(); i++){
        for(size_t j = 0; j < fh.size2(); j++){
            std::cout << fh(i,j) << "     ";
        }
        std::cout << '\n';
    }
     */
