//
//  OperatorsAndPotentials.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef OperatorsAndPotentials_hpp
#define OperatorsAndPotentials_hpp

#include "legendrebasis.hpp"
#include "quadTables.hpp"
#include "geometry.hpp"
#include <boost/numeric/ublas/matrix.hpp>

// mass matrices
void massMatrixXhXh(geometry& g, int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh);
void massMatrixXhYh(geometry& g, int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh);
void massMatrixYhYh(geometry& g, int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh);

// operators

// potentials

#endif /* OperatorsAndPotentials_hpp */