//
//  testsAndProjections.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/22/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef testsAndProjections_hpp
#define testsAndProjections_hpp

#include "geometry.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>


// test against the BEM spaces
void testXh(geometry& g, double (*f)(double,double), int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh);
void testYh(geometry& g, double (*f)(double,double), int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh);
//void testYh(geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh);

// project into the BEM spaces
void projectionXh();
void projectionYh();


#endif /* testsAndProjections_hpp */
