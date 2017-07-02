//
//  legendrebasis.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/18/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef legendrebasis_hpp
#define legendrebasis_hpp

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>

void legendreBasis(int n, const std::vector<double> &x, int type, boost::numeric::ublas::matrix<double> &y);

#endif /* legendrebasis_hpp */
