//
//  legendrebasis.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/18/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef legendrebasis_hpp
#define legendrebasis_hpp

#include <Eigen/Dense>

void legendrebasis(int n, Eigen::VectorXd &x, int type, Eigen::MatrixXd &y);

#endif /* legendrebasis_hpp */
