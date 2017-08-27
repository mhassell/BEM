//
//  quadTables.hpp
//  BEM
//
//  Created by Matthew Hassell on 6/24/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef quadTables_h
#define quadTables_h

#include <Eigen/Dense>

Eigen::MatrixXd tableGauss(int k);

Eigen::MatrixXd tableLogGauss(int k);

#endif /* quadTables_h */
