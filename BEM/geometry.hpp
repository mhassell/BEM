//
//  geometry.hpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef geometry_hpp
#define geometry_hpp

#include <Eigen/Dense>

class geometry{
    
public:
    
    // methods
    geometry(Eigen::MatrixXd &coords, Eigen::MatrixXd &elts);
    void enhance();
    void refine();
    void refine(Eigen::MatrixXd tag);
    
    // attributes
    bool enhanced;
    Eigen::MatrixXd &elements;
    Eigen::MatrixXd &coordinates;
    Eigen::MatrixXd normals;
    Eigen::MatrixXd lengths;
    Eigen::MatrixXd prev;
    Eigen::MatrixXd next;
    size_t nElts;
    
};

#endif /* geometry_hpp */
