//
//  geometry.hpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#ifndef geometry_hpp
#define geometry_hpp

#include <vector>

class geometry {
    
public:
    // methods
    geometry(std::vector<double>& coordx,
             std::vector<double>& coordy, std::vector<int>& elts);
    
    void enhance();
    void refine();
    
    // attributes
    bool enhanced;
    std::vector<int>& elements;
    std::vector<double>& xcoords;
    std::vector<double>& ycoords;
    std::vector<double> normalx;
    std::vector<double> normaly;
    std::vector<double> lengths;
    std::vector<int> prev;
    std::vector<int> next;
    size_t nElts;
    
};

#endif /* geometry_hpp */
