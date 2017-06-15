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

class geometry{
    
public:
    
    // methods
    geometry(std::vector<std::vector<double> >& coords, std::vector<std::vector<int> >& elts);
    void enhance();
    void refine();
    void refine(std::vector<int> tag);
    
    // attributes
    bool enhanced;
    std::vector<std::vector<int> >& elements;
    std::vector<std::vector<double> >& coordinates;
    std::vector<std::vector<double> > normals;
    std::vector<double> lengths;
    std::vector<int> prev;
    std::vector<int> next;
    size_t nElts;
    
};

#endif /* geometry_hpp */
