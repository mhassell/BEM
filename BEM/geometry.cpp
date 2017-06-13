//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>

class geometry{
    
public:
    
    // methods
    geometry(std::vector<std::vector<double> >& coords, std::vector<std::vector<int> >& elts);
    void enhance();
    void refine();
    
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

geometry::geometry(std::vector<std::vector<double> >& coords, std::vector<std::vector<int> >& elts)
: coordinates(coords), elements(elts)
{
    // use the constructor list above to make the needed basic arrays
    enhanced = false;
    nElts = elts.size();
    
    // zero all the other fields before enhancing
    for(size_t i = 0; i < nElts; i++){
        normals[i].assign(2, 0);
    }
    
    lengths.assign(nElts, 0);
    prev.assign(nElts, 0);
    next.assign(nElts, 0);
    
    geometry::enhance();
    
}


void geometry::enhance(){
    
    
    enhanced = true;
    
}

void geometry::refine(){
    
    
    geometry::enhance();
    
}


