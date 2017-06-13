//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>
#include <math.h>

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
    
    normals.assign(nElts,std::vector<double>(2));
    
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
    
    std::vector<double> x1(nElts);
    std::vector<double> x2(nElts);
    std::vector<double> y1(nElts);
    std::vector<double> y2(nElts);
    
    double norm;
    
    // fill in all the fields
    for(size_t i = 0; i<nElts; i++){
        x1[i] = coordinates[elements[i][0]][0];
        x2[i] = coordinates[elements[i][1]][0];
        y1[i] = coordinates[elements[i][0]][1];
        y2[i] = coordinates[elements[i][1]][1];
        
        lengths[i] = sqrt(pow(x2[i]-x1[i],2)+pow(y2[i]-y1[i],2));
        
        normals[i][0] = y2[i]-y1[i];
        normals[i][1] = x1[i]-x2[i];
        norm = sqrt(pow(normals[i][0],2)+pow(normals[i][1],2));
        normals[i][0]/= norm;
        normals[i][1]/= norm;
        
        
        
    }
    
    enhanced = true;
    
}

void geometry::refine(){
    
    
    geometry::enhance();
    
}


