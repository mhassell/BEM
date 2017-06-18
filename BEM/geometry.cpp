//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <vector>
#include <math.h>
#include <iostream>

class geometry{
    
public:
    
    // methods
    geometry(std::vector<std::vector<double> >& coords, std::vector<std::vector<int> >& elts);
    void enhance();
    void refine();
    void refine(std::vector<int> tags);
    
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
    geometry::enhance();
    
}


void geometry::enhance(){
    
    nElts = elements.size();
    
    normals.assign(nElts,std::vector<double>(2));
    lengths.assign(nElts, 0);
    prev.assign(nElts, 0);
    next.assign(nElts, 0);
    
    std::vector<double> x1(nElts);
    std::vector<double> x2(nElts);
    std::vector<double> y1(nElts);
    std::vector<double> y2(nElts);
    
    std::vector<int> tmp(nElts);
    
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
        normals[i][0] /= norm;
        normals[i][1] /= norm;
        
        tmp[elements[i][0]] = (int) i;
        
    }
    
    // it would be great to be able to do this in one felswoop
    for(size_t i = 0; i<nElts; i++){
        next[i] = tmp[elements[i][1]];
        prev[next[i]] = (int) i;
    }
    
    enhanced = true;
    
}

// uniform refinement of all elements
void geometry::refine(){
    
    std::vector<std::vector<double> > allCoord(2*nElts);
    allCoord.assign(2*nElts,std::vector<double>(2));
    std::vector<std::vector<int> > allElts(2*nElts);
    allElts.assign(2*nElts, std::vector<int>(2));
    
    for(size_t i = 0; i < nElts; i++){
        // put the new coordinates at the beginning of the vector
        allCoord[i][0] = coordinates[i][0];
        allCoord[i][1] = coordinates[i][1];
        
        // and put the new coordinates at the end of the vector
        allCoord[nElts+i][0] =
            0.5*(coordinates[elements[i][0]][0] + coordinates[elements[i][1]][0]);
        allCoord[nElts+i][1] =
            0.5*(coordinates[elements[i][0]][1] + coordinates[elements[i][1]][1]);
        
        allElts[2*i+1][0] = (int) (nElts+i);
        allElts[2*i][0] = (int) elements[i][0];
        
        allElts[2*i][1] = (int) (nElts+i);
        allElts[2*i+1][1] = (int) elements[i][1];

    }
    
    coordinates = allCoord;
    elements = allElts;
    
    geometry::enhance();
    
}

// bisect only the tagged elements
void geometry::refine(std::vector<int> tag){
    
    // new elements
    std::vector<std::vector<double> > newElts(tag.size());
    
    geometry::enhance();
    
}


