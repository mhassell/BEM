//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

class geometry{
    
public:
    
    // methods
    geometry(boost::numeric::ublas::matrix<double>& coord, boost::numeric::ublas::matrix<int>& elts);
    void enhance();
    void refine();
    
    // attributes
    bool enhanced;
    boost::numeric::ublas::matrix<int>& elements;
    boost::numeric::ublas::matrix<double>& coordinates;
    boost::numeric::ublas::matrix<double> normals;
    boost::numeric::ublas::vector<double> lengths;
    boost::numeric::ublas::vector<int> prev;
    boost::numeric::ublas::vector<int> next;
    size_t nElts;
    
};

geometry::geometry(boost::numeric::ublas::matrix<double>& coord, boost::numeric::ublas::matrix<int>& elts)
: xcoords(coordx), ycoords(coordy), elements(elts)
{
    // use the constructor list above to make the needed basic arrays
    enhanced = false;
    nElts = elts.size();
    
    // zero all the other fields before enhancing
    normalx.assign(nElts, 0);
    normaly.assign(nElts, 0);
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


