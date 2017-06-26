//
//  testsAndProjections.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/22/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include "testsAndProjections.hpp"

#include "geometry.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "legendrebasis.hpp"

void testXh(geometry& g, double *f(double,double), int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
{

    // these  contain the quadrature points mapped to the phyical elements
    std::vector<std::vector<double> > P1t, P2t;
    size_t Nqd = q1d.size();
    
    for(size_t i = 0; i < g.nElts; i++){
        for(size_t j = 0; j < Nqd; j++){
            P1t[j][i] = 0.5*(1 - q1d[j][0])*g.coordinates[g.elements[i][0]][0] + 0.5*(1 + q1d[j][0])*g.coordinates[g.elements[i][1]][0];
            P2t[j][i] = 0.5*(1 - q1d[j][0])*g.coordinates[g.elements[i][0]][1] + 0.5*(1 + q1d[j][0])*g.coordinates[g.elements[i][1]][1];
        }
    }
    
    std::cout << std::endl;
    
}


void testYh(geometry& g, double *f(double,double), int k, std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh);
