//
//  testsAndProjections.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/22/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "geometry.hpp"

void testXh(geometry& g, double (*f)(double,double), int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
{
    size_t Nqd = q1d.size();
    size_t Nelts = g.nElts;

    // these  contain the quadrature points mapped to the phyical elements
    std::vector<std::vector<double> > P1t, P2t;
    P1t.assign(Nqd, std::vector<double>(Nelts));
    P2t.assign(Nqd, std::vector<double>(Nelts));
    
    // the function evaluated at the physical quadrature points
    boost::numeric::ublas::matrix<double> F(Nqd, Nelts);
    
    
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t[i][j] = 0.5*(1 - q1d[i][0])*g.coordinates[g.elements[j][0]][0] + 0.5*(1 + q1d[i][0])*g.coordinates[g.elements[j][1]][0];
            P2t[i][j] = 0.5*(1 - q1d[i][0])*g.coordinates[g.elements[j][0]][1] + 0.5*(1 + q1d[i][0])*g.coordinates[g.elements[j][1]][1];
            F(i,j) = f(P1t[i][j], P2t[i][j]);
        }
    }
    
    // polynomial bits
    
    // basis on physical element
    std::vector<double> x;
    x.assign(Nqd, 0);
    for(size_t i = 0; i < Nqd; i++){
        x[i] = q1d[i][0];
    }
    
    boost::numeric::ublas::matrix<double> P(Nqd, Nelts);
    legendreBasis(k, x, 1, P);
    
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P(i,j) = P(i,j)*q1d[i][1];
        }
    }
    
    P = boost::numeric::ublas::trans(P);
    fh = boost::numeric::ublas::prod(P, F);
    
    for(size_t i = 0; i < k+1; i++){
        for(size_t j = 0; j < Nelts; j++){
            fh(i,j) *= 0.5*g.lengths[j];
        }
    }

}

// test against yh for scalar valued functions
void testYh(geometry& g, double (*f)(double,double), int k, std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
    size_t Nqd = q1d.size();
    size_t Nelts = g.nElts;
    
    // these  contain the quadrature points mapped to the phyical elements
    std::vector<std::vector<double> > P1t, P2t;
    P1t.assign(Nqd, std::vector<double>(Nelts));
    P2t.assign(Nqd, std::vector<double>(Nelts));
    
    
    // the function evaluated at the physical quadrature points
    boost::numeric::ublas::matrix<double> F(Nqd, Nelts);
    
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelts; j++){
            P1t[i][j] = 0.5*(1 - q1d[i][0])*g.coordinates[g.elements[j][0]][0] + 0.5*(1 + q1d[i][0])*g.coordinates[g.elements[j][1]][0];
            P2t[i][j] = 0.5*(1 - q1d[i][0])*g.coordinates[g.elements[j][0]][1] + 0.5*(1 + q1d[i][0])*g.coordinates[g.elements[j][1]][1];
            F(i,j) = f(P1t[i][j], P2t[i][j]);
        }
    }
    
    // polynomial bits
    
    // basis on physical element
    std::vector<double> x;
    x.assign(Nqd, 0);
    for(size_t i = 0; i < Nqd; i++){
        x[i] = q1d[i][0];
    }
    
    boost::numeric::ublas::matrix<double> P(Nqd, k+2);
    legendreBasis(k+1, x, 2, P);
    
    // tests
//    for(size_t i = 0; i < Nqd; i++){
//        for(size_t j = 0; j < k+2; j++){
//            std::cout << P(i,j) << "      ";
//        }
//        std::cout << '\n';
//    }
//    
//    std::cout << std::endl;
    
    // attach to quad weights
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < k+1; j++){
            P(i,j) = P(i,j)*q1d[i][1];
        }
    }
    
    P = boost::numeric::ublas::trans(P);
    boost::numeric::ublas::matrix<double> V = boost::numeric::ublas::prod(P, F);
    
    // scale by element lengths
    for(size_t i = 0; i < k+1; i++){
        for(size_t j = 0; j < Nelts; j++){
            V(i,j) *= 0.5*g.lengths[j];
        }
    }
    
//    // tests
//    for(size_t i = 0; i < k+2; i++){
//        for(size_t j = 0; j < Nelts; j++){
//            std::cout << V(i,j) << "     ";
//        }
//        std::cout << std::endl;
//    }
    
    // assemble the nodal DoFs
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][0]) = V(0,j);
    }
    
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][1]) += V(1,j);
    }
    
    // dump what's left in V into fh
    for(size_t j = 0; j < Nelts; j++){
        for(size_t i = 0; i < k+1; i++){
            fh(i,j) = V(i+1,j);
        }
    }
    
}

// test vector valued functions (dotted with n) against Yh
void testYh(geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{

    
    
}

