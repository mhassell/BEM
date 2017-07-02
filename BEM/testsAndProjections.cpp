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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

// test a scalar function against Xh
void testXh(const geometry& g, double (*f)(double,double), int k, const std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
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
    
    boost::numeric::ublas::matrix<double> P(Nqd, k+1);
    legendreBasis(k, x, 1, P);
    
    for(size_t i = 0; i < P.size1(); i++){
        for(size_t j = 0; j < P.size2(); j++){
            P(i,j) = P(i,j)*q1d[i][1];
        }
    }
    
    P = boost::numeric::ublas::trans(P);
    fh = boost::numeric::ublas::prod(P, F);
    
    for(size_t i = 0; i < fh.size1(); i++){
        for(size_t j = 0; j < fh.size2(); j++){
            fh(i,j) *= 0.5*g.lengths[j];
        }
    }

}

// test a scalar valued function against Yh
void testYh(const geometry& g, double (*f)(double,double), int k, const std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
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

    boost::numeric::ublas::matrix<double> Psi(Nqd, k+2);
    legendreBasis(k+1, x, 2, Psi);
    
    // attach to quad weights
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Psi.size2(); j++){
            Psi(i,j) = Psi(i,j)*q1d[i][1];
        }
    }
    
    // transpose and multiply
    Psi = boost::numeric::ublas::trans(Psi);
    boost::numeric::ublas::matrix<double> V = boost::numeric::ublas::prod(Psi, F);
    
    // scale by element lengths
    for(size_t i = 0; i < V.size1(); i++){
        for(size_t j = 0; j < V.size2(); j++){
            V(i,j) *= 0.5*g.lengths[j];
        }
    }
    
    // assemble the nodal DoFs
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][0]) = V(0,j);
    }
    
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][1]) += V(1,j);
    }
    
    // dump what's left in V into fh
    for(size_t i = 1; i < V.size1()-1; i++){
        for(size_t j = 0; j < V.size2(); j++){
            fh(i,j) = V(i+1,j);
        }
    }
    
}

// test vector valued functions (dotted with n) against Yh
void testYh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
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
            F(i,j) = f1(P1t[i][j], P2t[i][j])*g.normals[j][0] + f2(P1t[i][j], P2t[i][j])*g.normals[j][1];
        }
    }
    
    // polynomial bits
    
    // basis on physical element
    std::vector<double> x;
    x.assign(Nqd, 0);
    for(size_t i = 0; i < Nqd; i++){
        x[i] = q1d[i][0];
    }
    
    boost::numeric::ublas::matrix<double> Psi(Nqd, k+2);
    legendreBasis(k+1, x, 2, Psi);
    
    // attach to quad weights
    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Psi.size2(); j++){
            Psi(i,j) = Psi(i,j)*q1d[i][1];
        }
    }
    
    // transpose and multiply
    Psi = boost::numeric::ublas::trans(Psi);
    boost::numeric::ublas::matrix<double> V = boost::numeric::ublas::prod(Psi, F);
    
    // scale by element lengths
    for(size_t i = 0; i < V.size1(); i++){
        for(size_t j = 0; j < V.size2(); j++){
            V(i,j) *= 0.5*g.lengths[j];
        }
    }
    
    // assemble the nodal DoFs
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][0]) = V(0,j);
    }
    
    for(size_t j = 0; j < Nelts; j++){
        fh(0,g.elements[j][1]) += V(1,j);
    }
    
    // dump what's left in V into fh
    for(size_t i = 1; i < V.size1()-1; i++){
        for(size_t j = 0; j < V.size2(); j++){
            fh(i,j) = V(i+1,j);
        }
    }
    
}

// project a scalar onto Xh
void projectXh(const geometry& g, double (*f)(double,double), int k, const std::vector<std::vector<double> > q1d,
               boost::numeric::ublas::matrix<double>& fh)
{
 
    size_t Nqd = q1d.size();
    boost::numeric::ublas::matrix<double> P(Nqd, k+1);
    
    // basis on physical element
    std::vector<double> x;
    x.assign(Nqd, 0);
    for(size_t i = 0; i < Nqd; i++){
        x[i] = q1d[i][0];
    }
    
    legendreBasis(k, x, 1, P);
    
    boost::numeric::ublas::matrix<double> Pt = boost::numeric::ublas::trans(P);
    
    // attach quad weights
    for(size_t i = 0; i < P.size1(); i++){
        for(size_t j = 0; j < P.size2(); j++){
            P(i,j)*=q1d[i][1]/2;
        }
    }
    
    boost::numeric::ublas::matrix<double> PP = boost::numeric::ublas::prod(Pt, P);
    boost::numeric::ublas::matrix<double> PPinv = boost::numeric::ublas::identity_matrix<float>(PP.size1());
    boost::numeric::ublas::permutation_matrix<size_t> pm(PP.size1());
    boost::numeric::ublas::lu_factorize(PP, pm);
    boost::numeric::ublas::lu_substitute(PP, pm, PPinv);
    boost::numeric::ublas::matrix<double> ffh(k+1, g.nElts);
    
    testXh(g, f, k, q1d, ffh);
    
    fh = boost::numeric::ublas::prod(PPinv, ffh);
    
    for(size_t i = 0; i < fh.size1(); i++){
        for(size_t j = 0; j < fh.size2(); j++){
            fh(i,j) /= g.lengths[j];
        }
    }
    
}


// project a vector field along the normal vector into Xh
void projectXh(const geometry& g, double (*f1)(double,double), double(*f2)(double,double), int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{

    size_t Nqd = q1d.size();
    boost::numeric::ublas::matrix<double> P(Nqd, k+1);
    
    // basis on physical element
    std::vector<double> x(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x[i] = q1d[i][0];
    }
    
    legendreBasis(k, x, 1, P);
    
    boost::numeric::ublas::matrix<double> Pt = boost::numeric::ublas::trans(P);
    
    // attach quad weights
    for(size_t i = 0; i < P.size1(); i++){
        for(size_t j = 0; j < P.size2(); j++){
            P(i,j)*=q1d[i][1]/2;
        }
    }
    
    boost::numeric::ublas::matrix<double> PP = boost::numeric::ublas::prod(Pt, P);
    boost::numeric::ublas::matrix<double> PPinv = boost::numeric::ublas::identity_matrix<float>(PP.size1());
    boost::numeric::ublas::permutation_matrix<size_t> pm(PP.size1());
    boost::numeric::ublas::lu_factorize(PP, pm);
    boost::numeric::ublas::lu_substitute(PP, pm, PPinv);
    boost::numeric::ublas::matrix<double> fhx(k+1, g.nElts);
    boost::numeric::ublas::matrix<double> fhy(k+1, g.nElts);
    
    testXh(g, f1, k, q1d, fhx);
    testXh(g, f2, k, q1d, fhy);
    
    fhx = boost::numeric::ublas::prod(PPinv, fhx);
    fhy = boost::numeric::ublas::prod(PPinv, fhy);
    
    for(size_t i = 0; i < fh.size1(); i++){
        for(size_t j = 0; j < fh.size2(); j++){
            fh(i,j) = (g.normals[j][0]/g.lengths[j])*fhx(i,j)
                + (g.normals[j][1]/g.lengths[j])*fhy(i,j);
        }
    }
    
}

// vanilla Yh projection for a scalar function
void projectYh(const geometry& g, double (*f)(double,double), int k, const std::vector<std::vector<double> > q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
    // TBD: needs Yh x Yh mass matrix
    
}





