//
//  OperatorsAndPotentials.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

// Includes methods for computing: Mass matrices, BEM matrices, and discretized layer potentials

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "OperatorsAndPotentials.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"

void massMatrixXhXh(const geometry& g, int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
}

void massMatrixXhYh(const geometry& g, int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
}

void massMatrixYhYh(const geometry& g, int k, const std::vector<std::vector<double> >& q1d, boost::numeric::ublas::matrix<double>& fh)
{
    
    size_t Nelt = g.nElts;
    boost::numeric::ublas::matrix<double> Psi(k+2,Nelt);
    
    std::vector<double> x(q1d.size());
    for(size_t i = 0; i < q1d.size(); i++){
        x[i] = q1d[i][0];   // dump quad nodes into x
    }
    
    // generate the legendre basis and compute the inner products on the reference element
    legendreBasis(k+1, x, 2, Psi);
    boost::numeric::ublas::matrix<double> PsiPsi(k+2,k+2);
    boost::numeric::ublas::matrix<double> PsiQd(k+2, Nelt);
    
    // attach quad weights
    for(size_t i = 0; i < PsiQd.size1(); i++){
        for(size_t j = 0; j < PsiQd.size2(); j++){
            PsiQd(i,j) = Psi(i,j)*q1d[i][0];
        }
    }
    
    // and inner products are done
    PsiQd = boost::numeric::ublas::trans(PsiQd);
    PsiPsi = boost::numeric::ublas::prod(PsiQd, Psi);
    
    // which is more appropriate? mapped_matrix or compressed_matrix?
    boost::numeric::ublas::mapped_matrix<double> lengths(g.nElts,g.nElts,g.nElts);
    
    // map to physical elements with a sparse matrix of weights and a kronecker product
    for(size_t i = 0; i < g.nElts; i++){
        lengths(i,i) = g.lengths[i];
    }
    
}
