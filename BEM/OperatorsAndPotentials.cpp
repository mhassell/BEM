//
//  OperatorsAndPotentials.cpp
//  BEM
//
//  Created by Matthew Hassell on 6/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

// Includes methods for computing: Mass matrices, BEM matrices, and discretized layer potentials

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "OperatorsAndPotentials.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"
#include "geometry.hpp"
#include "quadTables.hpp"

// XhXh identity
Eigen::MatrixXd massMatrixXhXh(const geometry& g, int k, const Eigen::MatrixXd& q1d)
{

	size_t Nelt = g.nElts;
	size_t Nqd = q1d.rows();

	Eigen::MatrixXd P(Nqd,k+1);
	P.setZero();	

	Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);
    }
	
	legendrebasis(k, x, 1, P);
	
	Eigen::MatrixXd PQd(Nqd,k+1);
	for(size_t i = 0; i < P.rows(); i++){
		for(size_t j = 0; j < P.cols(); j++){
			PQd(i,j) = P(i,j)*q1d(i,1);
		}
	}

	P.transposeInPlace();
	Eigen::MatrixXd PP = P*PQd;

	Eigen::MatrixXd lengths(Nelt,Nelt);
	lengths.setZero();	

	for(size_t i = 0; i < Nelt; i++){
        lengths(i,i) = 0.5*g.lengths(i);
    }
 
    Eigen::MatrixXd M = kron(lengths, PP);

	return M;
    
}

//YhXh Identity
Eigen::MatrixXd massMatrixXhYh(const geometry& g, int k, const Eigen::MatrixXd& q1d)
{

	size_t Nelt = g.nElts;
	size_t Nnd = g.coordinates.rows();
    size_t Nqd = q1d.rows();
    	
    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);   // dump quad nodes into x
    }
    
    // generate the legendre basis and compute the inner products on the reference element

	Eigen::MatrixXd Psi(Nqd, k+2);
	Eigen::MatrixXd P(Nqd,k+1);
	    
	legendrebasis(k+1, x, 2, Psi);
	legendrebasis(k, x, 1, P);	

    Eigen::MatrixXd PsiQd = Eigen::MatrixXd(Nqd, k+2);
    
    // attach quad weights
    for(size_t i = 0; i < PsiQd.rows(); i++){
        for(size_t j = 0; j < PsiQd.cols(); j++){
            PsiQd(i,j) = Psi(i,j)*q1d(i,1);
        }
    }

    // and inner products are done
    PsiQd.transposeInPlace();
    Eigen::MatrixXd PsiP = PsiQd*P;
    
    // this should be sparse!
    Eigen::MatrixXd lengths = Eigen::MatrixXd(Nelt,Nelt);
	lengths.setZero();
    
    // map to physical elements with a sparse matrix of weights and a kronecker product
    for(size_t i = 0; i < g.nElts; i++){
        lengths(i,i) = 0.5*g.lengths(i);
    }
    
    Eigen::MatrixXd M = kron(lengths, PsiP);
    
    // nodalDOF arrays
    Eigen::VectorXd nodalDOF1 = Eigen::VectorXd(Nelt);
	Eigen::VectorXd nodalDOF2 = Eigen::VectorXd(Nelt);    

    for(size_t i = 0; i < Nelt; i++){
        nodalDOF1(i) =     (k+2)*i;
        nodalDOF2(i) = 1 + (k+2)*i;
    }

    // internalDOF array
    Eigen::VectorXd internalDOF((k+2)*Nelt - 2*Nelt);
	for(size_t i = 0; i < Nelt; i++){
		// each element has 2 nodal DOFs and k internal DOFs (for a total of k+2 DOFs)
		for(size_t j = 0; j < k; j++){
			internalDOF(i*k + j) = 2+(k+2)*i + j;	
		}		
	}

    // assemble M in to MM by rows
    Eigen::MatrixXd MM((k+1)*Nelt, (k+1)*Nelt);
	MM.setZero();
	
	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < MM.cols(); j++){
			MM(g.elements(i,0),j) = M(nodalDOF1(i),j); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < MM.cols(); j++){
			MM(g.elements(i,1),j) += M(nodalDOF2(i),j);   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < MM.cols(); j++){
			MM(Nnd+i,j) = M(internalDOF(i),j);
		}
	}

	return MM;    
    
}

// YhYh identity
Eigen::MatrixXd massMatrixYhYh(const geometry& g, int k, const Eigen::MatrixXd& q1d)
{
    
    size_t Nelt = g.nElts;
	size_t Nnd = g.coordinates.rows();
    size_t Nqd = q1d.rows();
    Eigen::MatrixXd Psi(Nqd, k+2);

    Eigen::VectorXd x = Eigen::VectorXd(Nqd);
    for(size_t i = 0; i < Nqd; i++){
        x(i) = q1d(i,0);   // dump quad nodes into x
    }
    
    // generate the legendre basis and compute the inner products on the reference element
    legendrebasis(k+1, x, 2, Psi);
	
    Eigen::MatrixXd PsiPsi = Eigen::MatrixXd(k+2,k+2);
    Eigen::MatrixXd PsiQd = Eigen::MatrixXd(Nqd, k+2);
    
    // attach quad weights
    for(size_t i = 0; i < PsiQd.rows(); i++){
        for(size_t j = 0; j < PsiQd.cols(); j++){
            PsiQd(i,j) = Psi(i,j)*q1d(i,1);
        }
    }

    // and inner products are done
    PsiQd.transposeInPlace();
    PsiPsi = PsiQd*Psi;
    
    // this should be sparse!
    Eigen::MatrixXd lengths = Eigen::MatrixXd(g.nElts,g.nElts);
	lengths.setZero();
    
    // map to physical elements with a sparse matrix of weights and a kronecker product
    for(size_t i = 0; i < g.nElts; i++){
        lengths(i,i) = 0.5*g.lengths(i);
    }
    
    Eigen::MatrixXd M = kron(lengths, PsiPsi);
    
    // nodalDOF arrays
    Eigen::VectorXd nodalDOF1 = Eigen::VectorXd(Nelt);
	Eigen::VectorXd nodalDOF2 = Eigen::VectorXd(Nelt);    

    for(size_t i = 0; i < Nelt; i++){
        nodalDOF1(i) =     (k+2)*i;
        nodalDOF2(i) = 1 + (k+2)*i;
    }

    // internalDOF array
    Eigen::VectorXd internalDOF((k+2)*Nelt - 2*Nelt);
	for(size_t i = 0; i < Nelt; i++){
		// each element has 2 nodal DOFs and k internal DOFs (for a total of k+2 DOFs)
		for(size_t j = 0; j < k; j++){
			internalDOF(i*k + j) = 2+(k+2)*i + j;	
		}		
	}

    // assemble M in to MM by rows
    Eigen::MatrixXd MM((k+1)*Nelt, (k+2)*Nelt);
	MM.setZero();
	
	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < MM.cols(); j++){
			MM(g.elements(i,0),j) = M(nodalDOF1(i),j); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < MM.cols(); j++){
			MM(g.elements(i,1),j) += M(nodalDOF2(i),j);   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < MM.cols(); j++){
			MM(Nnd+i,j) = M(internalDOF(i),j);
		}
	}
	
	// assemble MM into M3 by columns
	Eigen::MatrixXd M3((k+1)*Nelt, (k+1)*Nelt);
	M3.setZero();

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < M3.cols(); j++){
			M3(j,g.elements(i,0)) = MM(j,nodalDOF1(i)); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < M3.cols(); j++){
			M3(j,g.elements(i,1)) += MM(j,nodalDOF2(i));   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < M3.rows(); j++){
			M3(j,Nnd+i) = MM(j,internalDOF(i));
		}
	}

	return M3;   
    
}

// differentiation matrix
Eigen::MatrixXd differentiationMatrix(const geometry &g, int k){

	size_t Nelt = g.nElts;
	size_t Nnd = g.coordinates.rows();

	Eigen::MatrixXd D1(k+1,k+2);
	D1(0,0) = -0.5;
	D1(0,1) = 0.5;
	
	// it would be nice to do this with D1.bottomRightCorner()
	for(size_t i = 1; i < k+1; i++){
		for(size_t j = 2; j < k+2; j++){
			if(i == j - 1){
				double val = (double) 2*i+1;
				D1(i,j) = pow(val/2,0.5);
			}	
		}
	}

	Eigen::MatrixXd lengths(Nelt,Nelt);
	for(size_t i = 0; i < Nelt; i++){
		lengths(i,i) = 2/g.lengths(i);
	}

	Eigen::MatrixXd D2 = kron(lengths,D1);

	// nodalDOF arrays
    Eigen::VectorXd nodalDOF1 = Eigen::VectorXd(Nelt);
	Eigen::VectorXd nodalDOF2 = Eigen::VectorXd(Nelt);    

    for(size_t i = 0; i < Nelt; i++){
        nodalDOF1(i) =     (k+2)*i;
        nodalDOF2(i) = 1 + (k+2)*i;
    }

    // internalDOF array
    Eigen::VectorXd internalDOF((k+2)*Nelt - 2*Nelt);
	for(size_t i = 0; i < Nelt; i++){
		// each element has 2 nodal DOFs and k internal DOFs (for a total of k+2 DOFs)
		for(size_t j = 0; j < k; j++){
			internalDOF(i*k + j) = 2+(k+2)*i + j;	
		}		
	}

	// assemble
	Eigen::MatrixXd D((k+1)*Nelt, (k+1)*Nelt);
	D.setZero();

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < D.cols(); j++){
			D(j,g.elements(i,0)) = D2(j,nodalDOF1(i)); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < D.cols(); j++){
			D(j,g.elements(i,1)) += D2(j,nodalDOF2(i));   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < D.rows(); j++){
			D(j,Nnd+i) = D2(j,internalDOF(i));
		}
	}

	return D;
}


