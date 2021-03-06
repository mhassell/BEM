#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#include "geometry.hpp"
#include "legendrebasis.hpp"
#include "matrixRoutines.hpp"

Eigen::MatrixXd testPotentialXh(const geometry& g, double (*kernel)(double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d){

	size_t Nelt = g.nElts;
	size_t Nqd = q1d.rows();
	size_t Nobs = obs.rows();

	// these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelt);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// reshape P1t and P2t
	Eigen::MatrixXd tmp(1,Nqd*Nelt);
	
	P1t.transposeInPlace();   // these appear to be expensive calls
	P1t.resize(1,Nqd*Nelt);
	
	P2t.transposeInPlace();   // these appear to be expensive calls
	P2t.resize(1,Nqd*Nelt);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd*Nelt);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd*Nelt);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){		
			Z1minusY1(i,j) = obs(i,0) - P1t(0,j);
			Z2minusY2(i,j) = obs(i,1) - P2t(0,j);
		}
	}

	Z1minusY1.resize(Nobs*Nelt, Nqd);
	Z2minusY2.resize(Nobs*Nelt, Nqd);
	
	Eigen::MatrixXd r(Nobs*Nelt, Nqd);

	for(size_t i = 0; i < Nobs*Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			r(i,j) = pow(Z1minusY1(i,j),2) + pow(Z2minusY2(i,j),2);
			r(i,j) = pow(r(i,j), 0.5);
		}
	}

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	Eigen::MatrixXd lb(Nqd,k+1);
	legendrebasis(k,x,1,lb);

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < lb.cols(); j++){
			lb(i,j) *= q1d(i,1);
		}
	}

	Eigen::MatrixXd SL(r.rows(), r.cols());

	for(size_t i = 0; i < r.rows(); i++){
		for(size_t j = 0; j < r.cols(); j++){
			SL(i,j) = kernel(r(i,j));		
		}	
	}

	SL = SL*lb;

	Eigen::MatrixXd lens(Nelt,1);

	for(size_t i = 0; i < Nelt; i++){
		lens(i) = 0.5*g.lengths(i);
	}

	Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(Nobs,1,1);
	
	lens = kron(lens,ones);

	for(size_t i = 0; i < SL.rows(); i++){
		for(size_t j = 0; j < SL.cols(); j++){
			SL(i,j) *= lens(i);
		}
	}

	Eigen::MatrixXd tmp2 = Eigen::MatrixXd::Zero(Nobs,(k+1)*Nelt);

	// use SL.block(i,j,m,n) here to simplify greatly 
	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < k+1; j++){
			tmp2.block(0,j+i*(k+1),Nobs,1) = SL.block(i*Nobs,j,Nobs,1);
		}	
	}

	SL = tmp2;

	return SL;

}


Eigen::MatrixXd testPotentialYh(const geometry& g, double (*kernel)(double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d){

	size_t Nelt = g.nElts;
	size_t Nqd = q1d.rows();
	size_t Nobs = obs.rows();

	// these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelt);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// reshape P1t and P2t
	Eigen::MatrixXd tmp(1,Nqd*Nelt);
	
	P1t.transposeInPlace();    // these appear to be expensive calls
	P1t.resize(1,Nqd*Nelt);
	
	P2t.transposeInPlace();    // these appear to be expensive calls
	P2t.resize(1,Nqd*Nelt);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd*Nelt);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd*Nelt);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){		
			Z1minusY1(i,j) = obs(i,0) - P1t(0,j);
			Z2minusY2(i,j) = obs(i,1) - P2t(0,j);
		}
	}

	Eigen::MatrixXd N(Nobs,Nqd*Nelt);
	
	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){
			N(i,j) = Z1minusY1(i,j)*g.normals(j%Nelt,0)+Z2minusY2(i,j)*g.normals(j%Nelt,1);
		}
	}

	Z1minusY1.resize(Nobs*Nelt, Nqd);
	Z2minusY2.resize(Nobs*Nelt, Nqd);
	N.resize(Nobs*Nelt, Nqd);
	
	Eigen::MatrixXd r(Nobs*Nelt, Nqd);

	for(size_t i = 0; i < Nobs*Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			r(i,j) = pow(Z1minusY1(i,j),2) + pow(Z2minusY2(i,j),2);
			r(i,j) = pow(r(i,j), 0.5);
		}
	}

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	Eigen::MatrixXd Psi(Nqd,k+2);
	legendrebasis(k+1,x,2,Psi);

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Psi.cols(); j++){
			Psi(i,j) *= q1d(i,1);
		}
	}

	Eigen::MatrixXd M(r.rows(), r.cols());

	for(size_t i = 0; i < r.rows(); i++){
		for(size_t j = 0; j < r.cols(); j++){
			M(i,j) = kernel(r(i,j))*N(i,j);		
		}	
	}

	M = M*Psi;

	Eigen::MatrixXd lens(Nelt,1);

	for(size_t i = 0; i < Nelt; i++){
		lens(i) = 0.5*g.lengths(i);
	}

	Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(Nobs,1,1);
	
	lens = kron(lens,ones);

	for(size_t i = 0; i < M.rows(); i++){
		for(size_t j = 0; j < M.cols(); j++){
			M(i,j) *= lens(i);
		}
	}

	Eigen::MatrixXd tmp2 = Eigen::MatrixXd::Zero(Nobs,(k+2)*Nelt);

	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < k+2; j++){
			tmp2.block(0,j+i*(k+2),Nobs,1) = M.block(i*Nobs,j,Nobs,1);
		}	
	}

	// assembly

    Eigen::VectorXd nodalDOF1 = Eigen::VectorXd(Nelt);
	Eigen::VectorXd nodalDOF2 = Eigen::VectorXd(Nelt);    

    for(size_t i = 0; i < Nelt; i++){
        nodalDOF1(i) =     (k+2)*i;
        nodalDOF2(i) = 1 + (k+2)*i;
    }

    Eigen::VectorXd internalDOF((k+2)*Nelt - 2*Nelt);
	for(size_t i = 0; i < Nelt; i++){
		// each element has 2 nodal DOFs and k internal DOFs (for a total of k+2 DOFs)
		for(size_t j = 0; j < k; j++){
			internalDOF(i*k + j) = 2+(k+2)*i + j;	
		}		
	}

	// assemble by columns

	Eigen::MatrixXd K(Nobs, (k+1)*Nelt);
	size_t Nnd = g.nElts;
	
	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < tmp2.rows(); j++){
			K(j, g.elements(i,0)) = tmp2(j,nodalDOF1(i)); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < tmp2.rows(); j++){			
			K(j,g.elements(i,1)) += tmp2(j,nodalDOF2(i));   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < tmp2.rows(); j++){
			K(j,Nnd+i) = tmp2(j,internalDOF(i));
		}
	}
	
	return K;

}


Eigen::MatrixXd testPotentialYhSL(const geometry& g, double (*kernel)(double), const Eigen::MatrixXd& obs, int k, const Eigen::MatrixXd& q1d){

	size_t Nelt = g.nElts;
	size_t Nqd = q1d.rows();
	size_t Nobs = obs.rows();

	// these  contain the quadrature points mapped to the phyical elements
    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);
    // the function evaluated at the physical quadrature points
    Eigen::MatrixXd F = Eigen::MatrixXd(Nqd, Nelt);

    for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - q1d(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + q1d(i,0))*g.coordinates(g.elements(j,1),1);
        }
    }

	// reshape P1t and P2t
	Eigen::MatrixXd tmp(1,Nqd*Nelt);
	
	P1t.transposeInPlace();
	P1t.resize(1,Nqd*Nelt);
	
	P2t.transposeInPlace();
	P2t.resize(1,Nqd*Nelt);

	Eigen::MatrixXd Z1minusY1(Nobs,Nqd*Nelt);
	Eigen::MatrixXd Z2minusY2(Nobs,Nqd*Nelt);

	for(size_t i = 0; i < Nobs; i++){
		for(size_t j = 0; j < Nqd*Nelt; j++){		
			Z1minusY1(i,j) = obs(i,0) - P1t(0,j);
			Z2minusY2(i,j) = obs(i,1) - P2t(0,j);
		}
	}

	Z1minusY1.resize(Nobs*Nelt, Nqd);
	Z2minusY2.resize(Nobs*Nelt, Nqd);
	
	Eigen::MatrixXd r(Nobs*Nelt, Nqd);

	for(size_t i = 0; i < Nobs*Nelt; i++){
		for(size_t j = 0; j < Nqd; j++){
			r(i,j) = pow(Z1minusY1(i,j),2) + pow(Z2minusY2(i,j),2);
			r(i,j) = pow(r(i,j), 0.5);
		}
	}

	Eigen::VectorXd x(Nqd);
	for(size_t i = 0; i < Nqd; i++){
		x(i) = q1d(i,0);
	}
	
	Eigen::MatrixXd Psi(Nqd,k+2);
	legendrebasis(k+1,x,2,Psi);

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Psi.cols(); j++){
			Psi(i,j) *= q1d(i,1);
		}
	}

	Eigen::MatrixXd M(r.rows(), r.cols());

	

	for(size_t i = 0; i < r.rows(); i++){
		for(size_t j = 0; j < r.cols(); j++){
			M(i,j) = kernel(r(i,j));		
		}	
	}

	M = M*Psi;

	Eigen::MatrixXd lens(Nelt,1);

	for(size_t i = 0; i < Nelt; i++){
		lens(i) = 0.5*g.lengths(i);
	}

	Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(Nobs,1,1);
	
	lens = kron(lens,ones);

	for(size_t i = 0; i < M.rows(); i++){
		for(size_t j = 0; j < M.cols(); j++){
			M(i,j) *= lens(i);
		}
	}

	Eigen::MatrixXd tmp2 = Eigen::MatrixXd::Zero(Nobs,(k+2)*Nelt);

	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < k+2; j++){
			tmp2.block(0,j+i*(k+2),Nobs,1) = M.block(i*Nobs,j,Nobs,1);
		}	
	}

	// assembly

    Eigen::VectorXd nodalDOF1 = Eigen::VectorXd(Nelt);
	Eigen::VectorXd nodalDOF2 = Eigen::VectorXd(Nelt);    

    for(size_t i = 0; i < Nelt; i++){
        nodalDOF1(i) =     (k+2)*i;
        nodalDOF2(i) = 1 + (k+2)*i;
    }

    Eigen::VectorXd internalDOF((k+2)*Nelt - 2*Nelt);
	for(size_t i = 0; i < Nelt; i++){
		// each element has 2 nodal DOFs and k internal DOFs (for a total of k+2 DOFs)
		for(size_t j = 0; j < k; j++){
			internalDOF(i*k + j) = 2+(k+2)*i + j;	
		}		
	}

	// assemble by columns

	Eigen::MatrixXd K(Nobs, (k+1)*Nelt);
	size_t Nnd = g.nElts;
	
	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < tmp2.rows(); j++){
			K(j, g.elements(i,0)) = tmp2(j,nodalDOF1(i)); 
        }
    }

	for(size_t i = 0; i < Nnd; i++){
        for(size_t j = 0; j < tmp2.rows(); j++){			
			K(j,g.elements(i,1)) += tmp2(j,nodalDOF2(i));   
        }
    }

	for(size_t i = 0; i < internalDOF.rows(); i++){
		for(size_t j = 0; j < tmp2.rows(); j++){
			K(j,Nnd+i) = tmp2(j,internalDOF(i));
		}
	}
	
	return K;

}
