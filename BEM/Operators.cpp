#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#include "matrixRoutines.hpp"
#include "geometry.hpp"
#include "legendrebasis.hpp"

Eigen::MatrixXd WeaklySingularXh(const geometry& g, double (*kernel)(double), int k, const Eigen::MatrixXd& quadf, const Eigen::MatrixXd& quads, const Eigen::MatrixXd& quadss){

	size_t Nelt = g.nElts;
	size_t Nqd = quadf.rows();

    Eigen::MatrixXd P1t = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2t = Eigen::MatrixXd(Nqd,Nelt);

	Eigen::MatrixXd P1tau = Eigen::MatrixXd(Nqd,Nelt);
    Eigen::MatrixXd P2tau = Eigen::MatrixXd(Nqd,Nelt);

	for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - quadf(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadf(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - quadf(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadf(i,0))*g.coordinates(g.elements(j,1),1);
            P1tau(i,j) = 0.5*(1 - quadf(i,1))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadf(i,1))*g.coordinates(g.elements(j,1),0);
            P2tau(i,j) = 0.5*(1 - quadf(i,1))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadf(i,1))*g.coordinates(g.elements(j,1),1);
        }
    }

	Eigen::MatrixXd Polt = Eigen::MatrixXd::Zero(Nqd,k+1);
	Eigen::MatrixXd Poltau = Eigen::MatrixXd::Zero(Nqd,k+1);

	Eigen::VectorXd x1(Nqd);
	Eigen::VectorXd x2(Nqd);	
	
	for(size_t i = 0; i < Nqd; i++){
		x1(i) = quadf(i,0);
		x2(i) = quadf(i,1);
	}

	legendrebasis(k,x1,1,Polt);
	legendrebasis(k,x2,1,Poltau);

	Eigen::MatrixXd lengthlength = Eigen::MatrixXd::Zero(Nelt,Nelt);
	for(size_t i = 0; i < Nelt; i++){
		for(size_t j = 0; j < Nelt; j++){
			lengthlength(i,j) = 0.25*g.lengths(i)*g.lengths(j);
		}
	}

	// emulate the neighbor sparse matrix here by hand (csc format)
	for(size_t i = 0; i < Nelt; i++){
		lengthlength(i,i) = 0;
		lengthlength(i,g.next(i)) = 0;
		lengthlength(i,g.prev(i)) = 0;
	}

	//printMatrix(lengthlength);
	
	Eigen::MatrixXd X1minusY1(Nelt,Nelt);
	Eigen::MatrixXd X2minusY2(Nelt,Nelt);
	Eigen::MatrixXd R(Nelt,Nelt);
	Eigen::MatrixXd PolPol = Eigen::MatrixXd::Zero(k+1,k+1);
	Eigen::MatrixXd ker(Nelt,Nelt); 

	Eigen::MatrixXd K = Eigen::MatrixXd::Zero((k+1)*Nelt,(k+1)*Nelt);

	for(size_t q = 0; q < Nqd; q++){		
		PolPol = quadf(q,2)*Polt.block(q,0,1,k+1).transpose()*Poltau.block(q,0,1,k+1);
		// compute pairwise diffs (not vectorized in eigen)
		for(size_t i = 0; i < Nelt; i++){
			for(size_t j = 0; j < Nelt; j++){
				X1minusY1(i,j) = P1t(q,i) - P1tau(q,j);
				X2minusY2(i,j) = P2t(q,i) - P2tau(q,j);
				R(i,j) = pow(X1minusY1(i,j),2) + pow(X2minusY2(i,j),2);
				R(i,j) = pow(R(i,j),0.5);
				ker(i,j) = kernel(R(i,j))*lengthlength(i,j);					
			}
		}
		K += kron(ker, PolPol);
	}

	// clean up the infs/nans
	for(size_t i = 0; i < K.rows(); i++){
		for(size_t j = 0; j < K.cols(); j++){
			if(std::isinf(K(i,j)) || std::isnan(K(i,j))){
				K(i,j) = 0;			
			}		
		}
	}

	// diagonal element interactions
	
	Nqd = quadss.rows();

	Polt.resize(Nqd,k+1);
	Poltau.resize(Nqd,k+1);
	
	x1.resize(Nqd);
	x2.resize(Nqd);
		
	Polt.setZero();
	Poltau.setZero();

	x1.setZero();
	x2.setZero();

	for(size_t i = 0; i < Nqd; i++){
		x1(i) = quadss(i,0);
		x2(i) = quadss(i,1);
	}

	legendrebasis(k,x1,1,Polt);
	legendrebasis(k,x2,1,Poltau);
	
	P1t.resize(Nqd,Nelt);
	P2t.resize(Nqd,Nelt);

	P1t.setZero();
	P2t.setZero();

	P1tau.resize(Nqd,Nelt);
	P2tau.resize(Nqd,Nelt);

	P1tau.setZero();
	P2tau.setZero();

	Eigen::MatrixXd Kdiag(Nqd,Nelt);

	for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
            P1t(i,j) = 0.5*(1 - quadss(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadss(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 - quadss(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadss(i,0))*g.coordinates(g.elements(j,1),1);
            P1tau(i,j) = 0.5*(1 - quadss(i,1))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quadss(i,1))*g.coordinates(g.elements(j,1),0);
            P2tau(i,j) = 0.5*(1 - quadss(i,1))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quadss(i,1))*g.coordinates(g.elements(j,1),1);
        }
    }			

	double r;	
		
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			r = pow(P1t(i,j)-P1tau(i,j),2) + pow(P2t(i,j)-P2tau(i,j),2);	
			Kdiag(i,j) = 0.25*kernel(pow(r,0.5))*pow(g.lengths(j),2);
		}
	}

	Eigen::MatrixXd Kdq = Eigen::MatrixXd::Zero(Nelt,Nelt);

	for(size_t q = 0; q < Nqd; q++){
		PolPol.setZero();		
		PolPol = quadss(q,2)*Polt.block(q,0,1,k+1).transpose()*Poltau.block(q,0,1,k+1);		
		for(size_t i = 0; i < Nelt; i++){
				Kdq(i,i) = Kdiag(q,i);
		}		
		K += kron(Kdq,PolPol);
	}

	// bottom corner singularity
	Nqd = quads.rows();

	Polt.resize(Nqd,k+1);
	Poltau.resize(Nqd,k+1);

	Polt.setZero();
	Poltau.setZero();

	x1.resize(Nqd);
	x2.resize(Nqd);

	x1.setZero();
	x2.setZero();	
		
	for(size_t i = 0; i < Nqd; i++){
		x1(i) = -quads(i,0);
		x2(i) = quads(i,1); 
	} 		

	legendrebasis(k,x1,1,Polt);
	legendrebasis(k,x2,1,Poltau);

	P1t.resize(Nqd,Nelt);
	P2t.resize(Nqd,Nelt);

	P1t.setZero();
	P2t.setZero();

	P1tau.resize(Nqd,Nelt);
	P2tau.resize(Nqd,Nelt);

	P1tau.setZero();
	P2tau.setZero();

	for(size_t i = 0; i < Nqd; i++){
        for(size_t j = 0; j < Nelt; j++){
			// careful of the sign on P1t and P2t!
            P1t(i,j) = 0.5*(1 + quads(i,0))*g.coordinates(g.elements(j,0),0) + 0.5*(1 - quads(i,0))*g.coordinates(g.elements(j,1),0);
            P2t(i,j) = 0.5*(1 + quads(i,0))*g.coordinates(g.elements(j,0),1) + 0.5*(1 - quads(i,0))*g.coordinates(g.elements(j,1),1);
            P1tau(i,j) = 0.5*(1 - quads(i,1))*g.coordinates(g.elements(j,0),0) + 0.5*(1 + quads(i,1))*g.coordinates(g.elements(j,1),0);
            P2tau(i,j) = 0.5*(1 - quads(i,1))*g.coordinates(g.elements(j,0),1) + 0.5*(1 + quads(i,1))*g.coordinates(g.elements(j,1),1);
        }
    }			
		
	Eigen::MatrixXd Knext = Eigen::MatrixXd::Zero(Nqd,Nelt);
	Eigen::MatrixXd Kprev = Eigen::MatrixXd::Zero(Nqd,Nelt);

	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			r = pow(P1t(i,j) - P1tau(i,g.next(j)),2) + pow(P2t(i,j) - P2tau(i,g.next(j)),2);
			Knext(i,j) = kernel(pow(r,0.5));
			r = pow(P1t(i,g.prev(j)) - P1tau(i,j),2) + pow(P2t(i,g.prev(j)) - P2tau(i,j),2);
			Kprev(i,j) = kernel(pow(r,0.5));					
		}
	}

	// why do I have to do this separately?
	for(size_t i = 0; i < Nqd; i++){
		for(size_t j = 0; j < Nelt; j++){
			Knext(i,j) *= 0.25*g.lengths(j)*g.lengths(g.next(j));
			Kprev(i,j) *= 0.25*g.lengths(j)*g.lengths(g.prev(j));
		}
	}

	for(size_t q = 0; q < Nqd; q++){

		PolPol.setZero();
		PolPol = quads(q,2)*Polt.block(q,0,1,k+1).transpose()*Poltau.block(q,0,1,k+1);	
		Kdq.setZero();	

		for(size_t i = 0; i < Nelt; i++){
			Kdq(i,g.next(i)) = Knext(q,i);	
		}

		K += kron(Kdq,PolPol);
		Kdq.setZero();

		for(size_t i = 0; i < Nelt; i++){
			Kdq(i,g.prev(i)) = Kprev(q,i);	
		}
		
		K += kron(Kdq,PolPol.transpose());		
	
	}

	return K;
	
}

Eigen::MatrixXd DipoleXhYh(const geometry& g, double (*kernel)(double), int k, const Eigen::MatrixXd& quadf, const Eigen::MatrixXd& quads){

		

}
