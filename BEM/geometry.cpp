//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <Eigen/Dense>
#include <math.h>

class geometry{
    
public:
    
    // methods
    geometry(Eigen::MatrixXd &coords, Eigen::MatrixXi &elts);
    void enhance();
    void refine();
    void refine(Eigen::MatrixXi tags);
    
    // attributes
    bool enhanced;
    Eigen::MatrixXi & elements;
    Eigen::MatrixXd & coordinates;
    Eigen::MatrixXd normals;
    Eigen::MatrixXd lengths;
    Eigen::MatrixXi prev;
    Eigen::MatrixXi next;
    size_t nElts;
    
};

geometry::geometry(Eigen::MatrixXd &coords, Eigen::MatrixXi &elts): coordinates(coords), elements(elts)
{
    // use the constructor list above to make the needed basic arrays
    enhanced = false;
    size_t nElts = this->elements.rows();
    geometry::enhance();
    
}


void geometry::enhance(){
    
    nElts = elements.rows();
    
    this->normals = Eigen::MatrixXd::Zero(nElts, 2);
    this->lengths = Eigen::VectorXd::Zero(nElts);
    this->prev = Eigen::VectorXi::Zero(nElts);
    this->next = Eigen::VectorXi::Zero(nElts);
    
    Eigen::VectorXd x1 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd x2 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd y1 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd y2 = Eigen::VectorXd::Zero(nElts);
    
    Eigen::VectorXi tmp = Eigen::VectorXi::Zero(nElts);
    
    double norm;
    
    // fill in all the fields
    for(size_t i = 0; i<nElts; i++){
        x1(i) = coordinates(elements(i,0), 0);
        x2(i) = coordinates(elements(i,1), 0);
        y1(i) = coordinates(elements(i,0), 1);
        y2(i) = coordinates(elements(i,1), 1);
  
        lengths(i) = sqrt(pow(x2(i)-x1(i),2)+pow(y2(i)-y1(i),2));
        
        this->normals(i,0) = y2(i)-y1(i);
        this->normals(i,1) = x1(i)-x2(i);
        norm = sqrt(pow(normals(i,0),2)+pow(normals(i,1),2));
        this->normals(i,0) /= norm;
        this->normals(i,1) /= norm;

		tmp(elements(i,0)) = (int) i;
		
    }

    // it would be great to be able to do this in one felswoop
    for(size_t i = 0; i<nElts; i++){
        this->next(i) = tmp(elements(i,1));       
		this->prev(next(i)) = (int) i;
    }
    
    enhanced = true;
    
}

// uniform refinement of all elements
void geometry::refine(){
    
    Eigen::MatrixXd allCoord = Eigen::MatrixXd::Zero(2*nElts,2);
  	Eigen::MatrixXi allElts = Eigen::MatrixXi::Zero(2*nElts,2);
    
    for(size_t i = 0; i < nElts; i++){
        // put the new coordinates at the beginning of the vector
        allCoord(i,0) = coordinates(i,0);
        allCoord(i,1) = coordinates(i,1);
        
        // and put the new coordinates at the end of the vector
        allCoord(nElts+i,0) =
            0.5*(coordinates(elements(i,0),0) + coordinates(elements(i,1),0));
        allCoord(nElts+i,1) =
            0.5*(coordinates(elements(i,0),1) + coordinates(elements(i,1),1));
        
        allElts(2*i+1,0) = (int) (nElts+i);
        allElts(2*i,0) = (int) elements(i,0);
        
        allElts(2*i,1) = (int) (nElts+i);
        allElts(2*i+1,1) = (int) elements(i,1);

    }
    
    this->coordinates = allCoord;
    this->elements = allElts;
    
    geometry::enhance();
    
}

/*
// bisect only the tagged elements
void geometry::refine(std::vector<int> tag){
    
    // new elements
    std::vector<std::vector<double> > newElts(tag.size());
    
    geometry::enhance();
    
}
*/
