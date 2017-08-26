//
//  geometry.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <Eigen/Dense>
#include <math.h>
#include <iostream>

class geometry{
    
public:
    
    // methods
    geometry(Eigen::MatrixXd &coords, Eigen::MatrixXd &elts);
    void enhance();
    void refine();
    void refine(Eigen::MatrixXd tags);
    
    // attributes
    bool enhanced;
    Eigen::MatrixXd & elements;
    Eigen::MatrixXd & coordinates;
    Eigen::MatrixXd normals;
    Eigen::MatrixXd lengths;
    Eigen::MatrixXd prev;
    Eigen::MatrixXd next;
    size_t nElts;
    
};

geometry::geometry(Eigen::MatrixXd &coords, Eigen::MatrixXd &elts)
: coordinates(coords), elements(elts)
{
    // use the constructor list above to make the needed basic arrays
    enhanced = false;
    geometry::enhance();
    
}


void geometry::enhance(){
    
    nElts = elements.rows();
    
    Eigen::MatrixXd normals = Eigen::MatrixXd::Zero(nElts, 2);
    Eigen::MatrixXd lengths = Eigen::MatrixXd::Zero(nElts, 2);
    Eigen::MatrixXd prev = Eigen::MatrixXd::Zero(nElts, 1);
    Eigen::MatrixXd next = Eigen::MatrixXd::Zero(nElts, 1);
    
    Eigen::VectorXd x1 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd x2 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd y1 = Eigen::VectorXd::Zero(nElts);
    Eigen::VectorXd y2 = Eigen::VectorXd::Zero(nElts);
    
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nElts);
    
    double norm;
    
    // fill in all the fields
    for(size_t i = 0; i<nElts; i++){
        x1(i) = coordinates(elements(i,0), 0);
        x2(i) = coordinates(elements(i,1), 0);
        y1(i) = coordinates(elements(i,0), 1);
        y2(i) = coordinates(elements(i,1), 1);
  
        lengths(i) = sqrt(pow(x2(i)-x1(i),2)+pow(y2(i)-y1(i),2));
        
        normals(i,0) = y2(i)-y1(i);
        normals(i,1) = x1(i)-x2(i);
        norm = sqrt(pow(normals(i,0),2)+pow(normals(i,1),2));
        normals(i,0) /= norm;
        normals(i,1) /= norm;

		std::cout << tmp((int) elements(i,0)) << std::endl;

        // std::cout << tmp(elements(i,0)) << std::endl;// = (int) i;
      
		// tmp(elements(i,0)) = (int) i;
		
    }
    
/*

    // it would be great to be able to do this in one felswoop
    for(size_t i = 0; i<nElts; i++){
        // next(i) = tmp(elements(i,1));
		// int n = next(i);
		// prev(n) = (int) i;        
		// prev(next(i)) = (int) i;
    }
    
    enhanced = true;

*/
    
}

/*

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

*/
