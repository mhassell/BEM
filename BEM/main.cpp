//
//  main.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include "geometry.hpp"
#include "legendrebasis.hpp"
#include "testsAndProjections.hpp"
#include "quadTables.hpp"
#include "testLegendre.hpp"

double ff(double x, double y){
    
    return pow(x,2) + pow(y,2);
    
};

int main(){
    
    // geometry mess
    std::vector<std::vector<double> > coords;
    std::vector<std::vector<int> > elts;
    coords = {{0,0},{1,0},{0.8,0.8},{0.2,1}};
    elts = {{0,1},{2,3},{1,2},{3,0}};
    geometry g(coords,elts);
    
    // function stuff
    double (*f)(double,double);
    f = ff;
    
    // quadrature
    int k = 1;
    std::vector<std::vector<double> > q1d;
    q1d = tableGauss(9);
    
    // solution array
    boost::numeric::ublas::matrix<double> fh(k+1, g.nElts);
    
    //testXh(g, f, k, q1d, fh);
    testYh(g, f, k, q1d, fh);
    
    for(int i = 0; i < k+1; i++){
        for(int j = 0; j < g.nElts; j++){
            std::cout << fh(i,j) << "     ";
        }
        std::cout << '\n' << std::endl;
    }
 
    return 0;
    
}

/*
 
 int degree = 2;
 int type = 2;
 std::vector<std::vector<double> > q1d = tableGauss(9);
 std::vector<double> pts(q1d.size());
 
 // check the points we're evaluating
 std::cout << "quadrature nodes: " << std::endl;
 for(size_t i = 0; i < q1d.size(); i++){
 pts[i] = q1d[i][0];
 std::cout << pts[i] << std::endl;
 }
 std::cout << "\n\n";
 
 // now spit out the polynomial values
 std::cout << "polynomial values: " << std::endl;
 testLegendre(degree, pts, type);
 
*/

/*
 
 // geometry mess
 std::vector<std::vector<double> > coords;
 std::vector<std::vector<int> > elts;
 coords = {{0,0},{1,0},{0.8,0.8},{0.2,1}};
 elts = {{0,1},{2,3},{1,2},{3,0}};
 geometry g(coords,elts);
 
 // function stuff
 double (*f)(double,double);
 f = ff;
 
 // quadrature
 int k = 3;
 std::vector<std::vector<double> > q1d;
 q1d = tableGauss(11);
 
 // solution array
 boost::numeric::ublas::matrix<double> fh(k+1, g.nElts);
 
 //testXh(g, f, k, q1d, fh);
 testYh(g, f, k, q1d, fh);
 
 for(int i = 0; i < k+1; i++){
 for(int j = 0; j < g.nElts; j++){
 std::cout << fh(i,j) << "     ";
 }
 std::cout << '\n' << std::endl;
 }
 
 */


/*
 
 double (*square_ptr) (double);
 square_ptr = square;
 
 double z = square_ptr(2);
 
 std::cout << z << std::endl;
 
 double square(double x){
 return pow(x,2);
 }
 
 
 */

/*
 
 size_t nPts = 5;
 int deg = 2;
 int type = 2;
 
 
 std::vector<double> *xp = new std::vector<double>(nPts);
 std::vector<double>& x = *xp;
 if(nPts == 5)
 x = {-1, -0.5, 0, 0.5, 1};
 else if (nPts == 6)
 x = {-1, -.6, -.2, .2, .6, 1};
 
 
 boost::numeric::ublas::matrix<double>* yp = new boost::numeric::ublas::matrix<double>(nPts,deg+1);
 boost::numeric::ublas::matrix<double>& y = *yp;
 
 legendreBasis(deg, x, type, y);
 
 for(int i = 0; i < nPts; i++){
 for(int j = 0; j < deg+1; j++){
 std::cout << y(i,j) << "      ";
 }
 std::cout << std::endl;
 }
 
 */

/*

 // first, test the basic constructor
 int N = 4;
 
 std::vector<std::vector<double> > *coords = new std::vector<std::vector<double> >(N);
 std::vector<std::vector<double> >& coordsref = *coords;
 std::vector<std::vector<int> >*elts = new std::vector<std::vector<int> >(N);
 std::vector<std::vector<int> >& eref = *elts;
 
 coordsref = {{0,0},{1,0},{0.8,0.8},{0.2,1}};
 eref = {{0,1},{2,3},{1,2},{3,0}};
 
 geometry g(coordsref,eref);
 
 std::cout << "Testing the normal vector \n";
 for(size_t i = 0; i<g.normals.size(); i++ ){
 for(int j = 0; j<2; j++){
 std::cout << g.normals[i][j] << ' ';
 }
 std::cout << std::endl;
 }
 
 std::cout << "Testing number of elements \n";
 std::cout << g.nElts << std::endl;
 
 
 std::cout << "Testing the lengths vector \n";
 for(size_t i = 0; i<g.lengths.size(); i++ ){
 std::cout << g.lengths[i] << ' ';
 }
 std::cout << std::endl;
 
 std::cout << "Testing the prev vector \n";
 for(size_t i = 0; i<g.prev.size(); i++ ){
 std::cout << g.prev[i] << ' ';
 }
 std::cout << std::endl;
 
 std::cout << "Testing the next vector \n";
 for(size_t i = 0; i<g.next.size(); i++ ){
 std::cout << g.next[i] << ' ';
 }
 std::cout << std::endl;
 
 std::cout << "Testing refinement: new coordinates \n";
 g.refine();
 for(size_t i = 0; i < g.nElts; i++){
 for(size_t j = 0; j < 2; j++){
 std::cout << g.coordinates[i][j] << ' ';
 }
 std::cout << std::endl;
 }
 std::cout << std::endl;
 
 std::cout << "Testing refinement: new elements \n";
 for(size_t i = 0; i < g.nElts; i++){
 for(size_t j = 0; j < 2; j++){
 std::cout << g.elements[i][j] << ' ';
 }
 std::cout << std::endl;
 }
 std::cout << std::endl;
 
 std::cout << "Number of new elements: \n";
 std::cout << g.nElts << std::endl;
 
 
 // delete the new'd stuff
 delete coords;
 delete elts;
 
*/
